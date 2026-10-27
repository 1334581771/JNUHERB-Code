"""Offline enrichment analysis utilities for JNUHERB.

Uses local GMT gene-set files stored in ``data/enrichment/`` and performs
hypergeometric tests via gseapy.  No external network calls are required
after the GMT files are bundled with the package.

Enrichment results are prioritised by *topic relevance* (see
:func:`run_enrichment`): instead of showing the plain top-N most
significant terms, terms are tiered — Core (matches the research topic),
Related (other non-disease mechanisms) and Review (KEGG ``Human
Diseases`` terms that usually crowd out mechanistic pathways) — and the
top-N is picked per tier.  GO results are additionally de-replicated by
gene-set overlap (Jaccard), mirroring ``enrichplot::pairwise_termsim``.
The full, untiered result table is always preserved for download.
"""

from __future__ import annotations

import re
from functools import lru_cache
from importlib.resources import files
from typing import Optional

import gseapy as gp
import numpy as np
import pandas as pd


_GMT_FILES = {
    "kegg": "kegg.gmt",
    "go_bp": "go_bp.gmt",
    "go_cc": "go_cc.gmt",
    "go_mf": "go_mf.gmt",
    # CTD / Open Targets disease associations replace the old 90-term
    # Disease Ontology GMT (14k+ disease terms, top 500 genes each)
    "do": "disease_ctdot.gmt",
}

# KEGG BRITE hierarchy (br08901): pathway name -> top-level category.
# Bundled as a TSV so no network access is needed at runtime.
_BRITE_TSV = "kegg_brite.tsv"


def _gmt_path(filename: str) -> str:
    """Return the absolute filesystem path to a bundled GMT file."""
    return str(files("jnuherb.data.enrichment").joinpath(filename))


def _read_gmt_genes(gmt_file: str) -> set[str]:
    """Read all gene symbols present in a GMT file."""
    genes = set()
    try:
        lib = gp.get_library(name=gmt_file)
        for term_genes in lib.values():
            genes.update(term_genes)
    except Exception:
        pass
    return genes


# ---------------------------------------------------------------------------
# Topic-relevance prioritisation
# ---------------------------------------------------------------------------

def _topic_tokens(topic: Optional[str]) -> list[tuple[str, bool]]:
    """Split a free-text topic (e.g. a disease name) into matchable tokens.

    Returns ``(token, is_prefix)`` pairs.  Only alphabetic tokens of >= 4
    characters are used so that short words ("of", "the", "II") cannot
    cause false-positive matches.  Long tokens (>= 7 chars) additionally
    contribute a 5-character *prefix* — matched at word boundaries only —
    so that e.g. *osteoarthritis* also surfaces *osteoclast
    differentiation* without matching mid-word fragments like the
    "arthr" inside *dysarthria*.
    """
    if not topic:
        return []
    fulls = [t for t in re.split(r"[^a-zA-Z]+", str(topic).lower())
             if len(t) >= 4]
    pairs = [(t, False) for t in fulls]
    pairs += [(t[:5], True) for t in fulls if len(t) >= 7]
    return pairs


def _topic_match(term: str, tokens: list[tuple[str, bool]]) -> bool:
    """Case-insensitive match of any topic token in the term.

    Full tokens match anywhere (so "arthritis" hits "osteoarthritis");
    prefix tokens must start a word (so "osteo" hits "osteoclast" but
    "arthr" does not hit "dysarthria").
    """
    t = str(term).lower()
    for tok, is_prefix in tokens:
        if is_prefix:
            if re.search(r"\b" + re.escape(tok), t):
                return True
        elif tok in t:
            return True
    return False


@lru_cache(maxsize=1)
def _kegg_brite_map() -> dict[str, str]:
    """Lowercase KEGG pathway name -> BRITE top-level category.

    Falls back to an empty map when the bundled TSV is missing (e.g. an
    incomplete deployment): every pathway then classifies as "Other" and
    relevance tiering degrades to significance order instead of crashing.
    """
    try:
        path = str(files("jnuherb.data.enrichment").joinpath(_BRITE_TSV))
        df = pd.read_csv(path, sep="\t")
        return dict(zip(df["term"].str.lower(), df["class1"]))
    except Exception:  # pragma: no cover - missing data file
        return {}


_TIER_RANK = {"Core": 0, "Related": 1, "Review": 2}


def _top_tiered(df: pd.DataFrame, n: int) -> pd.DataFrame:
    """Pick the n display terms from a tier-ordered result table.

    ``df`` must already be sorted by relevance tier (Core -> Related ->
    Review), then P-value within each tier — i.e. the output of
    :func:`_prioritise`.  The first ``n`` rows therefore implement the
    display rule directly: every Core term is shown first, Related terms
    fill the remainder, and Review terms (KEGG Human-Diseases entries)
    only appear when Core + Related together do not reach ``n``.
    """
    return df.head(n)


def _prioritise(
    df: pd.DataFrame,
    topic: Optional[str],
    collection: str,
) -> pd.DataFrame:
    """Add a ``Relevance`` tier column and re-sort by tier, then P-value.

    Tiers (per the enrichment-relevance-filtering workflow):
      * Core    — term matches the research topic (disease keywords).
      * Related — mechanistic term, no topic match.  For KEGG this is any
        term outside the BRITE ``Human Diseases`` branch.
      * Review  — KEGG ``Human Diseases`` terms without a topic match
        (cancer / infection pathways that share signalling genes with the
        study disease).  They are kept in the full table but pushed below
        mechanistic terms in the display top-N.

    Without a topic every term is Related and the ordering is unchanged
    (falling back to plain significance ranking).
    """
    if df.empty:
        return df
    out = df.copy()
    tokens = _topic_tokens(topic)
    if collection == "kegg":
        out["KEGG Class"] = out["Term"].str.lower().map(_kegg_brite_map())
        out["KEGG Class"] = out["KEGG Class"].fillna("Other")

        def _tier(row) -> str:
            if tokens and _topic_match(row["Term"], tokens):
                return "Core"
            if row["KEGG Class"] == "Human Diseases":
                return "Review"
            return "Related"
        out["Relevance"] = out.apply(_tier, axis=1)
    else:
        # GO and disease collections: only the topic distinguishes tiers.
        if tokens:
            out["Relevance"] = out["Term"].apply(
                lambda t: "Core" if _topic_match(t, tokens) else "Related")
        else:
            out["Relevance"] = "Related"
    out["_rank"] = out["Relevance"].map(_TIER_RANK)
    out = out.sort_values(["_rank", "P-value"]).drop(columns="_rank")
    return out.reset_index(drop=True)


def _dedup_gene_overlap(df: pd.DataFrame, cutoff: float = 0.7) -> pd.DataFrame:
    """Drop redundant GO terms by gene-set Jaccard similarity.

    Greedy: walking from most to least significant, a term is dropped when
    its gene set overlaps a kept term (same ontology category) with
    Jaccard >= ``cutoff``.  Mirrors ``enrichplot::pairwise_termsim`` with
    method = "JC" — redundant ``leukocyte migration`` /
    ``leukocyte chemotaxis`` style families collapse to their most
    significant representative.
    """
    if df.empty or "Genes" not in df.columns:
        return df

    def _gene_set(genes_str: str) -> frozenset:
        return frozenset(str(genes_str).split("/"))

    keep_mask = []
    kept: dict[str, list[frozenset]] = {}  # Category -> gene sets kept
    for _, row in df.iterrows():
        gs = _gene_set(row["Genes"])
        cat = row.get("Category", "")
        redundant = False
        for other in kept.get(cat, []):
            union = len(gs | other)
            jaccard = len(gs & other) / union if union else 0.0
            if jaccard >= cutoff:
                redundant = True
                break
        if not redundant:
            kept.setdefault(cat, []).append(gs)
        keep_mask.append(not redundant)
    return df[keep_mask].reset_index(drop=True)


def _tidy_enr_result(
    df: Optional[pd.DataFrame],
    gene_list: list[str],
    category: str,
) -> pd.DataFrame:
    """Post-process a gseapy enrich() result DataFrame.

    Sorts by raw ``P-value`` and computes Count, GeneRatio,
    ``-log10(P)`` and ``-log10(Adj P)``. Adds a ``Category`` column.
    """
    # gseapy returns a plain list (not a DataFrame) when no gene set has
    # any overlap with the input — treat that as "no enrichment" instead
    # of crashing on ``df.empty``.
    if not isinstance(df, pd.DataFrame) or df.empty:
        return pd.DataFrame()

    out = df.copy()
    out["Count"] = out["Overlap"].astype(str).str.extract(r"^(\d+)")[0].astype(int)
    n_input = len(gene_list)
    out["GeneRatio"] = out["Count"] / n_input if n_input else out["Count"]
    out["-log10(P)"] = -np.log10(out["P-value"].astype(float).apply(lambda x: max(x, 1e-300)))
    out["-log10(Adj P)"] = -np.log10(out["Adjusted P-value"].astype(float).apply(lambda x: max(x, 1e-300)))
    out["Category"] = category
    # Disease terms carry their prefixed id in the GMT name ("name
    # (MONDO:0002492)"); the id is already in the description column,
    # so strip it here for display and sentence-case the name.
    if category == "Disease" and "Term" in out.columns:
        out["Term"] = out["Term"].str.replace(
            r"\s*\([A-Za-z]+:[A-Za-z0-9]+\)$", "", regex=True)
        out["Term"] = out["Term"].str[:1].str.upper() + out["Term"].str[1:]
    # Sort by raw P-value
    out = out.sort_values("P-value")
    return out.reset_index(drop=True)


def _safe_enrich(
    gene_list: list[str],
    gmt_file: str,
    category: str,
) -> tuple[pd.DataFrame, int, int]:
    """Run offline enrichment against a local GMT file.

    Returns a tuple of (full_result_df, n_input_genes, n_mapped_genes).
    """
    n_input = len(gene_list)
    try:
        background_genes = _read_gmt_genes(gmt_file)
        n_mapped = len(set(g.upper() for g in gene_list) & set(g.upper() for g in background_genes))
        enr = gp.enrich(
            gene_list=gene_list,
            gene_sets=gmt_file,
            outdir=None,
            no_plot=True,
        )
        return _tidy_enr_result(enr.results, gene_list, category), n_input, n_mapped
    except Exception as exc:  # pragma: no cover
        empty = pd.DataFrame()
        empty.attrs["error"] = str(exc)
        return empty, n_input, 0


def run_enrichment(
    gene_list: list[str],
    which: str = "all",
    topic: Optional[str] = None,
) -> dict[str, dict]:
    """Run offline KEGG, GO (BP/CC/MF) and disease enrichment on a gene list.

    Parameters
    ----------
    gene_list:
        List of gene symbols.
    which:
        ``"all"`` computes every collection (default, backward compatible).
        ``"kegg"``, ``"go"`` or ``"disease"`` computes only that collection
        so the UI can run enrichment lazily per selected tab instead of
        computing all three on every page view.
    topic:
        Optional research topic (typically the study disease name).  When
        given, terms matching the topic are tiered ``Core`` and shown
        first; KEGG ``Human Diseases`` terms without a topic match are
        tiered ``Review`` and pushed below mechanistic terms.  The full
        result table always keeps every term for download.

    Returns
    -------
    Dictionary with keys ``kegg``, ``go``, ``do``.
    ``kegg`` and ``do`` contain ``full`` and ``top`` DataFrames plus
    ``n_input`` and ``n_mapped``.
    ``go`` contains the same keys, where ``full`` / ``top`` combine
    BP/CC/MF with an additional ``Category`` column.
    Collections not requested by ``which`` come back empty.
    """
    genes = [g.strip() for g in gene_list if g and str(g).strip()]
    empty = {"full": pd.DataFrame(), "top": pd.DataFrame(), "n_input": 0, "n_mapped": 0}
    out: dict[str, dict] = {"kegg": empty.copy(), "go": empty.copy(), "do": empty.copy()}
    if not genes:
        return out
    if which not in {"all", "kegg", "go", "disease"}:
        raise ValueError(f"Invalid which='{which}'; use all/kegg/go/disease.")

    # KEGG
    if which in ("all", "kegg"):
        kegg_full, kegg_n_in, kegg_n_map = _safe_enrich(genes, _gmt_path(_GMT_FILES["kegg"]), "KEGG")
        kegg_full = _prioritise(kegg_full, topic, "kegg")
        out["kegg"] = {"full": kegg_full, "top": _top_tiered(kegg_full, 10),
                       "n_input": kegg_n_in, "n_mapped": kegg_n_map}

    # GO: combine BP/CC/MF (de-replicated by gene overlap per ontology)
    if which in ("all", "go"):
        go_parts = []
        go_n_in, go_n_map = 0, 0
        for sub, label in [
            ("go_bp", "BP"),
            ("go_cc", "CC"),
            ("go_mf", "MF"),
        ]:
            part, n_in, n_map = _safe_enrich(genes, _gmt_path(_GMT_FILES[sub]), label)
            part = _dedup_gene_overlap(part)
            part = _prioritise(part, topic, "go")
            go_parts.append(part)
            go_n_in = max(go_n_in, n_in)
            go_n_map = max(go_n_map, n_map)

        if go_parts:
            go_full = pd.concat(go_parts, ignore_index=True)
            # Top 3 per category, relevance tier first then P-value
            go_top = pd.concat([
                go_full[go_full["Category"] == cat].head(3)
                for cat in ["BP", "CC", "MF"]
            ]).sort_values(["Category", "P-value"]).reset_index(drop=True)
        else:
            go_full = pd.DataFrame()
            go_top = pd.DataFrame()
        out["go"] = {"full": go_full, "top": go_top,
                     "n_input": go_n_in, "n_mapped": go_n_map}

    # Disease (CTD / Open Targets gene sets)
    if which in ("all", "disease"):
        do_full, do_n_in, do_n_map = _safe_enrich(genes, _gmt_path(_GMT_FILES["do"]), "Disease")
        do_full = _prioritise(do_full, topic, "disease")
        out["do"] = {"full": do_full, "top": _top_tiered(do_full, 10),
                     "n_input": do_n_in, "n_mapped": do_n_map}

    return out


def extract_genes_from_result(result: Optional[dict]) -> list[str]:
    """Extract unique gene symbols from a JNUHERB analysis result."""
    if result is None:
        return []
    prot_df = result.get("protein_result", {}).get("proteins")
    if prot_df is None or prot_df.empty:
        return []

    if "gene_name" in prot_df.columns:
        genes = prot_df["gene_name"].dropna().astype(str).tolist()
    elif "Gene" in prot_df.columns:
        genes = prot_df["Gene"].dropna().astype(str).tolist()
    else:
        genes = []
        for col in prot_df.columns:
            if "gene" in col.lower() or "name" in col.lower():
                genes = prot_df[col].dropna().astype(str).tolist()
                break

    seen = set()
    out = []
    for g in genes:
        if g and g not in seen:
            seen.add(g)
            out.append(g)
    return out
