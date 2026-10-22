"""Offline enrichment analysis utilities for JNUHERB.

Uses local GMT gene-set files stored in ``data/enrichment/`` and performs
hypergeometric tests via gseapy.  No external network calls are required
after the GMT files are bundled with the package.
"""

from __future__ import annotations

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
        out["kegg"] = {"full": kegg_full, "top": kegg_full.head(10),
                       "n_input": kegg_n_in, "n_mapped": kegg_n_map}

    # GO: combine BP/CC/MF
    if which in ("all", "go"):
        go_parts = []
        go_n_in, go_n_map = 0, 0
        for sub, label in [
            ("go_bp", "BP"),
            ("go_cc", "CC"),
            ("go_mf", "MF"),
        ]:
            part, n_in, n_map = _safe_enrich(genes, _gmt_path(_GMT_FILES[sub]), label)
            go_parts.append(part)
            go_n_in = max(go_n_in, n_in)
            go_n_map = max(go_n_map, n_map)

        if go_parts:
            go_full = pd.concat(go_parts, ignore_index=True).sort_values("P-value").reset_index(drop=True)
            # Top 3 per category, preserving overall P-value order within each
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
        out["do"] = {"full": do_full, "top": do_full.head(10),
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
