"""Core network pharmacology analysis functions.

This module provides the three central analysis workflows:

* ``from_tcm_or_formula`` — **forward** analysis: TCM/formula → chemicals → proteins.
* ``from_proteins`` — **reverse** analysis: proteins → chemicals → TCMs → formulas.
* ``resolve_disease_to_targets`` — disease lookup that seeds Reverse mode,
  with LLM translation for Chinese disease names (the link table only
  carries English MONDO/EFO/MeSH/Orphanet vocabulary).
* ``dfs_filter`` — depth-first search filter that keeps only nodes with
  complete pathways through the network.
"""

from __future__ import annotations

import re
from typing import Any, Optional

import pandas as pd

from jnuherb.query import (
    get_chem_protein_links,
    get_chemicals,
    get_formula,
    get_formula_tcm_links,
    get_protein_disease_links,
    get_proteins,
    get_tcm,
    get_tcm_chem_links,
)

_CJK_RE = re.compile(r"[一-鿿]")


def resolve_disease_to_targets(
    disease: str,
    min_score: float = 0.1,
    max_targets: int = 200,
    sources: Optional[list[str]] = None,
    translate: bool = True,
    api_key: Optional[str] = None,
    disease_id: Optional[str] = None,
) -> dict[str, Any]:
    """Resolve a disease name/id to seed target proteins for Reverse mode.

    Looks the disease up in the protein–disease link table (CTD / Open
    Targets), ranks associated proteins by evidence score and returns the
    top Ensembl IDs, ready to feed :func:`from_proteins`.

    The link table only carries English vocabulary (MONDO/EFO/Orphanet
    from OT, MeSH from CTD).  When a query contains Chinese characters
    and matches nothing directly, it is translated through the LLM
    (``jnuherb.chat.translate_disease_name``) and the English term is
    looked up instead.  Translation is best-effort: any API failure just
    yields the normal "unknown disease" result.

    Args:
        disease: Disease name (English, case-insensitive) or prefixed id
            (``EFO:0000180``, ``MONDO:0008903``, ``MESH:D006973``), or a
            Chinese disease name when ``translate`` is enabled.
        min_score: Evidence score floor; CTD curated rows (empty score)
            always pass.
        max_targets: Cap on the returned seed target list.
        sources: Restrict to link sources (e.g. ``["OT", "CTD"]``).
        translate: Allow the LLM Chinese→English fallback.
        api_key: DeepSeek API key for the translation fallback (falls back
            to env / Streamlit secrets when None).
        disease_id: When given (disease disambiguation), restrict the
            lookup to this exact vocabulary ID instead of letting a
            name match span multiple diseases.

    Returns:
        ``{"disease_id", "disease_name", "targets", "n_total",
        "min_score", "top_genes", "translated_from"}`` — ``targets`` is
        empty when the disease is unknown; ``translated_from`` is the
        original Chinese query when a translation was used, else None.
    """
    # With a disambiguated pick, query by the exact vocabulary ID — the
    # by_disease name filter only matches whole names, so a short name
    # like "stroke" would never reach its "stroke disorder" entry.
    lookup = disease_id if disease_id is not None else disease
    links = get_protein_disease_links(by_disease=lookup, min_score=min_score,
                                      sources=sources)
    if disease_id is not None:
        # Belt-and-braces: keep only the chosen vocabulary entry so a
        # name match never mixes targets from several same-named diseases.
        links = links[links["disease_id"] == disease_id]
    translated_from: Optional[str] = None
    if links.empty and translate and _CJK_RE.search(disease):
        # Lazy import: keeps chat.py (and its LLM deps) out of the
        # import graph unless a Chinese query actually needs it.
        from jnuherb.chat import translate_disease_name

        en = translate_disease_name(disease, api_key=api_key)
        if en:
            tlinks = get_protein_disease_links(
                by_disease=en, min_score=min_score, sources=sources)
            if disease_id is not None:
                tlinks = tlinks[tlinks["disease_id"] == disease_id]
            if not tlinks.empty:
                links = tlinks
                translated_from = disease

    out: dict[str, Any] = {
        "disease_id": "", "disease_name": disease,
        "targets": [], "n_total": 0, "min_score": min_score,
        "top_genes": [], "translated_from": translated_from,
    }
    if links.empty:
        return out

    best = (links.groupby(["Ensembl_ID", "gene_name", "disease_id",
                           "disease_name"], as_index=False)
                 .agg(score=("score", "max")))
    best = best.sort_values("score", ascending=False, na_position="last")
    out["disease_id"] = best["disease_id"].iloc[0]
    out["disease_name"] = best["disease_name"].iloc[0] or disease
    out["n_total"] = len(best)
    out["targets"] = best["Ensembl_ID"].head(max_targets).tolist()
    out["top_genes"] = best["gene_name"].head(max_targets).tolist()
    return out


# ---------------------------------------------------------------------------
# dfs_filter
# ---------------------------------------------------------------------------


def dfs_filter(
    formula: Optional[pd.DataFrame],
    formula_tcm_links: Optional[pd.DataFrame],
    tcm: pd.DataFrame,
    tcm_chem_links: pd.DataFrame,
    chem: pd.DataFrame,
    chem_protein_links: pd.DataFrame,
    proteins: pd.DataFrame,
) -> dict[str, Optional[pd.DataFrame]]:
    """Filter nodes to keep only complete pathways.

    Keeps only nodes that participate in at least one complete path
    ending at a valid target protein. Implemented as a backward pass
    from proteins → chemicals → TCMs → formulas using set operations,
    which is much faster than the original nested-loop forward traversal.

    Args:
        formula: Formula information (or None if not in formula mode).
        formula_tcm_links: Formula-TCM links (or None).
        tcm: TCM information.
        tcm_chem_links: TCM-chemical links.
        chem: Chemical information.
        chem_protein_links: Chemical-protein links.
        proteins: Protein information.

    Returns:
        A dict with keys ``formula``, ``formula_tcm_links``, ``tcm``,
        ``tcm_chem_links``, ``chem``, ``chem_protein_links``, ``proteins``
        containing filtered DataFrames (or None where appropriate).
    """
    target_protein_ids = set(proteins["Ensembl_ID"])

    # Backward pass: start from target proteins
    valid_chem_protein_links = chem_protein_links[
        chem_protein_links["Ensembl_ID"].isin(target_protein_ids)
    ]
    chem_id_set = set(valid_chem_protein_links["HVCID"])

    valid_tcm_chem_links = tcm_chem_links[
        tcm_chem_links["HVCID"].isin(chem_id_set)
    ]
    tcm_id_set = set(valid_tcm_chem_links["HVMID"])

    formula_id_set: Optional[set[str]] = None
    if formula_tcm_links is not None and formula is not None:
        valid_formula_tcm_links = formula_tcm_links[
            formula_tcm_links["HVMID"].isin(tcm_id_set)
        ]
        formula_id_set = set(valid_formula_tcm_links["HVPID"])

    # Filter each DataFrame
    if formula is not None and formula_id_set is not None:
        formula = formula[formula["HVPID"].isin(formula_id_set)].copy()
    if formula_tcm_links is not None:
        if formula_id_set is not None:
            formula_tcm_links = formula_tcm_links[
                formula_tcm_links["HVPID"].isin(formula_id_set)
                & formula_tcm_links["HVMID"].isin(tcm_id_set)
            ].copy()
        else:
            formula_tcm_links = formula_tcm_links[
                formula_tcm_links["HVMID"].isin(tcm_id_set)
            ].copy()

    tcm = tcm[tcm["HVMID"].isin(tcm_id_set)].copy()
    chem = chem[chem["HVCID"].isin(chem_id_set)].copy()

    # Proteins: keep only those actually connected to surviving chemicals
    connected_protein_ids = set(valid_chem_protein_links["Ensembl_ID"])
    proteins = proteins[proteins["Ensembl_ID"].isin(connected_protein_ids)].copy()

    tcm_chem_links = tcm_chem_links[
        tcm_chem_links["HVMID"].isin(tcm_id_set)
        & tcm_chem_links["HVCID"].isin(chem_id_set)
    ].copy()
    chem_protein_links = chem_protein_links[
        chem_protein_links["HVCID"].isin(chem_id_set)
        & chem_protein_links["Ensembl_ID"].isin(connected_protein_ids)
    ].copy()

    # Reset indices
    for df in [tcm_chem_links, chem_protein_links, proteins]:
        df.reset_index(drop=True, inplace=True)

    return {
        "formula": formula,
        "formula_tcm_links": formula_tcm_links,
        "tcm": tcm,
        "tcm_chem_links": tcm_chem_links,
        "chem": chem,
        "chem_protein_links": chem_protein_links,
        "proteins": proteins,
    }


# ---------------------------------------------------------------------------
# from_tcm_or_formula  (forward analysis)
# ---------------------------------------------------------------------------


def from_tcm_or_formula(
    tcm_or_formula_id: list[str],
    proteins_id: Optional[list[str]] = None,
    score: float = 990,
    out_for_cytoscape: bool = True,
    out_graph: bool = True,
    re: bool = True,
    path: str = "results",
    id_type: str = "HVMID",
) -> Optional[dict[str, Any]]:
    """Forward network pharmacology analysis from TCMs or formulas.

    Traces the pathway **TCM/formula → chemicals → protein targets**.

    Args:
        tcm_or_formula_id: TCM IDs (HVMID) or formula IDs (HVPID).
            Also supports ``cn_name``, ``pinyin_name``, ``en_name``,
            ``latin_name`` when ``id_type`` is set accordingly.
        proteins_id: Optional protein Ensembl IDs to restrict to.
            If None, all proteins connected to the chemicals are used.
        score: Minimum ``Combined_score`` threshold for chemical-protein
            links (default 990, i.e. 0.99 after conversion).
        out_for_cytoscape: Whether to write Cytoscape-compatible CSV files.
        out_graph: Whether to generate an interactive HTML network graph.
        re: Whether to return the analysis result.
        path: Directory path for output files.
        id_type: Column name used to match ``tcm_or_formula_id``.
            One of ``"HVMID"``, ``"cn_name"``, ``"pinyin_name"``,
            ``"en_name"``, ``"latin_name"``.  For formulas only
            ``"HVPID"`` is used automatically.

    Returns:
        If ``re`` is True, a dict with keys ``formula``,
        ``formula_tcm_links``, ``tcm``, ``tcm_chem_links``, ``chem``,
        ``chem_protein_links``, ``proteins``.  The ``formula`` keys
        are present only when the input is a formula (HVPID).
        Returns None if ``re`` is False.
    """
    # --- validate id_type ---
    valid_id_types = {"HVMID", "cn_name", "pinyin_name", "en_name", "latin_name"}
    if id_type not in valid_id_types:
        raise ValueError(
            f"Invalid id_type '{id_type}'. Must be one of: "
            f"{', '.join(sorted(valid_id_types))}"
        )

    # --- detect formula vs TCM ---
    if len(tcm_or_formula_id) == 0:
        raise ValueError("tcm_or_formula_id must not be empty.")
    first_id = tcm_or_formula_id[0]
    is_formula = first_id.startswith("HVP")

    if is_formula:
        fml = get_formula("HVPID", tcm_or_formula_id)
        fml_tcm_links = get_formula_tcm_links("HVPID", list(fml["HVPID"]))
        tcm_df = get_tcm("HVMID", list(fml_tcm_links["HVMID"]))
    else:
        fml = None
        fml_tcm_links = None
        tcm_df = get_tcm(id_type, tcm_or_formula_id)

    tcm_chem_links = get_tcm_chem_links("HVMID", list(tcm_df["HVMID"]))
    chem = get_chemicals("HVCID", list(tcm_chem_links["HVCID"]))
    chem_prot_links = get_chem_protein_links("HVCID", list(chem["HVCID"]), score)

    if proteins_id is None:
        proteins = get_proteins("Ensembl_ID", list(chem_prot_links["Ensembl_ID"]))
    else:
        proteins = get_proteins("Ensembl_ID", proteins_id)

    # --- DFS filter ---
    filtered = dfs_filter(
        fml, fml_tcm_links, tcm_df, tcm_chem_links,
        chem, chem_prot_links, proteins,
    )
    fml = filtered["formula"]
    fml_tcm_links = filtered["formula_tcm_links"]
    tcm_df = filtered["tcm"]
    tcm_chem_links = filtered["tcm_chem_links"]
    chem = filtered["chem"]
    chem_prot_links = filtered["chem_protein_links"]
    proteins = filtered["proteins"]

    # --- scores ---
    from jnuherb.score import herbiv_score  # local to avoid circular import

    score_result = herbiv_score(
        tcm_df, tcm_chem_links, chem, chem_prot_links,
        formula=fml, formula_tcm_links=fml_tcm_links,
    )
    tcm_df = score_result["tcm_and_score"]
    chem = score_result["chem_and_score"]
    fml = score_result["formula_and_score"]

    # --- output ---
    if out_for_cytoscape:
        from jnuherb.output import out_for_cyto

        out_for_cyto(tcm_df, tcm_chem_links, chem, chem_prot_links, proteins, path)

    if out_graph:
        from jnuherb.output import vis_network

        vis_network(tcm_df, tcm_chem_links, chem, chem_prot_links, proteins, path)

    if re:
        if fml is not None:
            return {
                "formula": fml,
                "formula_tcm_links": fml_tcm_links,
                "tcm": tcm_df,
                "tcm_chem_links": tcm_chem_links,
                "chem": chem,
                "chem_protein_links": chem_prot_links,
                "proteins": proteins,
            }
        else:
            return {
                "tcm": tcm_df,
                "tcm_chem_links": tcm_chem_links,
                "chem": chem,
                "chem_protein_links": chem_prot_links,
                "proteins": proteins,
            }

    return None


# ---------------------------------------------------------------------------
# from_proteins  (reverse analysis)
# ---------------------------------------------------------------------------


def from_proteins(
    proteins: list[str],
    score: float = 0,
    random_state: Optional[int] = None,
    num: int = 1000,
    tcm_component: bool = True,
    formula_component: bool = True,
    out_for_cytoscape: bool = True,
    re: bool = True,
    path: str = "result",
    id_type: str = "Ensembl_ID",
) -> Optional[dict[str, Any]]:
    """Reverse network pharmacology analysis from target proteins.

    Traces the pathway **proteins → chemicals → TCMs → formulas**
    and optionally optimizes TCM/formula combinations.

    Args:
        proteins: Protein IDs (Ensembl IDs or gene names depending on
            ``id_type``).
        score: Minimum ``Combined_score`` threshold for chemical-protein
            links (default 0, meaning all links are included).
        random_state: Random seed for reproducible optimization.
        num: Number of solution sets in optimization (default 1000).
        tcm_component: Whether to run TCM combination optimization.
        formula_component: Whether to run formula combination optimization.
        out_for_cytoscape: Whether to write Cytoscape CSV files.
        re: Whether to return the analysis result.
        path: Directory path for output files.
        id_type: ``"Ensembl_ID"`` (default) or ``"gene_name"``.

    Returns:
        If ``re`` is True, a dict with keys ``formula``,
        ``formula_tcm_links``, ``tcm``, ``tcm_chem_links``, ``chem``,
        ``chem_protein_links``, ``proteins``, ``tcms`` (optimized
        combinations), ``formulas`` (optimized combinations).
    """
    # --- validate id_type ---
    valid_id_types = {"Ensembl_ID", "gene_name"}
    if id_type not in valid_id_types:
        raise ValueError(
            f"Invalid id_type '{id_type}'. Must be one of: "
            f"{', '.join(sorted(valid_id_types))}"
        )

    # --- trace: proteins → chemicals → TCMs → formulas ---
    prot_df = get_proteins(id_type, proteins)

    if len(prot_df) == 0:
        # Build helpful suggestion
        if id_type == "Ensembl_ID":
            hint = (
                "These Ensembl ID(s) were not found in the database. "
                "The database contains ~14,200 proteins with known "
                "TCM-related chemical interactions — it is NOT the full "
                "human proteome. Try using id_type='gene_name' with gene "
                "symbols instead (e.g. PDCD1, TP53, TNF, IL6), or verify "
                "your Ensembl IDs at https://useast.ensembl.org."
            )
        else:
            hint = (
                "These gene name(s) were not found. Try using "
                "id_type='Ensembl_ID' with Ensembl protein IDs instead "
                "(e.g. ENSP00000381588 for TMPRSS2)."
            )
        raise ValueError(
            f"No proteins found for the given IDs using id_type='{id_type}'. "
            f"{hint}"
        )

    chem_prot_links = get_chem_protein_links(
        "Ensembl_ID", list(prot_df["Ensembl_ID"]), score
    )

    if len(chem_prot_links) == 0:
        raise ValueError(
            f"No compound-protein links found for the given proteins "
            f"at score >= {score}. The protein(s) exist in the database "
            f"but have no known chemical interactions. "
            f"Try different proteins or check if the targets are in "
            f"the chemical-protein interaction database."
        )

    chem = get_chemicals("HVCID", list(chem_prot_links["HVCID"]))
    tcm_chem_links = get_tcm_chem_links("HVCID", list(chem["HVCID"]))
    tcm_df = get_tcm("HVMID", list(tcm_chem_links["HVMID"]))
    fml_tcm_links = get_formula_tcm_links("HVMID", list(tcm_df["HVMID"]))
    fml = get_formula("HVPID", list(fml_tcm_links["HVPID"]))

    # --- DFS filter ---
    filtered = dfs_filter(
        fml, fml_tcm_links, tcm_df, tcm_chem_links,
        chem, chem_prot_links, prot_df,
    )
    fml = filtered["formula"]
    fml_tcm_links = filtered["formula_tcm_links"]
    tcm_df = filtered["tcm"]
    tcm_chem_links = filtered["tcm_chem_links"]
    chem = filtered["chem"]
    chem_prot_links = filtered["chem_protein_links"]
    prot_df = filtered["proteins"]

    # --- scores ---
    from jnuherb.score import component, herbiv_score

    score_result = herbiv_score(
        tcm_df, tcm_chem_links, chem, chem_prot_links,
        formula=fml, formula_tcm_links=fml_tcm_links,
    )
    tcm_df = score_result["tcm_and_score"]
    chem = score_result["chem_and_score"]
    fml = score_result["formula_and_score"]

    # --- component optimization ---
    tcms_opt: Optional[pd.DataFrame] = None
    formulas_opt: Optional[pd.DataFrame] = None

    if tcm_component:
        tcm_subset = tcm_df[tcm_df["Importance Score"] != 1.0]
        if len(tcm_subset) > 0:
            tcms_opt = component(tcm_subset, random_state, num)

    if formula_component:
        fml_subset = fml[fml["Importance Score"] != 1.0]
        if len(fml_subset) > 0:
            formulas_opt = component(fml_subset, random_state, num)

    # --- output ---
    if out_for_cytoscape:
        from jnuherb.output import out_for_cyto

        out_for_cyto(tcm_df, tcm_chem_links, chem, chem_prot_links, prot_df, path)

    if re:
        return {
            "formula": fml,
            "formula_tcm_links": fml_tcm_links,
            "tcm": tcm_df,
            "tcm_chem_links": tcm_chem_links,
            "chem": chem,
            "chem_protein_links": chem_prot_links,
            "proteins": prot_df,
            "tcms": tcms_opt,
            "formulas": formulas_opt,
        }

    return None
