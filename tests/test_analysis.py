"""Tests for jnuherb.analysis."""

import pytest
from jnuherb.analysis import (
    dfs_filter,
    from_proteins,
    from_tcm_or_formula,
    suggest_formulas_from_tcms,
)
from jnuherb.query import (
    get_chem_protein_links,
    get_chemicals,
    get_formula,
    get_formula_tcm_links,
    get_proteins,
    get_tcm,
    get_tcm_chem_links,
)


def test_suggest_formulas_from_tcms():
    """The herbs of a known formula must re-suggest that formula at 100%
    coverage (e.g. the 6 herbs of 小柴胡汤 HVP1625 → 小柴胡汤 6/6)."""
    links = get_formula_tcm_links("HVPID", ["HVP1625"])
    tcm = get_tcm("HVMID", list(links["HVMID"]))
    sugg = suggest_formulas_from_tcms(tcm)
    assert len(sugg) > 0
    top = sugg.iloc[0]
    assert top["name"] == "小柴胡汤"
    assert top["matched"] == top["total"] == 6
    assert top["coverage"] == 1.0
    assert "柴胡" in top["matched_names"]
    # Sorted by matched count descending.
    assert list(sugg["matched"]) == sorted(sugg["matched"], reverse=True)


def test_suggest_formulas_from_tcms_empty():
    """No TCMs in → empty suggestions out, no exception."""
    import pandas as pd

    empty = suggest_formulas_from_tcms(pd.DataFrame())
    assert len(empty) == 0


def test_from_proteins_basic():
    """Reverse analysis with 2 proteins, no optimization."""
    result = from_proteins(
        ["ENSP00000381588", "ENSP00000252519"],
        score=0,
        tcm_component=False,
        formula_component=False,
        out_for_cytoscape=False,
        re=True,
        path="/tmp/jnuherb_test_basic",
    )
    assert result is not None
    assert "tcm" in result
    assert "chem" in result
    assert "proteins" in result
    assert len(result["tcm"]) > 0
    assert len(result["chem"]) > 0
    assert len(result["proteins"]) == 2


def test_from_proteins_with_optimization():
    """Reverse analysis with optimization."""
    result = from_proteins(
        ["ENSP00000381588", "ENSP00000252519"],
        score=0,
        random_state=138192,
        num=50,
        tcm_component=True,
        formula_component=True,
        out_for_cytoscape=False,
        re=True,
        path="/tmp/jnuherb_test_opt",
    )
    assert result is not None
    assert result["tcms"] is not None
    assert result["formulas"] is not None
    assert len(result["tcms"]) > 0
    assert "Boost" in result["tcms"].columns


def test_from_proteins_gene_name():
    """Reverse analysis using gene names."""
    result = from_proteins(
        ["PDCD1", "TP53"],
        score=0,
        tcm_component=False,
        formula_component=False,
        out_for_cytoscape=False,
        re=True,
        path="/tmp/jnuherb_test_gene",
        id_type="gene_name",
    )
    assert result is not None
    assert len(result["proteins"]) >= 1


def test_from_proteins_invalid_id_type():
    """Invalid id_type should raise ValueError."""
    with pytest.raises(ValueError, match="Invalid id_type"):
        from_proteins(
            ["PDCD1"],
            id_type="invalid_type",
            out_for_cytoscape=False,
            path="/tmp/jnuherb_test_invalid",
        )


def test_from_proteins_no_links_error():
    """Very high score should raise ValueError (no links found)."""
    with pytest.raises(ValueError, match="No compound-protein links"):
        from_proteins(
            ["ENSP00000381588"],
            score=9999,
            out_for_cytoscape=False,
            path="/tmp/jnuherb_test_empty",
        )


def test_from_tcm_or_formula_tcm():
    """Forward analysis from TCM IDs."""
    result = from_tcm_or_formula(
        ["HVM0367", "HVM1695"],
        score=900,
        out_for_cytoscape=False,
        out_graph=False,
        re=True,
        path="/tmp/jnuherb_test_fwd_tcm",
    )
    assert result is not None
    assert "tcm" in result
    assert "chem" in result
    assert "proteins" in result


def test_from_tcm_or_formula_formula():
    """Forward analysis from formula ID."""
    result = from_tcm_or_formula(
        ["HVP1625"],
        score=900,
        out_for_cytoscape=False,
        out_graph=False,
        re=True,
        path="/tmp/jnuherb_test_fwd_fml",
    )
    assert result is not None
    assert "formula" in result


def test_dfs_filter():
    """DFS filter should prune nodes without complete pathways."""
    proteins = get_proteins("Ensembl_ID", ["ENSP00000381588", "ENSP00000252519"])
    chem_prot = get_chem_protein_links("Ensembl_ID", list(proteins["Ensembl_ID"]), 0)
    chem = get_chemicals("PubChem_id", list(chem_prot["PubChem_id"]))
    tcm_chem = get_tcm_chem_links("PubChem_id", list(chem["PubChem_id"]))
    tcm = get_tcm("HVMID", list(tcm_chem["HVMID"]))
    fml_tcm = get_formula_tcm_links("HVMID", list(tcm["HVMID"]))
    fml = get_formula("HVPID", list(fml_tcm["HVPID"]))

    filtered = dfs_filter(fml, fml_tcm, tcm, tcm_chem, chem, chem_prot, proteins)

    # Filtered counts should be <= original
    assert len(filtered["proteins"]) <= len(proteins)
    assert len(filtered["tcm"]) <= len(tcm)
    assert len(filtered["chem"]) <= len(chem)

    # All filtered nodes should have complete pathways
    for _, row in filtered["proteins"].iterrows():
        assert row["Ensembl_ID"] in proteins["Ensembl_ID"].values
