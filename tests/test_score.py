"""Tests for jnuherb.score."""

import numpy as np
import pandas as pd
import pytest
from jnuherb.score import component, herbiv_score


def test_union_score():
    """The internal _union_score should compute probability union."""
    from jnuherb.score import _union_score

    s = pd.Series([0.5, 0.3])
    result = _union_score(s)
    expected = 1 - (1 - 0.5) * (1 - 0.3)  # 0.65
    assert abs(result - expected) < 1e-10

    # Empty series returns 0
    assert _union_score(pd.Series([], dtype=float)) == 0.0

    # Single value returns same value
    assert abs(_union_score(pd.Series([0.7])) - 0.7) < 1e-10


def test_herbiv_score_basic():
    """herbiv_score with synthetic data."""
    tcm = pd.DataFrame({"HVMID": ["M1", "M2"], "cn_name": ["TCM1", "TCM2"]})
    tcm_chem_links = pd.DataFrame(
        {"HVMID": ["M1", "M1", "M2"], "PubChem_id": ["C1", "C2", "C2"]}
    )
    chem = pd.DataFrame({"PubChem_id": ["C1", "C2"], "Name": ["Chem1", "Chem2"]})
    chem_prot = pd.DataFrame(
        {
            "PubChem_id": ["C1", "C2", "C1"],
            "Ensembl_ID": ["P1", "P1", "P2"],
            "Combined_score": [0.5, 0.3, 0.8],
        }
    )
    proteins = pd.DataFrame(
        {"Ensembl_ID": ["P1", "P2"], "gene_name": ["Prot1", "Prot2"]}
    )

    result = herbiv_score(tcm, tcm_chem_links, chem, chem_prot)

    tcm_scored = result["tcm_and_score"]
    chem_scored = result["chem_and_score"]

    assert "Importance Score" in tcm_scored.columns
    assert "Importance Score" in chem_scored.columns
    # Scores should be in [0, 1]
    assert (tcm_scored["Importance Score"] >= 0).all()
    assert (tcm_scored["Importance Score"] <= 1).all()

    # Sorted descending
    assert tcm_scored["Importance Score"].is_monotonic_decreasing


def test_herbiv_score_with_formula():
    """herbiv_score including formula level."""
    fml = pd.DataFrame({"HVPID": ["F1"], "name": ["Formula1"]})
    fml_tcm = pd.DataFrame({"HVPID": ["F1"], "HVMID": ["M1"]})
    tcm = pd.DataFrame({"HVMID": ["M1"], "cn_name": ["TCM1"]})
    tcm_chem_links = pd.DataFrame({"HVMID": ["M1"], "PubChem_id": ["C1"]})
    chem = pd.DataFrame({"PubChem_id": ["C1"], "Name": ["Chem1"]})
    chem_prot = pd.DataFrame(
        {"PubChem_id": ["C1"], "Ensembl_ID": ["P1"], "Combined_score": [0.7]}
    )

    result = herbiv_score(
        tcm, tcm_chem_links, chem, chem_prot,
        formula=fml, formula_tcm_links=fml_tcm,
    )

    assert result["formula_and_score"] is not None
    assert "Importance Score" in result["formula_and_score"].columns


def test_component_basic():
    """component() should produce valid combinations."""
    items = pd.DataFrame(
        {
            "HVMID": [f"M{i}" for i in range(20)],
            "cn_name": [f"TCM{i}" for i in range(20)],
            "Importance Score": np.random.RandomState(42).uniform(0.1, 0.9, 20),
        }
    )
    result = component(items, random_state=42, num=100)

    assert len(result) > 0
    assert "items" in result.columns
    assert "Importance Score" in result.columns
    assert "Boost" in result.columns
    # Sorted by Boost descending
    assert result["Boost"].is_monotonic_decreasing


def test_component_formula():
    """component() with formula IDs."""
    items = pd.DataFrame(
        {
            "HVPID": [f"F{i}" for i in range(10)],
            "name": [f"Formula{i}" for i in range(10)],
            "Importance Score": np.random.RandomState(99).uniform(0.1, 0.8, 10),
        }
    )
    result = component(items, random_state=99, num=50)
    assert len(result) > 0
    assert "Boost" in result.columns
