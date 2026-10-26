"""Tests for jnuherb.query — mirrors R tests/test-get.R."""

import pytest
from jnuherb.query import (
    get_chem_protein_links,
    get_chemicals,
    get_formula,
    get_formula_tcm_links,
    get_proteins,
    get_tcm,
    get_tcm_chem_links,
)


def test_get_formula():
    df = get_formula("HVPID", ["HVP1625"])
    assert len(df) >= 1
    assert "HVPID" in df.columns


def test_get_tcm():
    df = get_tcm("HVMID", ["HVM0367", "HVM1695"])
    assert len(df) == 2
    assert "HVMID" in df.columns


def test_get_chemicals():
    df = get_chemicals("PubChem_id", ["1032"])  # acetic acid
    assert len(df) == 1
    assert "PubChem_id" in df.columns
    assert "HVCID" in df.columns  # provenance column kept
    # Should be deduplicated on PubChem_id
    assert df["PubChem_id"].is_unique


def test_get_proteins():
    df = get_proteins("Ensembl_ID", ["ENSP00000381588", "ENSP00000252519"])
    assert len(df) == 2
    assert "Ensembl_ID" in df.columns
    assert df["Ensembl_ID"].is_unique


def test_get_formula_tcm_links():
    df = get_formula_tcm_links("HVPID", ["HVP1625"])
    assert len(df) >= 1


def test_get_tcm_chem_links():
    df = get_tcm_chem_links("HVMID", ["HVM0367"])
    assert len(df) >= 1


def test_get_chem_protein_links():
    df = get_chem_protein_links("Ensembl_ID", ["ENSP00000381588"], score=0)
    assert df.shape[1] >= 3
    # Combined_score should be in [0, 1]
    assert (df["Combined_score"].astype(float) >= 0).all()
    assert (df["Combined_score"].astype(float) <= 1).all()
