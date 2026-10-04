"""Unit tests for the docking module's pure-Python helpers.

These tests intentionally avoid rdkit/meeko/vina so they run anywhere.
"""
from pathlib import Path

import pytest

from jnuherb.docking import (
    DockingError,
    binding_box_from_residues,
    default_box_from_receptor,
    make_complex_pdb,
    parse_residue_list,
    pose_pdbqt_to_pdb,
)


# ---------------------------------------------------------------------------
# Synthetic receptor: residues 10 and 20 with CA atoms at known positions
# ---------------------------------------------------------------------------

def _atom_line(serial, name, resname, resid, x, y, z, elem):
    return (
        f"ATOM  {serial:5d} {name:<4s}{resname:>3s} A{resid:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {elem:>2s}"
    )


@pytest.fixture()
def receptor_pdb(tmp_path: Path) -> Path:
    lines = [
        "HEADER    SYNTHETIC TEST RECEPTOR",
        _atom_line(1, "N", "ALA", 10, 0.0, 0.0, 0.0, "N"),
        _atom_line(2, "CA", "ALA", 10, 1.0, 0.0, 0.0, "C"),
        _atom_line(3, "C", "ALA", 10, 2.0, 0.0, 0.0, "C"),
        _atom_line(4, "O", "ALA", 10, 3.0, 0.0, 0.0, "O"),
        _atom_line(5, "CB", "ALA", 10, 1.0, 1.0, 0.0, "C"),
        _atom_line(6, "CA", "GLY", 20, 5.0, 5.0, 5.0, "C"),
        "END",
    ]
    p = tmp_path / "rec.pdb"
    p.write_text("\n".join(lines) + "\n")
    return p


# ---------------------------------------------------------------------------
# parse_residue_list
# ---------------------------------------------------------------------------

def test_parse_residue_list_basic():
    assert parse_residue_list("163, 168, 194") == [163, 168, 194]


def test_parse_residue_list_range_and_dedup():
    assert parse_residue_list("163-165, 164, 168") == [163, 164, 165, 168]


def test_parse_residue_list_chinese_comma():
    assert parse_residue_list("163，168") == [163, 168]


def test_parse_residue_list_invalid():
    with pytest.raises(DockingError):
        parse_residue_list("abc")


# ---------------------------------------------------------------------------
# binding_box_from_residues
# ---------------------------------------------------------------------------

def test_box_sidechain_centroid(receptor_pdb):
    # residue 10: heavy non-backbone atoms = CB at (1,1,0) -> center (1,1,0)
    center, size = binding_box_from_residues(receptor_pdb, [10])
    assert center == (1.0, 1.0, 0.0)
    assert size == (20.0, 20.0, 20.0)


def test_box_glycine_falls_back_to_ca(receptor_pdb):
    # GLY 20 has no side chain -> CA at (5,5,5)
    center, _ = binding_box_from_residues(receptor_pdb, [20])
    assert center == (5.0, 5.0, 5.0)


def test_box_multiple_residues(receptor_pdb):
    # CB of ALA10 (1,1,0) and CA of GLY20 (5,5,5) -> midpoint (3,3,2.5)
    center, _ = binding_box_from_residues(receptor_pdb, [10, 20])
    assert center == (3.0, 3.0, 2.5)


def test_box_missing_residue_raises(receptor_pdb):
    with pytest.raises(DockingError, match="were found"):
        binding_box_from_residues(receptor_pdb, [999])


# ---------------------------------------------------------------------------
# default_box_from_receptor (blind-docking fallback)
# ---------------------------------------------------------------------------

def test_default_box_ca_centroid(receptor_pdb):
    # C-alpha atoms at (1,0,0) and (5,5,5) -> centroid (3.0, 2.5, 2.5)
    center, size = default_box_from_receptor(receptor_pdb)
    assert center == (3.0, 2.5, 2.5)
    assert size == (20.0, 20.0, 20.0)


# ---------------------------------------------------------------------------
# pose_pdbqt_to_pdb / make_complex_pdb
# ---------------------------------------------------------------------------

@pytest.fixture()
def pose_pdbqt(tmp_path: Path) -> Path:
    # Minimal Vina-style pose: 2 atoms, one lacking an element column
    lines = [
        "REMARK VINA RESULT: -7.5 0.0 0.0",
        "ROOT",
        "ATOM      1  C1  LIG A   1      10.000  11.000  12.000  0.00  0.00    +0.00 C ",
        "ATOM      2  O1  LIG A   1      11.000  12.000  13.000  0.00  0.00    -0.00 OA",
        "ENDROOT",
        "TORSDOF 0",
    ]
    p = tmp_path / "pose.pdbqt"
    p.write_text("\n".join(lines) + "\n")
    return p


def test_pose_to_pdb(pose_pdbqt):
    pdb = pose_pdbqt_to_pdb(pose_pdbqt, "aspirin")
    hetatms = [l for l in pdb.splitlines() if l.startswith("HETATM")]
    assert len(hetatms) == 2
    assert hetatms[0][30:54].split() == ["10.000", "11.000", "12.000"]
    assert "LIG" in hetatms[0][17:20]
    assert pdb.rstrip().endswith("TER")


def test_pose_to_pdb_infers_element(pose_pdbqt):
    pdb = pose_pdbqt_to_pdb(pose_pdbqt, "LIG")
    hetatms = [l for l in pdb.splitlines() if l.startswith("HETATM")]
    # first line has explicit "C" element; second has "OA" -> falls back to O
    assert hetatms[0][-2:].strip() == "C"
    assert hetatms[1][-2:].strip() == "O"


def test_make_complex(receptor_pdb, pose_pdbqt):
    complex_pdb = make_complex_pdb(receptor_pdb, pose_pdbqt,
                                   ligand_name="aspirin", affinity=-7.5)
    lines = complex_pdb.splitlines()
    assert any("VINA_AFFINITY_KCAL_MOL: -7.50" in l for l in lines)
    assert any("RECEPTOR: rec.pdb" in l for l in lines)
    n_atom = sum(1 for l in lines if l.startswith("ATOM"))
    n_het = sum(1 for l in lines if l.startswith("HETATM"))
    assert n_atom == 6   # receptor kept intact
    assert n_het == 2    # ligand added
    assert lines[-1] == "END"
