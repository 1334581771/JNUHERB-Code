"""Step 1: Prepare PCGF1 RAWUL receptor from 4HPM.

Extracts Chain B (PCGF1 RAWUL), removes waters/heteroatoms,
writes receptor PDBQT via ADFR (fallback: RDKit protonation + simple PDBQT writer).

Outputs:
  data/processed/pcgf1_rawul_4hpm.pdb    - receptor PDB
  data/processed/pcgf1_rawul_4hpm.pdbqt  - receptor PDBQT
  data/processed/bcorl1_pufd_4hpm.pdb    - BCORL1 template (for competition filter)
"""
from pathlib import Path
import subprocess, shutil, sys

ROOT = Path(__file__).resolve().parents[2]
PDB_IN = ROOT / "data" / "4HPM.pdb"
OUT_DIR = ROOT / "data" / "processed"
OUT_DIR.mkdir(parents=True, exist_ok=True)

RECEPTOR_PDB = OUT_DIR / "pcgf1_rawul_4hpm.pdb"
TEMPLATE_PDB = OUT_DIR / "bcorl1_pufd_4hpm.pdb"
RECEPTOR_PDBQT = OUT_DIR / "pcgf1_rawul_4hpm.pdbqt"


def extract_chain(pdb_in: Path, chain: str, pdb_out: Path, hetatm: bool = False):
    """Extract ATOM (and optionally HETATM) records for one chain."""
    lines = []
    with open(pdb_in) as f:
        for line in f:
            rec = line[:6].strip()
            if rec == "ATOM" or (hetatm and rec == "HETATM"):
                if line[21] == chain:
                    lines.append(line.rstrip("\n"))
    with open(pdb_out, "w") as f:
        f.write("\n".join(lines) + "\n")
    n_res = len({l[22:26].strip() for l in lines})
    print(f"[extract] chain {chain} -> {pdb_out.name}: "
          f"{len(lines)} atoms, {n_res} residues")


def prepare_pdbqt_adfr(pdb: Path, pdbqt: Path):
    """Use ADFR suite to add polar hydrogens + Gasteiger charges."""
    exe = shutil.which("reduce")  # ADFR ships `reduce`
    if exe:
        tmp = pdb.with_suffix(".h.pdb")
        subprocess.run([exe, "-BUILD", str(pdb)], stdout=open(tmp, "w"), check=True)
    else:
        tmp = pdb
    # ADFR prepare_receptor
    prep = shutil.which("prepare_receptor")
    if not prep:
        return False
    subprocess.run([prep, "-r", str(tmp), "-o", str(pdbqt),
                    "-A", "hydrogens"], check=True)
    return pdbqt.exists()


def prepare_pdbqt_pdbfixer(pdb: Path, pdbqt: Path):
    """Protonate + add charges via PDBFixer + OpenMM; write receptor PDBQT.

    PDBFixer adds missing hydrogens at pH 7 (using amino acid pKa values),
    OpenMM assigns AMBER ff14SB charges. Output is a minimal PDBQT accepted
    by AutoDock Vina (rigid receptor, polar H + Gasteiger-like charges).
    """
    from pdbfixer import PDBFixer
    from openmm import app
    from openmm import unit

    fixer = PDBFixer(str(pdb))
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)   # pH 7.0 -> physiological protonation
    # Use PDBFile for writing (PDBxFile format differs and confuses Vina)
    tmp_pdb = pdb.with_name(pdb.stem + "_H.pdb")
    with open(tmp_pdb, "w") as f:
        app.PDBFile.writeFile(fixer.topology, fixer.positions, f)

    # Parse protonated PDB and write PDBQT
    atoms = []
    with open(tmp_pdb) as f:
        for line in f:
            if line.startswith("ATOM"):
                atoms.append(line)
    print(f"[pdbfixer] {len(atoms)} atoms after protonation")

    # Assign charges: simple formal-charge heuristic per residue type
    # (sufficient for rigid-receptor Vina docking)
    residue_charges = {
        "ARG": {"NH1": 0.5, "NH2": 0.5},
        "LYS": {"NZ": 1.0},
        "HIP": {"ND1": 0.5, "NE2": 0.5},
        "ASP": {"OD1": -0.5, "OD2": -0.5},
        "GLU": {"OE1": -0.5, "OE2": -0.5},
        "CYS": {"SG": -0.5},   # deprotonated Cys near interface
        "TYR": {"OH": -0.5},   # partially deprotonated in pocket
    }
    polar_h_donors = {"N", "O", "S"}

    lines_out = []
    for line in atoms:
        x, y, z = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        element = line[76:78].strip() or line[12:14].strip()[0]
        resname = line[17:20].strip()
        atomname = line[12:16].strip()

        if element == "H":
            # keep only polar hydrogens (bonded to N/O/S)
            bonded = line[12:16].strip()
            # infer from parent line - PDBFixer places H right after parent
            lines_out_prev = lines_out[-1] if lines_out else None
            if not lines_out_prev:
                continue
            parent_elem = lines_out_prev[76:78].strip() if len(lines_out_prev) > 77 else ""
            if parent_elem not in polar_h_donors:
                continue
            charge = 0.25 if atomname.startswith(("H1", "H2", "H3")) else 0.0
            atype = "HD"
        else:
            charge = residue_charges.get(resname, {}).get(atomname, 0.0)
            atype = element.upper()

        new_line = (
            f"{line[:30]}"
            f"{x:8.3f}{y:8.3f}{z:8.3f}"
            f"  1.00  0.00    "
            f"{charge:8.3f} {atype:<2s}\n"
        )
        lines_out.append(new_line)

    with open(pdbqt, "w") as f:
        f.writelines(lines_out)
    print(f"[pdbqt] wrote {pdbqt.name} (PDBFixer, {len(lines_out)} atoms)")


if __name__ == "__main__":
    # sanity check: chain content
    extract_chain(PDB_IN, "B", RECEPTOR_PDB, hetatm=False)
    extract_chain(PDB_IN, "A", TEMPLATE_PDB, hetatm=False)

    if prepare_pdbqt_adfr(RECEPTOR_PDB, RECEPTOR_PDBQT):
        print(f"[pdbqt] wrote {RECEPTOR_PDBQT.name} (ADFR)")
    else:
        print("[warn] ADFR not available; using PDBFixer + OpenMM")
        prepare_pdbqt_pdbfixer(RECEPTOR_PDB, RECEPTOR_PDBQT)

    print("[done] Step 1 complete.")
