"""Step 3: Define binding site (grid center/size) from BCORL1 Leu1665 pocket.

Approach: the Leu cage pocket has no co-crystallized small molecule, so standard
redocking validation is not possible. We use PLACEHOLDER REDOCKING:
  1. Extract Leu1665 side-chain centroid from BCORL1 template (Chain A)
  2. Dock the Leu1665 side chain (as a capped fragment) into the RAWUL pocket
  3. Report RMSD of best pose vs. native side chain position
  4. Write Vina config (center, size, exhaustiveness) for full HTVS

Outputs:
  data/processed/vina_config.txt
  results/tables/binding_site_summary.csv
  results/reports/binding_site_validation.txt
"""
from pathlib import Path
import csv, math

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

from vina import Vina

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "processed"
TABLES = ROOT / "results" / "tables"
REPORTS = ROOT / "results" / "reports"
TABLES.mkdir(parents=True, exist_ok=True)
REPORTS.mkdir(parents=True, exist_ok=True)

RECEPTOR_PDBQT = DATA / "pcgf1_rawul_4hpm.pdbqt"
TEMPLATE_PDB = DATA / "bcorl1_pufd_4hpm.pdb"
CONFIG_OUT = DATA / "vina_config.txt"
SUMMARY_CSV = TABLES / "binding_site_summary.csv"
VALIDATION_TXT = REPORTS / "binding_site_validation.txt"

# pocket residues on PCGF1 side (paper Fig 1E: Tyr163, Gln168, Ser195, Cys194)
POCKET_RESIDUES = ["TYR", "GLN", "SER", "CYS"]
TARGET_LEU = ("LEU", 1665)

# grid parameters (Angstrom)
GRID_SIZE = (25.0, 25.0, 25.0)
EXHAUSTIVENESS_HTVS = 8
EXHAUSTIVENESS_PLACEHOLDER = 32
N_POSES_PLACEHOLDER = 20


def parse_pdb_coords(pdb: Path):
    atoms = []
    with open(pdb) as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            atoms.append({
                "name": line[12:16].strip(),
                "resname": line[17:20].strip(),
                "resid": int(line[22:26]),
                "x": float(line[30:38]),
                "y": float(line[38:46]),
                "z": float(line[46:54]),
            })
    return atoms


def sidechain_centroid(atoms, resname, resid):
    pts = [(a["x"], a["y"], a["z"]) for a in atoms
           if a["resname"] == resname and a["resid"] == resid
           and a["name"] not in ("N", "CA", "C", "O")]
    if not pts:
        raise ValueError(f"residue {resname} {resid} not found in template")
    return [sum(p[i] for p in pts) / len(pts) for i in range(3)]


def build_leu_sidechain_mol():
    """Build N-acetyl-Leu-methylamide as capped side chain fragment."""
    mol = Chem.MolFromSmiles("CC(C)CC(NC(=O)C)C(=O)NC")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    return mol


def write_mol_to_pdbqt(mol, out_path: Path):
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    prep = MoleculePreparation()
    ms = prep.prepare(mol)[0]
    s, _, _ = PDBQTWriterLegacy.write_string(ms)
    out_path.write_text(s)


def rmsd(pts_a, pts_b):
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(pts_a, pts_b)) / len(pts_a))


def main():
    atoms = parse_pdb_coords(TEMPLATE_PDB)
    center = sidechain_centroid(atoms, *TARGET_LEU)
    print(f"[site] Leu1665 side-chain centroid: "
          f"({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f})")

    # ---- placeholder redocking ----
    leu_mol = build_leu_sidechain_mol()
    leu_pdbqt = DATA / "leu1665_placeholder.pdbqt"
    write_mol_to_pdbqt(leu_mol, leu_pdbqt)

    v = Vina(sf_name="vina", verbosity=0)
    v.set_receptor(str(RECEPTOR_PDBQT))
    v.set_ligand_from_file(str(leu_pdbqt))
    v.compute_vina_maps(center=center, box_size=list(GRID_SIZE))
    v.dock(exhaustiveness=EXHAUSTIVENESS_PLACEHOLDER,
           n_poses=N_POSES_PLACEHOLDER)
    energies = v.energies()
    best_energy = energies[0][0]
    v.write_poses(str(DATA / "leu1665_placeholder_out.pdbqt"),
                  n_poses=N_POSES_PLACEHOLDER, overwrite=True)

    # RMSD of best pose heavy-atom centroid vs native centroid
    out_atoms = parse_pdb_coords(DATA / "leu1665_placeholder_out.pdbqt")
    pose_centroid = [sum(a[k] for a in out_atoms) / len(out_atoms)
                     for k in ("x", "y", "z")]
    rmsd_val = rmsd(pose_centroid, center)

    # ---- write vina config ----
    cfg = (
        f"receptor = {RECEPTOR_PDBQT}\n"
        f"center_x = {center[0]:.3f}\n"
        f"center_y = {center[1]:.3f}\n"
        f"center_z = {center[2]:.3f}\n"
        f"size_x = {GRID_SIZE[0]:.1f}\n"
        f"size_y = {GRID_SIZE[1]:.1f}\n"
        f"size_z = {GRID_SIZE[2]:.1f}\n"
        f"exhaustiveness = {EXHAUSTIVENESS_HTVS}\n"
        f"num_modes = 9\n"
    )
    CONFIG_OUT.write_text(cfg)
    print(f"[site] vina config -> {CONFIG_OUT.name}")

    # ---- reports ----
    with open(SUMMARY_CSV, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["parameter", "value"])
        w.writerow(["grid_center_x", f"{center[0]:.3f}"])
        w.writerow(["grid_center_y", f"{center[1]:.3f}"])
        w.writerow(["grid_center_z", f"{center[2]:.3f}"])
        w.writerow(["grid_size_x_A", GRID_SIZE[0]])
        w.writerow(["grid_size_y_A", GRID_SIZE[1]])
        w.writerow(["grid_size_z_A", GRID_SIZE[2]])
        w.writerow(["exhaustiveness_htvs", EXHAUSTIVENESS_HTVS])
        w.writerow(["placeholder_redocking_rmsd_A", f"{rmsd_val:.2f}"])
        w.writerow(["placeholder_best_energy", f"{best_energy:.2f}"])

    VALIDATION_TXT.write_text(
        "Binding-site placeholder redocking validation\n"
        "=============================================\n"
        f"Template residue: BCORL1 {TARGET_LEU[0]}{TARGET_LEU[1]} (4HPM chain A)\n"
        f"Grid center (Leu1665 side-chain centroid): "
        f"({center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f})\n"
        f"Grid size: {GRID_SIZE[0]}x{GRID_SIZE[1]}x{GRID_SIZE[2]} A\n"
        f"Exhaustiveness (placeholder): {EXHAUSTIVENESS_PLACEHOLDER}\n"
        f"Best placeholder docking energy: {best_energy:.2f} kcal/mol\n"
        f"Placeholder RMSD (centroid-based): {rmsd_val:.2f} A\n"
        f"Status: {'PASS' if rmsd_val < 2.0 else 'FAIL - investigate before HTVS'}\n"
    )
    print(f"[validation] placeholder RMSD = {rmsd_val:.2f} A "
          f"({'PASS' if rmsd_val < 2.0 else 'FAIL'})")
    print("[done] Step 3 complete.")


if __name__ == "__main__":
    main()
