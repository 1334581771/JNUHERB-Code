"""Step 2: Prepare ligand library (2D SDF -> 3D protonated SDF -> per-molecule PDBQT).

Two-phase design to survive WSL process lifecycle:
  Phase A: iterate SDF, embed 3D, write to intermediate SDF (checkpointed every
           CHECKPOINT_EVERY molecules so progress survives interruption).
  Phase B: convert all 3D molecules to per-molecule PDBQT via meeko.

Outputs:
  data/processed/ligands_3d.sdf
  data/processed/ligands_pdbqt/*.pdbqt
  data/processed/ligand_prep_report.csv
"""
from pathlib import Path
import csv
import sys

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[2]
SDF_IN = ROOT / "data" / "L6020.sdf"
OUT_DIR = ROOT / "data" / "processed"
PDBQT_DIR = OUT_DIR / "ligands_pdbqt"
PDBQT_DIR.mkdir(parents=True, exist_ok=True)

SDF_3D = OUT_DIR / "ligands_3d.sdf"
SDF_3D_TMP = OUT_DIR / "ligands_3d.tmp.sdf"
REPORT = OUT_DIR / "ligand_prep_report.csv"
REPORT_TMP = OUT_DIR / "ligand_prep_report.tmp.csv"

CHECKPOINT_EVERY = 200

# SDF property fields to use as compound ID, in priority order, when the
# SDF title line is blank (RDKit still sets _Name="" in that case, so a
# bare HasProp("_Name") check is not enough).
ID_PROPS = ("CAS", "Name", "ID", "IDNUMBER", "hitid", "catalog_no", "inchlkey")


def get_cid(mol, i):
    """Robust compound ID: title first, then property fields, then fallback."""
    if mol is None:
        return f"mol_{i:05d}"
    cid = mol.GetProp("_Name").strip() if mol.HasProp("_Name") else ""
    if not cid:
        for k in ID_PROPS:
            if mol.HasProp(k) and mol.GetProp(k).strip():
                cid = mol.GetProp(k).strip()
                break
    return cid or f"mol_{i:05d}"


def embed_3d(mol, max_tries=20):
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.maxIterations = 2000
    cid = AllChem.EmbedMolecule(mol, params)
    if cid == -1:
        params.useRandomCoords = True
        cid = AllChem.EmbedMolecule(mol, params)
    if cid == -1:
        return None
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    except Exception:
        pass
    return mol


def phase_a():
    """2D -> 3D embedding with checkpoint."""
    suppl = Chem.SDMolSupplier(str(SDF_IN), removeHs=False)
    writer = Chem.SDWriter(str(SDF_3D_TMP))
    report_f = open(REPORT_TMP, "w", newline="")
    w = csv.writer(report_f)
    w.writerow(["compound_id", "status", "note"])

    n_ok = n_fail = 0
    for i, mol in enumerate(suppl):
        cid = get_cid(mol, i)
        if mol is None:
            w.writerow([cid, "FAIL_parse", ""])
            n_fail += 1
        else:
            try:
                mol3d = embed_3d(mol)
                if mol3d is None:
                    w.writerow([cid, "FAIL_embed", ""])
                    n_fail += 1
                else:
                    mol3d.SetProp("_Name", cid)
                    writer.write(mol3d)
                    w.writerow([cid, "OK", ""])
                    n_ok += 1
            except Exception as e:
                w.writerow([cid, "FAIL_error", str(e)[:80]])
                n_fail += 1

        if (i + 1) % CHECKPOINT_EVERY == 0:
            writer.flush()
            report_f.flush()
            print(f"[phase A] {i+1} processed, {n_ok} OK, {n_fail} FAIL",
                  flush=True)

    writer.close()
    report_f.close()
    SDF_3D_TMP.rename(SDF_3D)
    REPORT_TMP.rename(REPORT)
    print(f"[phase A done] total={n_ok+n_fail}, OK={n_ok}, FAIL={n_fail}",
          flush=True)
    return n_ok


def phase_b():
    """3D SDF -> per-molecule PDBQT."""
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    ok = fail = 0
    suppl = Chem.SDMolSupplier(str(SDF_3D), removeHs=False)
    rows = list(csv.reader(open(REPORT)))
    status_map = {r[0]: r[1] for r in rows[1:]}

    for i, mol in enumerate(suppl):
        if mol is None:
            continue
        name = get_cid(mol, i)
        if status_map.get(name) != "OK":
            continue
        out_path = PDBQT_DIR / f"{name}.pdbqt"
        if out_path.exists():
            ok += 1
            continue
        try:
            prep = MoleculePreparation()
            mol_setup = prep.prepare(mol)[0]
            pdbqt_string, _, _ = PDBQTWriterLegacy.write_string(mol_setup)
            out_path.write_text(pdbqt_string)
            ok += 1
        except Exception:
            status_map[name] = "FAIL_pdbqt"
            fail += 1

    # rewrite report with updated statuses
    with open(REPORT, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["compound_id", "status", "note"])
        for r in rows[1:]:
            if r[0] in status_map:
                w.writerow([r[0], status_map[r[0]], r[2] if len(r) > 2 else ""])
            else:
                w.writerow(r)

    print(f"[phase B done] PDBQT OK={ok}, FAIL={fail}", flush=True)


if __name__ == "__main__":
    print("=== Phase A: 2D -> 3D ===", flush=True)
    phase_a()
    print("=== Phase B: 3D -> PDBQT ===", flush=True)
    phase_b()
    print("[done] Step 2 complete.", flush=True)
