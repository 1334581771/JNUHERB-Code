"""Step 4: Run high-throughput virtual screening with Vina.

Reads vina_config.txt (center, size, exhaustiveness, num_modes) and all
ligand PDBQT files, docks each into the binding pocket using multiprocessing.

Outputs:
  results/models/<compound_id>.pdbqt      docked pose file, ONE per ligand,
                                          containing only the top KEEP_TOP_POSES
  results/tables/docking_results_raw.csv  scores for ALL candidate poses
                                          (up to num_modes rows per ligand),
                                          checkpointed, resumable

compound_id is the unique ligand identifier (PDBQT filename stem, truncated
to 80 chars by step 2). ligand_name is the FULL untruncated name recovered
from the source SDF's property fields (Name > ID > hitid > CAS > title), and
is also embedded in each pose file as a "REMARK LIGAND_NAME:" line, so both
the CSV and the docked structures carry the correct ligand name.

Performance notes:
  - the vina package is CPU-only; this script never touches the GPU
  - parallel across ligands via multiprocessing (near-linear scaling);
    default workers = min(cpu_count, 28) -- RAM usage is ~tens of MB per
    worker, so memory is not the bottleneck on a 96G machine

On resume, existing pose files are skipped and their full multi-pose scores
are carried over from the previous CSV (the pose file itself only stores the
top pose, so the CSV is the authoritative score record).

Features:
  - pose files double as checkpoint: existing non-empty outputs are skipped
    on rerun, and their scores are re-parsed from the REMARK lines so the
    CSV stays complete
  - CSV is flushed every 200 ligands, so a crash loses at most the last batch
  - visual progress bar (no tqdm dependency)
"""
from pathlib import Path
import csv
import os
import re
import time
import multiprocessing as mp
from functools import partial

from vina import Vina

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "processed"
TABLES = ROOT / "results" / "tables"
MODELS = ROOT / "results" / "models"
TABLES.mkdir(parents=True, exist_ok=True)
MODELS.mkdir(parents=True, exist_ok=True)
RESULTS_CSV = TABLES / "docking_results_raw.csv"

LIGAND_DIR = DATA / "ligands_pdbqt"
CONFIG_FILE = DATA / "vina_config.txt"
# source of the full ligand names (must be the SDF step 2 was run on)
LIGAND_SDF = ROOT / "data" / "L6020.sdf"

# filename stem -> full ligand name, built in main() before the Pool forks
NAME_MAP = {}

# how many top poses to keep in each pose FILE (CSV always keeps all poses)
KEEP_TOP_POSES = 1

CSV_HEADER = ["compound_id", "ligand_name", "pose_rank", "vina_affinity_kcal_mol",
              "rmsd_lb", "rmsd_ub", "error"]


def clean_name(compound_id: str) -> str:
    """Strip the step-2 dedup suffix (_NNNNN) for a human-readable name."""
    return re.sub(r"_\d{5}$", "", compound_id)


def build_name_map(sdf_path: Path):
    """Map pdbqt filename stem -> FULL ligand name from the source SDF.

    Replays step 2's exact naming (sanitize + truncation + dedup) so every
    existing pdbqt file can be matched back to its SDF record; the full name
    comes from property fields (Name > ID > hitid > CAS > title), NOT the
    80-char-truncated filename.

    Uses RDKit (not pybel) to read the SDF: pybel's SWIG data iterator can
    abort the whole process with 'swig::stop_iteration', especially under
    the GC pressure of a 16k-molecule loop.
    """
    from rdkit import Chem
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")

    # keep in sync with sdf2d_to_pdbqt.py -- CAS first, and NO "Name"
    # (the IUPAC Name field is only used for the human-readable ligand_name)
    ID_KEYS = ("CAS", "ID", "Catalog_NO", "IDNUMBER",
               "hitid", "catalog_no")
    FULL_NAME_KEYS = ("Name", "ID", "hitid", "CAS")

    def sanitize(name, idx, max_len=80):
        cleaned = re.sub(r'[\\/*?:"<>|]', "_", str(name)).strip()
        if len(cleaned) > max_len:
            cleaned = cleaned[:max_len - 6] + f"_{idx:05d}"
        return cleaned

    name_map = {}
    if not sdf_path.exists():
        print(f"[warn] {sdf_path} not found; ligand_name falls back to "
              f"truncated filename stem")
        return name_map

    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=False,
                               strictParsing=False)
    seen = set()
    for idx, mol in enumerate(suppl):
        if mol is None:
            continue
        lower_map = {p.lower(): p for p in mol.GetPropNames()
                     if p != "_Name"}

        def safe_get(prop):
            """Some vendor SDFs carry non-UTF8 bytes (e.g. GBK comments);
            RDKit decodes strictly and would raise UnicodeDecodeError."""
            try:
                return mol.GetProp(prop).strip()
            except UnicodeDecodeError:
                return ""

        # compound id exactly as sdf2d_to_pdbqt.py computed it:
        # CAS > other ID fields > SDF title > mol_XXXXX
        cid = ""
        for key in ID_KEYS:
            real = lower_map.get(key.lower())
            if real is not None and safe_get(real):
                cid = safe_get(real)
                break
        if not cid:
            cid = safe_get("_Name") if mol.HasProp("_Name") else ""
        base = cid
        cid = sanitize(base or f"mol_{idx:05d}", idx)
        while cid in seen:
            cid = sanitize(f"{base}_{idx:05d}", idx)
        seen.add(cid)
        # full human-readable name for CSV + pose REMARK
        full = ""
        for key in FULL_NAME_KEYS:
            real = lower_map.get(key.lower())
            if real is not None and safe_get(real):
                full = safe_get(real)
                break
        name_map[cid] = full or cid
    return name_map


def load_config(path: Path):
    cfg = {}
    for line in path.read_text().splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            cfg[k.strip()] = v.strip()
    return cfg


def pose_file_rows(name):
    """Rebuild the rank-1 score row from an existing pose file.

    Used on resume when a ligand's pose file exists but its CSV row was
    lost (previous run interrupted between completion and the 200-ligand
    CSV flush). The affinity is recovered from the REMARK VINA RESULT line.
    """
    full_name = NAME_MAP.get(name, clean_name(name))
    pf = MODELS / f"{name}.pdbqt"
    rows = []
    try:
        for line in pf.read_text(errors="replace").splitlines():
            if line.startswith("REMARK VINA RESULT:"):
                aff = round(float(line.split()[3]), 3)
                rows.append([name, full_name, len(rows) + 1, aff,
                             "", "", ""])
    except Exception:
        pass
    return rows


def load_carried_rows(path: Path):
    """Read a previous run's CSV so resumed ligands keep their full score
    history (the pose file alone only stores the top pose). Older 6-column
    CSVs (without ligand_name) are normalized on the fly."""
    carried = {}
    if path.exists():
        with open(path, newline="") as f:
            for row in csv.reader(f):
                if not row or row[0] == CSV_HEADER[0]:
                    continue
                if len(row) == len(CSV_HEADER) - 1:  # legacy format
                    row = row[:1] + [clean_name(row[0])] + row[1:]
                carried.setdefault(row[0], []).append(row)
    return carried


def dock_one(ligand_path: Path, cfg: dict):
    """Dock a single ligand.

    Writes only the top KEEP_TOP_POSES pose(s) to the pose file, but returns
    score rows for ALL candidate poses (up to n_poses).

    Returns (name, rows, status) where status is "OK" | "SKIP" | "FAIL".
    """
    name = ligand_path.stem
    full_name = NAME_MAP.get(name, clean_name(name))
    out_path = MODELS / f"{name}.pdbqt"
    n_poses = int(cfg.get("num_modes", 9))

    # checkpoint/resume: keep existing poses; scores are recovered by main()
    # from the previous CSV (pose file only holds the top pose)
    if out_path.exists() and out_path.stat().st_size > 0:
        return name, [], "SKIP"

    try:
        v = Vina(sf_name="vina", verbosity=0, seed=42)
        v.set_receptor(cfg["receptor"])
        v.set_ligand_from_file(str(ligand_path))
        v.compute_vina_maps(
            center=[float(cfg["center_x"]), float(cfg["center_y"]),
                    float(cfg["center_z"])],
            box_size=[float(cfg["size_x"]), float(cfg["size_y"]),
                      float(cfg["size_z"])])
        v.dock(exhaustiveness=int(cfg["exhaustiveness"]), n_poses=n_poses)
        energies = v.energies()          # must read before poses are written
        v.write_poses(str(out_path), n_poses=KEEP_TOP_POSES, overwrite=True)
        # embed the full ligand name in the pose file
        out_path.write_text(
            f"REMARK LIGAND_NAME: {full_name}\n" + out_path.read_text())
        rows = [
            [name, full_name, pose_rank, round(float(e[0]), 3),
             round(float(e[1]), 3), round(float(e[2]), 3), ""]
            for pose_rank, e in enumerate(energies, start=1)
        ]
        return name, rows, "OK"
    except Exception as e:
        return name, [[name, full_name, -1, None, None, None,
                       str(e)[:100]]], "FAIL"


# --------------------------------------------------------------------------- #
# progress bar (zero-dependency)
# --------------------------------------------------------------------------- #
def draw_progress(done: int, total: int, t0: float, stats: dict, width: int = 40):
    el = time.time() - t0
    frac = done / total if total else 1.0
    filled = int(width * frac)
    bar = "=" * filled + ">" + "-" * (width - filled - 1) if filled < width else "=" * width
    rate = done / el if el > 0 else 0
    eta = (total - done) / rate if rate > 0 else 0
    print(f"\r[{bar}] {frac*100:5.1f}% {done}/{total} | "
          f"OK={stats['OK']} SKIP={stats['SKIP']} FAIL={stats['FAIL']} | "
          f"{rate:.1f}/s ETA {eta/60:.1f}min   ", end="", flush=True)


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Step 4: HTVS with Vina")
    parser.add_argument("-w", "--workers", type=int,
                        default=min(mp.cpu_count(), 28),
                        help="parallel docking workers "
                             "(default: min(cpu_count, 28))")
    parser.add_argument("-l", "--ligand-dir", default=str(LIGAND_DIR),
                        help="directory of prepared ligand PDBQT files")
    args = parser.parse_args()

    cfg = load_config(CONFIG_FILE)
    ligand_dir = Path(args.ligand_dir)
    ligands = sorted(ligand_dir.glob("*.pdbqt"))
    print(f"[htvs] {len(ligands)} ligands from {ligand_dir}; "
          f"exhaustiveness={cfg['exhaustiveness']}; "
          f"num_modes={cfg.get('num_modes', 9)}; "
          f"workers={args.workers}")
    print(f"[htvs] poses -> {MODELS}")

    # each worker is already a separate process; pin OpenMP to 1 thread so
    # vina's internals don't oversubscribe the cores
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")

    # full ligand names for CSV + pose-file REMARKs (workers inherit via fork)
    global NAME_MAP
    NAME_MAP = build_name_map(LIGAND_SDF)
    print(f"[htvs] name map: {len(NAME_MAP)} ligands from {LIGAND_SDF.name}")

    worker = partial(dock_one, cfg=cfg)
    n_proc = max(1, min(args.workers, len(ligands) or 1))
    stats = {"OK": 0, "SKIP": 0, "FAIL": 0}
    carried = load_carried_rows(RESULTS_CSV)
    for cid, rows in carried.items():          # backfill full names on resume
        if cid in NAME_MAP:
            for r in rows:
                r[1] = NAME_MAP[cid]
    rows_all = []
    t0 = time.time()

    def flush_csv():
        with open(RESULTS_CSV, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(CSV_HEADER)
            w.writerows(rows_all)

    with mp.Pool(n_proc) as pool:
        for i, (name, rows, status) in enumerate(
                pool.imap_unordered(worker, ligands), 1):
            stats[status] += 1
            if status == "SKIP":
                # keep full score history; fall back to parsing the pose
                # file if the previous CSV row was never flushed
                rows = carried.get(name) or pose_file_rows(name)
            else:
                carried.pop(name, None)          # drop stale rows on re-dock
            rows_all.extend(rows)
            if i % 200 == 0 or i == len(ligands):
                flush_csv()  # checkpoint: a crash loses at most the last batch
            draw_progress(i, len(ligands), t0, stats)
        print()  # newline after the bar

    n_ok_ligands = stats["OK"] + stats["SKIP"]
    n_poses = sum(1 for _ in MODELS.glob("*.pdbqt"))
    print(f"[done] {n_ok_ligands}/{len(ligands)} ligands docked "
          f"(OK={stats['OK']}, resumed={stats['SKIP']}, FAIL={stats['FAIL']})")
    print(f"[done] {n_poses} pose files -> {MODELS}")
    print(f"[done] scores -> {RESULTS_CSV.name}")


if __name__ == "__main__":
    main()
