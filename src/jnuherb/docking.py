"""Molecular docking pipeline for the JNUHERB web app.

Single-ligand / single-receptor AutoDock Vina docking, built for the
website's "Molecular Docking" panel. Two entry points share one engine:

* **Automatic mode** — the receptor is a target protein selected from a
  network-pharmacology analysis result; its structure is downloaded from
  AlphaFold DB (Ensembl ID -> UniProt accession via the packaged mapping
  table ``ensembl_uniprot.csv``), and the ligand is the SMILES of a
  database chemical.
* **Upload mode** — the user provides a ligand file (SDF/MOL) or SMILES
  string and a receptor PDB file.

The docking box is defined either by user-supplied binding-site residue
numbers (centroid of their heavy atoms) or by an explicit center/size.

Heavy third-party dependencies (rdkit, meeko, vina) are imported lazily so
that the rest of the package — and the test suite — works without them.

Pipeline:
    ligand  -> RDKit 3D embed -> meeko PDBQT
    receptor -> mk_prepare_receptor (meeko CLI) -> PDBQT
    box     -> residue centroid -> Vina maps
    Vina    -> best pose -> PDBQT -> PDB (HETATM records)
    complex -> receptor PDB + ligand pose -> downloadable PDB text
"""

from __future__ import annotations

import csv
import json
import re
import subprocess
import sys
import tempfile
import urllib.request
from importlib import resources
from pathlib import Path

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction"
ALPHAFOLD_PDB = "https://alphafold.ebi.ac.uk/files/{file}"

DEFAULT_BOX_SIZE = (20.0, 20.0, 20.0)
DEFAULT_EXHAUSTIVENESS = 8
DEFAULT_NUM_POSES = 1

PACKAGE = "jnuherb"


class DockingError(Exception):
    """Raised for any user-facing docking failure (bad input, network, engine)."""


# ---------------------------------------------------------------------------
# Dependency helpers
# ---------------------------------------------------------------------------

def _require(module: str, hint: str | None = None):
    """Import *module* or raise a friendly DockingError."""
    try:
        return __import__(module, fromlist=["*"])
    except ImportError as e:
        msg = hint or f"pip install {module}"
        raise DockingError(
            f"Docking dependency '{module}' is not installed ({e}). "
            f"Install it with: {msg}"
        ) from e


def _require_vina():
    try:
        from vina import Vina  # type: ignore
        return Vina
    except ImportError as e:
        raise DockingError(
            "The AutoDock Vina Python package is not installed. "
            "It requires a C++ toolchain to build; see DEPLOY_GUIDE.md."
        ) from e


# ---------------------------------------------------------------------------
# Ensembl -> UniProt mapping (packaged CSV)
# ---------------------------------------------------------------------------

def load_uniprot_map() -> dict[str, str]:
    """Return {Ensembl_ID: UniProt_ID} from the packaged mapping table."""
    try:
        ref = resources.files(PACKAGE) / "data" / "ensembl_uniprot.csv"
        with ref.open("r", newline="", encoding="utf-8") as f:
            return {r["Ensembl_ID"]: r["UniProt_ID"] for r in csv.DictReader(f)}
    except FileNotFoundError as e:
        raise DockingError(
            "ensembl_uniprot.csv is missing. Run "
            "scripts/build_uniprot_mapping.py to generate it."
        ) from e


def ensembl_to_uniprot(ensembl_id: str) -> str:
    """Map one Ensembl protein ID to a UniProt accession."""
    up = load_uniprot_map().get(ensembl_id)
    if not up:
        raise DockingError(
            f"No UniProt mapping found for {ensembl_id}. "
            "Use upload mode and provide the protein structure (PDB) directly."
        )
    return up


# ---------------------------------------------------------------------------
# AlphaFold structure download (disk-cached; callers may wrap with st.cache_data)
# ---------------------------------------------------------------------------

def fetch_alphafold_pdb(uniprot_id: str, cache_dir: str | Path) -> Path:
    """Download the AlphaFold model for *uniprot_id* into *cache_dir*.

    Uses the AlphaFold DB prediction API to resolve the current model file
    URL, then downloads the PDB. Already-downloaded files are reused.
    """
    cache = Path(cache_dir)
    cache.mkdir(parents=True, exist_ok=True)
    dest = cache / f"AF-{uniprot_id}.pdb"
    if dest.exists() and dest.stat().st_size > 0:
        return dest

    url = f"{ALPHAFOLD_API}/{uniprot_id}"
    try:
        with urllib.request.urlopen(url, timeout=60) as resp:
            entries = json.loads(resp.read().decode())
    except Exception as e:  # noqa: BLE001 - network errors are all user-facing
        raise DockingError(
            f"Failed to query AlphaFold DB for {uniprot_id}: {e}"
        ) from e
    if not entries:
        raise DockingError(
            f"No AlphaFold model available for UniProt {uniprot_id}."
        )

    pdb_url = entries[0].get("pdbUrl")
    if not pdb_url:
        raise DockingError(f"AlphaFold DB returned no PDB URL for {uniprot_id}.")
    try:
        with urllib.request.urlopen(pdb_url, timeout=300) as resp:
            dest.write_bytes(resp.read())
    except Exception as e:  # noqa: BLE001
        raise DockingError(
            f"Failed to download AlphaFold structure for {uniprot_id}: {e}"
        ) from e
    return dest


# ---------------------------------------------------------------------------
# Binding box from residue numbers
# ---------------------------------------------------------------------------

_BACKBONE = {"N", "CA", "C", "O", "OXT"}


def parse_residue_list(text: str) -> list[int]:
    """Parse user input like '163,168,194-196' into residue numbers."""
    out: list[int] = []
    for part in re.split(r"[,，;\s]+", text.strip()):
        if not part:
            continue
        m = re.fullmatch(r"(\d+)(?:-(\d+))?", part)
        if not m:
            raise DockingError(f"Invalid residue number: {part!r}")
        a, b = int(m.group(1)), int(m.group(2) or m.group(1))
        if a > b:
            a, b = b, a
        out.extend(range(a, b + 1))
    seen, dedup = set(), []
    for r in out:
        if r not in seen:
            seen.add(r)
            dedup.append(r)
    return dedup


def _pdb_records(path: Path):
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            rec = line[:6].strip()
            if rec in ("ATOM", "HETATM") and len(line) >= 54:
                yield line


def binding_box_from_residues(
    pdb_path: str | Path,
    residues: list[int],
    box_size: tuple[float, float, float] = DEFAULT_BOX_SIZE,
) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    """Box center = centroid of side-chain heavy atoms of the given residues.

    Falls back to C-alpha atoms when a residue has no side-chain atoms,
    then to all atoms. Raises DockingError if none of the residue numbers
    exist in the structure.
    """
    wanted = set(residues)
    per_res: dict[int, tuple[list, list]] = {}
    for line in _pdb_records(Path(pdb_path)):
        resid = int(line[22:26])
        if resid not in wanted:
            continue
        side, ca = per_res.setdefault(resid, ([], []))
        atom = line[12:16].strip()
        xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        if atom == "CA":
            ca.append(xyz)
        elif atom not in _BACKBONE:
            side.append(xyz)
    if not per_res:
        raise DockingError(
            f"None of the residue numbers {residues} were found in the "
            "receptor structure. Check the numbering matches the PDB file."
        )
    # Per residue: side-chain heavy atoms if present, else C-alpha.
    coords: list[tuple[float, float, float]] = []
    for side, ca in per_res.values():
        coords.extend(side if side else ca)
    if not coords:
        raise DockingError(
            f"None of the residue numbers {residues} were found in the "
            "receptor structure. Check the numbering matches the PDB file."
        )
    n = len(coords)
    center = tuple(round(sum(c[i] for c in coords) / n, 3) for i in range(3))
    return center, box_size


def default_box_from_receptor(
    pdb_path: str | Path,
    box_size: tuple[float, float, float] = DEFAULT_BOX_SIZE,
) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    """Fallback box centered on the receptor's C-alpha centroid.

    Used when the user gives neither binding-site residues nor an explicit
    center (blind docking). Falls back to all heavy atoms for structures
    without C-alpha records.
    """
    cas, heavy = [], []
    for line in _pdb_records(Path(pdb_path)):
        atom = line[12:16].strip()
        if atom == "CA":
            cas.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
        elif atom != "H":
            heavy.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    coords = cas or heavy
    if not coords:
        raise DockingError(f"No ATOM records found in receptor: {pdb_path}")
    n = len(coords)
    center = tuple(round(sum(c[i] for c in coords) / n, 3) for i in range(3))
    return center, box_size


# ---------------------------------------------------------------------------
# Ligand preparation
# ---------------------------------------------------------------------------

def _mol_to_pdbqt(mol, ligand_name: str) -> str:
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    preparer = MoleculePreparation()
    setups = preparer.prepare(mol)
    if not setups:
        raise DockingError(
            f"Meeko could not prepare ligand {ligand_name!r} "
            "(unusual chemistry or missing 3D coordinates?)."
        )
    result = PDBQTWriterLegacy.write_string(setups[0])
    # meeko >= 0.6 returns (pdbqt, success, error); older versions return str
    if isinstance(result, tuple):
        pdbqt, success, error_msg = result
        if not success:
            raise DockingError(
                f"Meeko failed to write ligand {ligand_name!r}: {error_msg}"
            )
    else:
        pdbqt = result
    return f"REMARK LIGAND_NAME: {ligand_name}\n" + pdbqt


def ligand_pdbqt_from_smiles(smiles: str, ligand_name: str = "LIG") -> str:
    """SMILES -> 3D embedded, protonated, PDBQT string."""
    _require("rdkit", "pip install rdkit")
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit import RDLogger

    RDLogger.DisableLog("rdApp.*")
    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None:
        raise DockingError(f"Invalid SMILES string: {smiles!r}")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    if AllChem.EmbedMolecule(mol, params) == -1:
        raise DockingError(
            f"Failed to generate a 3D conformation for {ligand_name!r}."
        )
    AllChem.UFFOptimizeMolecule(mol, maxIters=500)
    return _mol_to_pdbqt(mol, ligand_name)


def ligand_pdbqt_from_bytes(data: bytes, filename: str,
                            ligand_name: str = "LIG") -> str:
    """SDF/MOL file content (bytes) -> PDBQT string."""
    _require("rdkit", "pip install rdkit")
    from rdkit import Chem
    from rdkit import RDLogger

    RDLogger.DisableLog("rdApp.*")
    suffix = Path(filename).suffix.lower()
    if suffix in (".mol", ".sdf"):
        mol = Chem.MolFromMolBlock(data.decode("utf-8", errors="replace"),
                                   removeHs=False, sanitize=False)
    elif suffix in (".pdb",):
        mol = Chem.MolFromPDBBlock(data.decode("utf-8", errors="replace"),
                                   removeHs=False, sanitize=False)
    else:
        # Last resort: try SMILES (one per line, first molecule wins)
        text = data.decode("utf-8", errors="replace").strip().splitlines()[0]
        return ligand_pdbqt_from_smiles(text, ligand_name)

    if mol is None:
        raise DockingError(f"Could not parse ligand file: {filename}")
    mol.UpdatePropertyCache(strict=False)
    mol = Chem.AddHs(mol, addCoords=True)
    return _mol_to_pdbqt(mol, ligand_name)


# ---------------------------------------------------------------------------
# Receptor preparation (meeko mk_prepare_receptor CLI)
# ---------------------------------------------------------------------------

def prepare_receptor(pdb_path: str | Path, workdir: str | Path) -> Path:
    """Receptor PDB -> PDBQT using meeko's mk_prepare_receptor."""
    pdb_path = Path(pdb_path)
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    out_base = workdir / "receptor"

    import shutil
    script = shutil.which("mk_prepare_receptor.py")
    if script is None:
        script = Path(sys.executable).parent / "mk_prepare_receptor.py"
    if not Path(script).exists():
        raise DockingError(
            "mk_prepare_receptor.py not found — install meeko with "
            "`pip install meeko gemmi scipy`."
        )
    # -p/--write_pdbqt is explicit in meeko >= 0.8 (-o alone writes nothing)
    proc = subprocess.run(
        [sys.executable, str(script), "-i", str(pdb_path), "-o", str(out_base),
         "-p"],
        capture_output=True, text=True, timeout=600,
    )
    pdbqt = out_base.with_suffix(".pdbqt")
    if proc.returncode != 0 or not pdbqt.exists():
        tail = (proc.stderr or proc.stdout or "").strip()[-500:]
        raise DockingError(
            f"Receptor preparation failed for {pdb_path.name}: {tail}"
        )
    return pdbqt


# ---------------------------------------------------------------------------
# Docking
# ---------------------------------------------------------------------------

def run_docking(
    receptor_pdbqt: str | Path,
    ligand_pdbqt: str,
    center: tuple[float, float, float],
    box_size: tuple[float, float, float],
    exhaustiveness: int = DEFAULT_EXHAUSTIVENESS,
    num_poses: int = DEFAULT_NUM_POSES,
    workdir: str | Path | None = None,
) -> tuple[float, Path]:
    """Run Vina. Returns (best affinity kcal/mol, pose PDBQT path)."""
    Vina = _require_vina()
    receptor_pdbqt = Path(receptor_pdbqt)
    if not receptor_pdbqt.exists():
        raise DockingError(f"Receptor PDBQT not found: {receptor_pdbqt}")

    tmp = None
    if workdir is None:
        tmp = tempfile.TemporaryDirectory(prefix="jnuherb_dock_")
        workdir = Path(tmp.name)
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    pose_path = workdir / "pose.pdbqt"

    v = Vina(sf_name="vina", verbosity=0, seed=42)
    v.set_receptor(str(receptor_pdbqt))
    v.set_ligand_from_string(ligand_pdbqt)
    v.compute_vina_maps(center=list(center), box_size=list(box_size))
    v.dock(exhaustiveness=exhaustiveness, n_poses=num_poses)
    energies = v.energies()
    if len(energies) == 0:
        raise DockingError("Vina produced no poses.")
    affinity = float(energies[0][0])
    v.write_poses(str(pose_path), n_poses=num_poses, overwrite=True)
    # Keep the temp dir alive until the caller has read the file.
    if tmp is not None:
        pose_path._tmpdir = tmp  # type: ignore[attr-defined]
    return affinity, pose_path


# ---------------------------------------------------------------------------
# PDBQT pose -> PDB, and complex assembly
# ---------------------------------------------------------------------------

# AutoDock atom types that differ from the element symbol in the PDBQT
# element column (vina writes types like OA/NA/HD, not true elements).
_AD_TYPE_TO_ELEMENT = {
    "OA": "O", "NA": "N", "SA": "S", "HD": "H",
    "CL": "Cl", "BR": "Br", "FE": "Fe", "ZN": "Zn",
    "CA": "Ca", "MG": "Mg", "MN": "Mn", "CU": "Cu",
}


def _pdbqt_element(line: str) -> str:
    raw = line[76:78].strip() if len(line) >= 78 else ""
    if raw in _AD_TYPE_TO_ELEMENT:
        return _AD_TYPE_TO_ELEMENT[raw]
    if raw:
        # True element symbols are 1-2 chars, e.g. "C", "N", "Cl", "Na".
        m = re.fullmatch(r"[A-Z][a-z]?", raw)
        if m:
            return raw
    # Infer from the atom name, PDB-style (leading letter, else first letter).
    atom = line[12:16]
    m = re.search(r"[A-Za-z]", atom)
    return m.group(0).upper() if m else "C"


def pose_pdbqt_to_pdb(pose_path: str | Path, ligand_name: str = "LIG") -> str:
    """Convert a (single-model) Vina pose PDBQT into PDB HETATM records.

    The ligand becomes residue LIG, chain Z, residue number 1. Coordinates
    come straight from the PDBQT ATOM/HETATM lines.
    """
    lines_out = [f"REMARK LIGAND {ligand_name} (Vina docked pose)"]
    serial = 1
    for line in _pdb_records(Path(pose_path)):
        atom = line[12:16]
        elem = _pdbqt_element(line)
        x, y, z = line[30:54].split()
        lines_out.append(
            f"HETATM{serial:5d} {atom}{'LIG':>4s} Z{1:4d}    "
            f"{float(x):8.3f}{float(y):8.3f}{float(z):8.3f}"
            f"  1.00  0.00          {elem:>2s}"
        )
        serial += 1
    lines_out.append("TER")
    return "\n".join(lines_out) + "\n"


def make_complex_pdb(
    receptor_pdb: str | Path,
    pose_pdbqt: str | Path,
    ligand_name: str = "LIG",
    affinity: float | None = None,
) -> str:
    """Merge receptor PDB and docked ligand pose into one complex PDB text."""
    receptor_text = Path(receptor_pdb).read_text(encoding="utf-8",
                                                 errors="replace")
    header = [f"REMARK JNUHERB molecular docking complex"]
    if affinity is not None:
        header.append(f"REMARK VINA_AFFINITY_KCAL_MOL: {affinity:.2f}")
    header.append(f"REMARK LIGAND: {ligand_name}")
    header.append("REMARK RECEPTOR: " + Path(receptor_pdb).name)

    ligand_pdb = pose_pdbqt_to_pdb(pose_pdbqt, ligand_name)
    body = receptor_text.rstrip("\n").split("\n")
    # Drop trailing END/TER from the receptor; ligand block ends with TER+END.
    while body and body[-1].startswith(("END", "TER")):
        body.pop()
    return "\n".join(header + body) + "\n" + ligand_pdb + "END\n"


def dock_smiles_to_protein(
    smiles: str,
    ligand_name: str,
    receptor_pdb: str | Path,
    residues: list[int] | None = None,
    center: tuple[float, float, float] | None = None,
    box_size: tuple[float, float, float] = DEFAULT_BOX_SIZE,
    exhaustiveness: int = DEFAULT_EXHAUSTIVENESS,
    workdir: str | Path | None = None,
) -> dict:
    """End-to-end convenience wrapper used by the web app.

    Returns dict with affinity, complex PDB text, pose path, and box used.
    """
    receptor_pdb = Path(receptor_pdb)
    if residues is not None:
        center, box_size = binding_box_from_residues(receptor_pdb, residues,
                                                     box_size)
    if center is None:
        raise DockingError("Either binding-site residues or an explicit "
                           "box center must be provided.")

    tmp = None
    if workdir is None:
        tmp = tempfile.TemporaryDirectory(prefix="jnuherb_dock_")
        workdir = Path(tmp.name)
    workdir = Path(workdir)

    try:
        ligand_pdbqt = ligand_pdbqt_from_smiles(smiles, ligand_name)
        receptor_pdbqt = prepare_receptor(receptor_pdb, workdir)
        affinity, pose_path = run_docking(
            receptor_pdbqt, ligand_pdbqt, center, box_size,
            exhaustiveness=exhaustiveness, num_poses=DEFAULT_NUM_POSES,
            workdir=workdir,
        )
        complex_pdb = make_complex_pdb(receptor_pdb, pose_path,
                                       ligand_name, affinity)
    finally:
        if tmp is not None:
            tmp.cleanup()

    return {
        "affinity": affinity,
        "center": center,
        "box_size": box_size,
        "complex_pdb": complex_pdb,
    }
