"""Internal data loading with lazy caching.

Tables are read from the bundled SQLite database (``jnuherb.db``) when it
is present — built from the CSV sources by ``python -m jnuherb.build_db`` —
and fall back to parsing the CSV files directly. All columns stay TEXT to
preserve the historical ``dtype=str, keep_default_na=False`` semantics, so
callers behave identically regardless of the backend.

Set ``JNUHERB_DATA_DIR`` to point at a different data directory (used by
tests and custom deployments).
"""

from __future__ import annotations

import os
import sqlite3
from importlib.resources import files
from pathlib import Path

import pandas as pd

_cache: dict[str, pd.DataFrame] = {}
_db_valid: bool | None = None

_TABLE_FILES: dict[str, str] = {
    "jnuherb_formula": "jnuherb_formula.csv",
    "jnuherb_tcm": "jnuherb_tcm.csv",
    "jnuherb_chemicals": "jnuherb_chemicals.csv",
    "jnuherb_proteins": "jnuherb_proteins.csv",
    "jnuherb_formula_tcm_links": "jnuherb_formula_tcm_links.csv",
    "jnuherb_tcm_chemical_links": "jnuherb_tcm_chemical_links.csv",
    "jnuherb_chemical_protein_links": "jnuherb_chemical_protein_links.csv",
    "jnuherb_protein_disease_links": "jnuherb_protein_disease_links.csv.gz",
}

# Columns indexed in jnuherb.db: primary IDs, join keys and the name
# columns used by analysis queries and AI retrieval.
_ID_INDEXES: dict[str, list[str]] = {
    "jnuherb_formula": ["HVPID", "name"],
    "jnuherb_tcm": ["HVMID", "cn_name", "pinyin_name", "en_name"],
    "jnuherb_chemicals": ["HVCID", "Name", "CAS_id", "PubChem_id"],
    "jnuherb_proteins": ["Ensembl_ID", "gene_name"],
    "jnuherb_formula_tcm_links": ["HVPID", "HVMID"],
    "jnuherb_tcm_chemical_links": ["HVMID", "HVCID"],
    "jnuherb_chemical_protein_links": ["HVCID", "Ensembl_ID"],
    "jnuherb_protein_disease_links": ["Ensembl_ID", "disease_id", "disease_name"],
}


def _data_dir() -> Path:
    """Filesystem directory holding the bundled data files."""
    override = os.environ.get("JNUHERB_DATA_DIR")
    if override:
        return Path(override)
    traversable = files("jnuherb.data")
    if isinstance(traversable, Path):
        return traversable
    # Namespace package (src tree on sys.path): resolve via the parent package.
    import jnuherb

    return Path(jnuherb.__file__).parent / "data"


def db_path() -> Path | None:
    """Path to jnuherb.db if it exists in the data directory, else None."""
    p = _data_dir() / "jnuherb.db"
    return p if p.exists() else None


def sqlite_connection() -> sqlite3.Connection:
    """Read-only connection to jnuherb.db for indexed SQL lookups.

    Raises:
        FileNotFoundError: if the database has not been built.
    """
    p = db_path()
    if p is None:
        raise FileNotFoundError(
            "jnuherb.db not found; run `python -m jnuherb.build_db`"
        )
    return sqlite3.connect(f"file:{p}?mode=ro", uri=True)


def _validate_db() -> bool:
    """Check _meta row counts against actual counts (once per process)."""
    global _db_valid
    if _db_valid is None:
        _db_valid = False
        try:
            conn = sqlite_connection()
            try:
                meta = pd.read_sql_query("SELECT table_name, nrows FROM _meta", conn)
                for row in meta.itertuples(index=False):
                    actual = conn.execute(
                        f'SELECT COUNT(*) FROM "{row.table_name}"'
                    ).fetchone()[0]
                    if actual != row.nrows:
                        return _db_valid  # stays False
                _db_valid = True
            finally:
                conn.close()
        except sqlite3.Error:
            _db_valid = False
    return _db_valid


def _load_csv(name: str) -> pd.DataFrame:
    path = _data_dir() / _TABLE_FILES[name]
    return pd.read_csv(str(path), dtype=str, keep_default_na=False)


def _load_sqlite(name: str) -> pd.DataFrame:
    conn = sqlite_connection()
    try:
        return pd.read_sql_query(f'SELECT * FROM "{name}"', conn, dtype=str)
    finally:
        conn.close()


def _load_table(name: str) -> pd.DataFrame:
    """Lazy-load a bundled table with in-memory caching.

    Prefers the SQLite database; falls back to the CSV source when the
    database is missing or fails validation.

    Args:
        name: Logical table name (key of _TABLE_FILES).

    Returns:
        DataFrame with all columns as strings.
    """
    if name not in _cache:
        if db_path() is not None and _validate_db():
            _cache[name] = _load_sqlite(name)
        else:
            _cache[name] = _load_csv(name)
    return _cache[name]


def _get_data_path(filename: str) -> str:
    """Get the absolute filesystem path to a bundled data file.

    Args:
        filename: CSV file name, e.g. ``"jnuherb_tcm.csv"``.

    Returns:
        Absolute path as a string.
    """
    return str(_data_dir() / filename)


def _reset_cache() -> None:
    """Clear the in-memory table cache and db validation flag (tests)."""
    global _db_valid
    _cache.clear()
    _db_valid = None
