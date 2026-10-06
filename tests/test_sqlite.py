"""Tests for the SQLite storage backend (jnuherb.db).

Verifies that build_db produces a database whose tables are identical
to the CSV sources, that _load_table prefers SQLite and falls back to
CSV, and that the read-only connection helper works.
"""

from pathlib import Path
import sqlite3

import pandas as pd
import pytest

from jnuherb import _data
from jnuherb._data import (
    _TABLE_FILES,
    _load_table,
    _reset_cache,
    db_path,
    sqlite_connection,
)
from jnuherb.build_db import build_db


@pytest.fixture()
def built_db(tmp_path: Path, monkeypatch) -> Path:
    """Build jnuherb.db from the real CSVs into a temp dir and use it."""
    src = _data._data_dir()
    for fname in _TABLE_FILES.values():
        (tmp_path / fname).write_bytes((src / fname).read_bytes())
    out = build_db(tmp_path)
    monkeypatch.setenv("JNUHERB_DATA_DIR", str(tmp_path))
    _reset_cache()
    yield out
    _reset_cache()


def test_build_db_creates_file(built_db: Path):
    assert built_db.exists() and built_db.stat().st_size > 0


def test_sqlite_tables_match_csv(built_db, tmp_path):
    """Every table loaded via SQLite must equal the CSV-loaded frame."""
    src = _data._data_dir()
    for name, fname in _TABLE_FILES.items():
        df_sql = _load_table(name)
        df_csv = pd.read_csv(str(src / fname), dtype=str, keep_default_na=False)
        assert df_sql.shape == df_csv.shape, name
        assert list(df_sql.columns) == list(df_csv.columns), name
        pd.testing.assert_frame_equal(
            df_sql.reset_index(drop=True),
            df_csv.reset_index(drop=True),
            check_dtype=False,
        )


def test_meta_row_counts_valid(built_db):
    """_validate_db must pass for a freshly built database."""
    assert _data._validate_db() is True


def test_csv_fallback_when_no_db(monkeypatch):
    """Without a db, _load_table must still serve tables from CSV."""
    monkeypatch.delenv("JNUHERB_DATA_DIR", raising=False)
    monkeypatch.setattr(_data, "db_path", lambda: None)
    _reset_cache()
    df = _load_table("jnuherb_formula")
    assert len(df) > 0 and "HVPID" in df.columns
    _reset_cache()


def test_readonly_connection(built_db):
    conn = sqlite_connection()
    try:
        n = conn.execute("SELECT COUNT(*) FROM jnuherb_formula").fetchone()[0]
        assert n == len(_load_table("jnuherb_formula"))
        with pytest.raises(sqlite3.OperationalError):
            conn.execute("DROP TABLE jnuherb_formula")
    finally:
        conn.close()
