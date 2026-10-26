"""Build the bundled SQLite database from the CSV table sources.

Usage:
    python -m jnuherb.build_db [--data-dir DIR]

Creates ``jnuherb.db`` in the package data directory with one TEXT-typed
table per CSV (preserving the ``dtype=str`` semantics the app relies on),
row-count metadata in ``_meta`` for startup validation, and indexes on the
ID, join and name columns used by analysis queries and AI retrieval.

The build is atomic: rows are written to a temporary file that replaces
the target only on success, so a failed rebuild never clobbers a good db.
"""

from __future__ import annotations

import argparse
import os
import sqlite3
from pathlib import Path

import pandas as pd

from jnuherb._data import _ID_INDEXES, _TABLE_FILES, _data_dir


def build_db(data_dir: Path | None = None) -> Path:
    """Build jnuherb.db from the CSV tables in ``data_dir``.

    Args:
        data_dir: Directory holding the CSV files. Defaults to the
            package data directory (or ``JNUHERB_DATA_DIR`` if set).

    Returns:
        Path of the created database file.
    """
    data_dir = Path(data_dir) if data_dir else _data_dir()
    db_path = data_dir / "jnuherb.db"
    tmp_path = data_dir / "jnuherb.db.tmp"
    if tmp_path.exists():
        tmp_path.unlink()

    conn = sqlite3.connect(tmp_path)
    try:
        meta_rows: list[tuple[str, int]] = []
        for name, fname in _TABLE_FILES.items():
            csv_path = data_dir / fname
            if not csv_path.exists():
                raise FileNotFoundError(f"missing source table: {csv_path}")
            df = pd.read_csv(str(csv_path), dtype=str, keep_default_na=False)
            text_dtype = {c: "TEXT" for c in df.columns}
            df.to_sql(name, conn, index=False, if_exists="replace", dtype=text_dtype)
            for col in _ID_INDEXES.get(name, []):
                if col in df.columns:
                    conn.execute(
                        f'CREATE INDEX IF NOT EXISTS "ix_{name}_{col}" '
                        f'ON "{name}"("{col}")'
                    )
            meta_rows.append((name, len(df)))
        # Case-insensitive index on the reverse-lookup name column: query.py
        # matches disease names with COLLATE NOCASE, which a BINARY index
        # cannot serve, and the table has ~1.7M rows.
        conn.execute(
            'CREATE INDEX IF NOT EXISTS '
            '"ix_jnuherb_protein_disease_links_disease_name_ci" '
            'ON "jnuherb_protein_disease_links"("disease_name" COLLATE NOCASE)'
        )
        # Pre-aggregated disease directory powering the reverse-lookup
        # suggestion box: one row per distinct (disease_id, disease_name),
        # so the UI's substring search scans ~17k rows instead of
        # aggregating 1.7M on every keystroke.
        conn.execute(
            'CREATE TABLE IF NOT EXISTS "jnuherb_disease_lookup" AS '
            'SELECT disease_id, disease_name, '
            'COUNT(DISTINCT "Ensembl_ID") AS n_links, '
            "MAX(CAST(NULLIF(score, '') AS REAL)) AS best_score, "
            'GROUP_CONCAT(DISTINCT source) AS sources '
            'FROM "jnuherb_protein_disease_links" '
            'GROUP BY disease_id, disease_name'
        )
        conn.execute(
            'CREATE INDEX IF NOT EXISTS "ix_jnuherb_disease_lookup_name" '
            'ON "jnuherb_disease_lookup"("disease_name")'
        )
        conn.execute(
            'CREATE INDEX IF NOT EXISTS "ix_jnuherb_disease_lookup_id" '
            'ON "jnuherb_disease_lookup"("disease_id")'
        )
        pd.DataFrame(meta_rows, columns=["table_name", "nrows"]).to_sql(
            "_meta", conn, index=False, if_exists="replace"
        )
        conn.commit()
    finally:
        conn.close()

    os.replace(tmp_path, db_path)
    return db_path


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=None,
        help="Directory containing the CSV tables "
        "(default: package data dir or JNUHERB_DATA_DIR)",
    )
    args = parser.parse_args(argv)
    out = build_db(args.data_dir)
    size_mb = out.stat().st_size / 1e6
    print(f"Built {out} ({size_mb:.1f} MB, {len(_TABLE_FILES)} tables)")


if __name__ == "__main__":
    main()
