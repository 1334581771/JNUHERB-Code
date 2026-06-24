"""Internal data loading with lazy caching.

Reads bundled CSV files from the package data directory using
importlib.resources, with in-memory caching so each table is
parsed at most once per session.
"""

from __future__ import annotations

from importlib.resources import files

import pandas as pd

_cache: dict[str, pd.DataFrame] = {}

_TABLE_FILES: dict[str, str] = {
    "jnuherb_formula": "jnuherb_formula.csv",
    "jnuherb_tcm": "jnuherb_tcm.csv",
    "jnuherb_chemicals": "jnuherb_chemicals.csv",
    "jnuherb_proteins": "jnuherb_proteins.csv",
    "jnuherb_formula_tcm_links": "jnuherb_formula_tcm_links.csv",
    "jnuherb_tcm_chemical_links": "jnuherb_tcm_chemical_links.csv",
    "jnuherb_chemical_protein_links": "jnuherb_chemical_protein_links.csv",
}


def _load_table(name: str) -> pd.DataFrame:
    """Lazy-load a bundled CSV table with in-memory caching.

    Args:
        name: Logical table name (key of _TABLE_FILES).

    Returns:
        DataFrame with all columns as strings.
    """
    if name not in _cache:
        path = files("jnuherb.data").joinpath(_TABLE_FILES[name])
        _cache[name] = pd.read_csv(str(path), dtype=str, keep_default_na=False)
    return _cache[name]


def _get_data_path(filename: str) -> str:
    """Get the absolute filesystem path to a bundled data file.

    Args:
        filename: CSV file name, e.g. ``"jnuherb_tcm.csv"``.

    Returns:
        Absolute path as a string.
    """
    return str(files("jnuherb.data").joinpath(filename))
