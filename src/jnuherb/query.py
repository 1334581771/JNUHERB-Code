"""Database query functions for jnuherb datasets.

These functions read from the bundled CSV data files and return filtered
DataFrames.  All seven exported functions follow the same pattern:

1. Load the full CSV table (cached after first access).
2. Filter rows where ``by`` column matches any value in ``items``.
3. Optionally deduplicate on a key column.
4. Return the result with a fresh index.

The underlying data tables are:

=============================  ========  ==================================
Table                          Records   Description
=============================  ========  ==================================
jnuherb_formula                ~6,200    TCM formulas with source literature
jnuherb_tcm                    ~4,700    Traditional Chinese medicines
jnuherb_chemicals              ~13,100   Chemical ingredients
jnuherb_proteins               ~14,200   Target proteins with gene names
jnuherb_formula_tcm_links      ~28,300   Formula-TCM connections
jnuherb_tcm_chemical_links     ~60,600   TCM-chemical connections
jnuherb_chemical_protein_links ~450,800  Chemical-protein connections
=============================  ========  ==================================
"""

from __future__ import annotations

from typing import Union

import pandas as pd

from jnuherb._data import _load_table


# Type alias for ``items`` parameter — a single string, list of strings, or a pandas Series.
ItemsLike = Union[str, list[str], "pd.Series[str]"]


def _ensure_list(items: ItemsLike) -> list[str]:
    """Normalize ``items`` to a list of strings."""
    if isinstance(items, str):
        return [items]
    if isinstance(items, pd.Series):
        return items.tolist()
    return list(items)


def get_formula(by: str, items: ItemsLike) -> pd.DataFrame:
    """Query formula information.

    Args:
        by: Column name to match against (e.g. ``"HVPID"``).
        items: Formula IDs to look up.

    Returns:
        DataFrame of matching formula records.
    """
    items_list = _ensure_list(items)
    df = _load_table("jnuherb_formula")
    result = df[df[by].isin(items_list)].copy()
    result.reset_index(drop=True, inplace=True)
    return result


def get_formula_tcm_links(by: str, items: ItemsLike) -> pd.DataFrame:
    """Query formula-TCM connection information.

    Args:
        by: Column name to match against (``"HVPID"`` or ``"HVMID"``).
        items: Formula or TCM IDs to look up.

    Returns:
        DataFrame of matching formula-TCM links.
    """
    items_list = _ensure_list(items)
    df = _load_table("jnuherb_formula_tcm_links")
    result = df[df[by].isin(items_list)].copy()
    result.reset_index(drop=True, inplace=True)
    return result


def get_tcm(by: str, items: ItemsLike) -> pd.DataFrame:
    """Query TCM (Traditional Chinese Medicine) information.

    Args:
        by: Column name to match against.  Supported values:
            ``"HVMID"``, ``"cn_name"``, ``"pinyin_name"``,
            ``"en_name"``, ``"latin_name"``.
        items: TCM IDs or names to look up.

    Returns:
        DataFrame of matching TCM records.
    """
    items_list = _ensure_list(items)
    df = _load_table("jnuherb_tcm")
    result = df[df[by].isin(items_list)].copy()
    result.reset_index(drop=True, inplace=True)
    return result


def get_tcm_chem_links(by: str, items: ItemsLike) -> pd.DataFrame:
    """Query TCM-chemical connection information.

    Args:
        by: Column name to match against (``"HVMID"`` or ``"HVCID"``).
        items: TCM or chemical IDs to look up.

    Returns:
        DataFrame of matching TCM-chemical links.
    """
    items_list = _ensure_list(items)
    df = _load_table("jnuherb_tcm_chemical_links")
    result = df[df[by].isin(items_list)].copy()
    result.reset_index(drop=True, inplace=True)
    return result


def get_chemicals(by: str, items: ItemsLike) -> pd.DataFrame:
    """Query chemical ingredient information.

    Results are deduplicated on ``HVCID`` so each chemical appears
    only once.

    Args:
        by: Column name to match against (``"HVCID"`` or ``"Name"``).
        items: Chemical IDs or names to look up.

    Returns:
        DataFrame of matching chemical records (unique by HVCID).
    """
    items_list = _ensure_list(items)
    df = _load_table("jnuherb_chemicals")
    result = df[df[by].isin(items_list)].copy()
    result = result.drop_duplicates(subset=["HVCID"])
    result.reset_index(drop=True, inplace=True)
    return result


def get_chem_protein_links(
    by: str, items: ItemsLike, score: float = 900
) -> pd.DataFrame:
    """Query chemical-protein connection information.

    Filters links where ``Combined_score >= score`` and converts
    ``Combined_score`` from the stored 0--1000 integer range to
    a 0--1 float.

    Args:
        by: Column name to match against (``"HVCID"`` or ``"Ensembl_ID"``).
        items: Chemical or protein IDs to look up.
        score: Minimum combined_score threshold (default 900,
            corresponding to 0.9 after conversion).

    Returns:
        DataFrame of matching chemical-protein links with
        ``Combined_score`` in [0, 1].
    """
    items_list = _ensure_list(items)
    df = _load_table("jnuherb_chemical_protein_links")
    mask = df[by].isin(items_list) & (pd.to_numeric(df["Combined_score"]) >= score)
    result = df[mask].copy()
    result["Combined_score"] = pd.to_numeric(result["Combined_score"]) / 1000.0
    result.reset_index(drop=True, inplace=True)
    return result


def get_proteins(by: str, items: ItemsLike) -> pd.DataFrame:
    """Query protein target information.

    Results are deduplicated on ``Ensembl_ID`` so each protein
    appears only once.

    Args:
        by: Column name to match against (``"Ensembl_ID"`` or ``"gene_name"``).
        items: Protein IDs or gene names to look up.

    Returns:
        DataFrame of matching protein records (unique by Ensembl_ID).
    """
    items_list = _ensure_list(items)
    df = _load_table("jnuherb_proteins")
    result = df[df[by].isin(items_list)].copy()
    result = result.drop_duplicates(subset=["Ensembl_ID"])
    result.reset_index(drop=True, inplace=True)
    return result
