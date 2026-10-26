"""Database query functions for jnuherb datasets.

Queries run against the bundled SQLite database (``jnuherb.db``) so that the
ID / name indexes do the filtering and only matching rows enter memory; the
multi-megabyte pandas table cache of the old whole-table-scan approach is
gone.  When the database is unavailable (CSV-only deployments, tests) every
getter transparently falls back to loading the full table, preserving the
historical behaviour.

Chinese disease names are served by a small bundled alias table
(``jnuherb_disease_aliases.csv``): common Chinese disease names map to the
exact vocabulary rows in the protein–disease link table, so reverse lookup
works offline, without any LLM translation step.

The underlying data tables are:

=============================  ============  ==================================
Table                          Records       Description
=============================  ============  ==================================
jnuherb_formula                ~6,200        TCM formulas with source literature
jnuherb_tcm                    ~4,700        Traditional Chinese medicines
jnuherb_chemicals              ~13,100       Chemical ingredients
jnuherb_proteins               ~14,200       Target proteins with gene names
jnuherb_formula_tcm_links      ~28,300       Formula-TCM connections
jnuherb_tcm_chemical_links     ~60,600       TCM-chemical connections
jnuherb_chemical_protein_links ~450,800      Chemical-protein connections
jnuherb_protein_disease_links  ~1,650,000    Protein–disease links (CTD / OT)
=============================  ============  ==================================
"""

from __future__ import annotations

import re
import sqlite3
from typing import Union

import pandas as pd

from jnuherb._data import _data_dir, _load_table, sqlite_connection


# Type alias for ``items`` parameter — a single string, list of strings, or a pandas Series.
ItemsLike = Union[str, list[str], "pd.Series[str]"]

_CJK_RE = re.compile(r"[一-鿿]")
_COL_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")

_DISEASE_EMPTY_COLUMNS = ["disease_id", "disease_name", "n_links",
                          "best_score", "sources"]

# Pre-aggregated disease directory built by build_db (one row per distinct
# disease).  Older databases built before it existed simply take the
# on-the-fly aggregation path below.
_LOOKUP_TABLE = "jnuherb_disease_lookup"


def _ensure_list(items: ItemsLike) -> list[str]:
    """Normalize ``items`` to a list of strings."""
    if isinstance(items, str):
        return [items]
    if isinstance(items, pd.Series):
        return items.tolist()
    return list(items)


def _sql_select(table: str, where: str | None = None,
                params: tuple = ()) -> pd.DataFrame | None:
    """Run a SELECT against jnuherb.db; None when the db file is missing.

    Only the missing-database case falls back to the pandas path; SQL
    errors propagate, because a silent fallback there would mask query
    bugs behind confusing whole-table results.

    All columns come back as strings with NULLs normalized to "" so callers
    see the same ``dtype=str, keep_default_na=False`` semantics as the CSV /
    pandas fallback path.
    """
    try:
        conn = sqlite_connection()
    except FileNotFoundError:
        return None
    try:
        sql = f'SELECT * FROM "{table}"'
        if where:
            sql += f" WHERE {where}"
        df = pd.read_sql_query(sql, conn, params=params, dtype=str)
    finally:
        conn.close()
    return df.fillna("")


def _check_col(col: str) -> str:
    """Guard against SQL injection via the ``by`` argument."""
    if not _COL_RE.match(col):
        raise ValueError(f"invalid column name: {col!r}")
    return col


def _in_clause(col: str, items: list[str], params: list) -> str:
    """Build a chunked ``col IN (...)`` clause, extending ``params``.

    Returns ``1=0`` for empty input so the clause never matches.
    """
    if not items:
        return "1=0"
    chunks = [items[i:i + 400] for i in range(0, len(items), 400)]
    parts = []
    for ch in chunks:
        parts.append(f'"{col}" IN ({",".join("?" * len(ch))})')
        params.extend(ch)
    return "(" + " OR ".join(parts) + ")"


def _like_escape(kw: str) -> str:
    """Escape LIKE wildcards in a user keyword."""
    return kw.replace("\\", "\\\\").replace("%", "\\%").replace("_", "\\_")


# ---------------------------------------------------------------------------
# Chinese disease aliases (offline reverse-lookup support)
# ---------------------------------------------------------------------------

_alias_cache: list[tuple[str, str, str]] | None = None


def _disease_aliases() -> list[tuple[str, str, str]]:
    """(cn_name, disease_id, disease_name) alias rows, lazily cached."""
    global _alias_cache
    if _alias_cache is None:
        rows: list[tuple[str, str, str]] = []
        try:
            path = _data_dir() / "jnuherb_disease_aliases.csv"
            if path.exists():
                adf = pd.read_csv(str(path), dtype=str, keep_default_na=False)
                for r in adf.itertuples(index=False):
                    rows.append((str(r[0]), str(r[1]), str(r[2])))
        except Exception:
            rows = []
        _alias_cache = rows
    return _alias_cache


def match_disease_aliases(keyword: str) -> list[str]:
    """Disease ids whose Chinese alias matches ``keyword``.

    Exact alias matches win outright; substring matches (in either
    direction) are only used when nothing matches exactly, so ``偏头痛``
    resolves to migraine rather than to the shorter alias ``头痛``
    (headache) it happens to contain.  Returns [] when nothing matches or
    the alias file is absent.
    """
    kw = keyword.strip()
    if not kw:
        return []
    exact: list[str] = []
    partial: list[str] = []
    for cn, did, _ in _disease_aliases():
        if did in exact or did in partial:
            continue
        if kw == cn:
            exact.append(did)
        elif kw in cn or cn in kw:
            partial.append(did)
    return exact or partial


# ---------------------------------------------------------------------------
# Getters
# ---------------------------------------------------------------------------

def get_formula(by: str, items: ItemsLike) -> pd.DataFrame:
    """Query formula information.

    Args:
        by: Column name to match against (e.g. ``"HVPID"``).
        items: Formula IDs to look up.

    Returns:
        DataFrame of matching formula records.
    """
    items_list = _ensure_list(items)
    params: list = []
    df = _sql_select("jnuherb_formula",
                     _in_clause(_check_col(by), items_list, params),
                     tuple(params))
    if df is None:
        df = _load_table("jnuherb_formula")
        df = df[df[by].isin(items_list)].copy()
    df = df.reset_index(drop=True)
    return df


def get_formula_tcm_links(by: str, items: ItemsLike) -> pd.DataFrame:
    """Query formula-TCM connection information.

    Args:
        by: Column name to match against (``"HVPID"`` or ``"HVMID"``).
        items: Formula or TCM IDs to look up.

    Returns:
        DataFrame of matching formula-TCM links.
    """
    items_list = _ensure_list(items)
    params: list = []
    df = _sql_select("jnuherb_formula_tcm_links",
                     _in_clause(_check_col(by), items_list, params),
                     tuple(params))
    if df is None:
        df = _load_table("jnuherb_formula_tcm_links")
        df = df[df[by].isin(items_list)].copy()
    df = df.reset_index(drop=True)
    return df


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
    params: list = []
    df = _sql_select("jnuherb_tcm",
                     _in_clause(_check_col(by), items_list, params),
                     tuple(params))
    if df is None:
        df = _load_table("jnuherb_tcm")
        df = df[df[by].isin(items_list)].copy()
    df = df.reset_index(drop=True)
    return df


def get_tcm_chem_links(by: str, items: ItemsLike) -> pd.DataFrame:
    """Query TCM-chemical connection information.

    Args:
        by: Column name to match against (``"HVMID"`` or ``"PubChem_id"``).
        items: TCM or chemical IDs to look up.

    Returns:
        DataFrame of matching TCM-chemical links.
    """
    items_list = _ensure_list(items)
    params: list = []
    df = _sql_select("jnuherb_tcm_chemical_links",
                     _in_clause(_check_col(by), items_list, params),
                     tuple(params))
    if df is None:
        df = _load_table("jnuherb_tcm_chemical_links")
        df = df[df[by].isin(items_list)].copy()
    df = df.reset_index(drop=True)
    return df


def get_chemicals(by: str, items: ItemsLike) -> pd.DataFrame:
    """Query chemical ingredient information.

    Results are deduplicated on ``PubChem_id`` so each chemical appears
    only once.

    Args:
        by: Column name to match against (``"PubChem_id"`` or ``"Name"``).
        items: Chemical IDs or names to look up.

    Returns:
        DataFrame of matching chemical records (unique by PubChem_id).
    """
    items_list = _ensure_list(items)
    params: list = []
    if by == "Name":
        # Case-insensitive matching for chemical names so that e.g.
        # "Deoxyshikonin" resolves to the database entry "deoxyshikonin".
        clause = f'LOWER("{_check_col(by)}") IN (%s)' % ",".join("?" * len(items_list))
        params = [str(i).lower() for i in items_list]
        if not items_list:
            clause = "1=0"
        df = _sql_select("jnuherb_chemicals", clause, tuple(params))
    else:
        df = _sql_select("jnuherb_chemicals",
                         _in_clause(_check_col(by), items_list, params),
                         tuple(params))
    if df is None:
        df = _load_table("jnuherb_chemicals")
        if by == "Name":
            items_lower = {str(i).lower() for i in items_list}
            mask = df[by].str.lower().isin(items_lower)
        else:
            mask = df[by].isin(items_list)
        df = df[mask].copy()
    df = df.drop_duplicates(subset=["PubChem_id"])
    df = df.reset_index(drop=True)
    return df


def get_chem_protein_links(
    by: str, items: ItemsLike, score: float = 900
) -> pd.DataFrame:
    """Query chemical-protein connection information.

    Filters links where ``Combined_score >= score`` and converts
    ``Combined_score`` from the stored 0--1000 integer range to
    a 0--1 float.

    Args:
        by: Column name to match against (``"PubChem_id"`` or ``"Ensembl_ID"``).
        items: Chemical or protein IDs to look up.
        score: Minimum combined_score threshold (default 900,
            corresponding to 0.9 after conversion).

    Returns:
        DataFrame of matching chemical-protein links with
        ``Combined_score`` in [0, 1].
    """
    items_list = _ensure_list(items)
    params: list = []
    clause = _in_clause(_check_col(by), items_list, params)
    clause += ' AND CAST("Combined_score" AS REAL) >= ?'
    params.append(score)
    df = _sql_select("jnuherb_chemical_protein_links", clause, tuple(params))
    if df is None:
        df = _load_table("jnuherb_chemical_protein_links")
        mask = df[by].isin(items_list) & (pd.to_numeric(df["Combined_score"]) >= score)
        df = df[mask].copy()
    df["Combined_score"] = pd.to_numeric(df["Combined_score"]) / 1000.0
    df = df.reset_index(drop=True)
    return df


def get_protein_disease_links(
    by_disease: ItemsLike | None = None,
    by_protein: ItemsLike | None = None,
    min_score: float = 0.0,
    sources: ItemsLike | None = None,
) -> pd.DataFrame:
    """Query protein–disease association links (CTD / Open Targets).

    At least one of ``by_disease`` / ``by_protein`` should be given; when
    both are None the whole (filtered) table is returned.

    Args:
        by_disease: disease_ids (``"EFO:0000180"``) or disease_names to match.
        by_protein: Ensembl protein IDs to match.
        min_score: keep rows with numeric ``score`` >= this; rows with an
            empty score (CTD curated) always pass.
        sources: restrict to sources (e.g. ``["OT", "CTD"]``).

    Returns:
        DataFrame of matching links, score as float in [0, 1].
    """
    clauses: list[str] = []
    params: list = []
    if by_disease is not None:
        items = _ensure_list(by_disease)
        # disease_id matches exactly (prefixed vocabulary ids); disease_name
        # matches case-insensitively ("Hypertension" vs "hypertension").
        # The COLLATE NOCASE index keeps this off the full table scan.
        id_clause = _in_clause("disease_id", items, params)
        name_params_start = len(params)
        name_clause = 'disease_name COLLATE NOCASE IN (%s)' % ",".join("?" * len(items)) if items else "1=0"
        params[name_params_start:name_params_start] = [str(i) for i in items]
        clauses.append(f"({id_clause} OR {name_clause})")
    if by_protein is not None:
        clauses.append(_in_clause("Ensembl_ID", _ensure_list(by_protein), params))
    if sources is not None:
        clauses.append(_in_clause("source", _ensure_list(sources), params))
    if min_score > 0:
        clauses.append("((score != '' AND CAST(score AS REAL) >= ?) OR score = '')")
        params.append(min_score)
    where = " AND ".join(clauses) if clauses else None
    df = _sql_select("jnuherb_protein_disease_links", where, tuple(params))
    if df is None:
        df = _pandas_protein_disease_links(by_disease, by_protein, min_score, sources)
    df["score"] = pd.to_numeric(df["score"], errors="coerce")
    df = df.reset_index(drop=True)
    return df


def _pandas_protein_disease_links(
    by_disease: ItemsLike | None,
    by_protein: ItemsLike | None,
    min_score: float,
    sources: ItemsLike | None,
) -> pd.DataFrame:
    """Whole-table pandas fallback (CSV deployments / tests)."""
    df = _load_table("jnuherb_protein_disease_links")
    if by_disease is not None:
        items = _ensure_list(by_disease)
        names_lower = df["disease_name"].str.lower().isin(
            {str(i).lower() for i in items})
        df = df[df["disease_id"].isin(items) | names_lower]
    if by_protein is not None:
        df = df[df["Ensembl_ID"].isin(_ensure_list(by_protein))]
    if sources is not None:
        df = df[df["source"].isin(_ensure_list(sources))]
    if min_score > 0:
        num = pd.to_numeric(df["score"], errors="coerce")
        df = df[(num >= min_score) | (df["score"] == "")]
    return df.copy()


def _disease_stats_rows(disease_ids: list[str]) -> pd.DataFrame:
    """Aggregate n_links / best_score / sources for the given disease ids.

    Uses the pre-aggregated lookup table when the database was built by a
    recent build_db; otherwise aggregates the links table on the fly
    (slower, but correct for databases built before the lookup existed).
    """
    params: list = []
    clause = _in_clause("disease_id", disease_ids, params)
    conn = None
    try:
        conn = sqlite_connection()
        if _has_table(conn, _LOOKUP_TABLE):
            sql = (f'SELECT disease_id, disease_name, n_links, best_score, '
                   f'sources FROM "{_LOOKUP_TABLE}" WHERE {clause} '
                   'ORDER BY n_links DESC')
        else:
            sql = (
                'SELECT disease_id, disease_name, '
                'COUNT(DISTINCT "Ensembl_ID") AS n_links, '
                "MAX(CAST(NULLIF(score, '') AS REAL)) AS best_score, "
                "GROUP_CONCAT(DISTINCT source) AS sources "
                'FROM "jnuherb_protein_disease_links" '
                f"WHERE {clause} "
                "GROUP BY disease_id, disease_name "
                "ORDER BY n_links DESC"
            )
        df = pd.read_sql_query(sql, conn, params=tuple(params), dtype=str)
    except FileNotFoundError:
        df = None
    finally:
        if conn is not None:
            conn.close()
    if df is None:
        # CSV fallback: restrict the pandas search to the alias ids.
        df = _pandas_disease_search("", 10_000)
        df = df[df["disease_id"].isin(disease_ids)]
    return df.fillna("")


def _has_table(conn, name: str) -> bool:
    row = conn.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?",
        (name,),
    ).fetchone()
    return row is not None


def get_diseases(keyword: str, limit: int = 20) -> pd.DataFrame:
    """Fuzzy-search disease names/ids for the Reverse disease lookup box.

    English keywords match as substring (case-insensitive) against
    ``disease_name`` and prefix against ``disease_id``.  Chinese keywords
    first go through the bundled alias table (offline, no LLM); only
    aliases that actually exist in the link table are returned, ranked by
    number of supporting protein links.
    """
    kw = keyword.strip()
    if not kw:
        return pd.DataFrame(columns=_DISEASE_EMPTY_COLUMNS)

    # Chinese query: alias table drives the result.
    if _CJK_RE.search(kw):
        ids = match_disease_aliases(kw)
        if not ids:
            return pd.DataFrame(columns=_DISEASE_EMPTY_COLUMNS)
        agg = _disease_stats_rows(ids)
        if agg.empty:
            return pd.DataFrame(columns=_DISEASE_EMPTY_COLUMNS)
        return _finalize_disease_suggestions(agg, limit)

    like_kw = _like_escape(kw)
    clauses = [
        "disease_name LIKE ? ESCAPE '\\'",
        "disease_id LIKE ? ESCAPE '\\'",
    ]
    params = [f"%{like_kw}%", f"{like_kw}%"]
    conn = None
    try:
        conn = sqlite_connection()
        if _has_table(conn, _LOOKUP_TABLE):
            # Pre-aggregated directory: substring scan over ~17k rows.
            sql = (
                f'SELECT disease_id, disease_name, n_links, best_score, sources '
                f'FROM "{_LOOKUP_TABLE}" '
                f"WHERE {' OR '.join(clauses)} "
                "ORDER BY n_links DESC LIMIT ?"
            )
        else:
            # Older database: aggregate the raw links table directly.
            sql = (
                'SELECT disease_id, disease_name, '
                'COUNT(DISTINCT "Ensembl_ID") AS n_links, '
                "MAX(CAST(NULLIF(score, '') AS REAL)) AS best_score, "
                "GROUP_CONCAT(DISTINCT source) AS sources "
                'FROM "jnuherb_protein_disease_links" '
                f"WHERE {' OR '.join(clauses)} "
                "GROUP BY disease_id, disease_name "
                "ORDER BY n_links DESC LIMIT ?"
            )
        params.append(int(limit))
        df = pd.read_sql_query(sql, conn, params=tuple(params), dtype=str)
    except FileNotFoundError:
        df = None
    finally:
        if conn is not None:
            conn.close()
    if df is None:
        df = _pandas_disease_search(kw, limit)
        return df
    return _finalize_disease_suggestions(df, limit)


def _finalize_disease_suggestions(agg: pd.DataFrame, limit: int) -> pd.DataFrame:
    """Normalize the SQL aggregation output to the public column semantics."""
    agg = agg.copy()
    agg["n_links"] = pd.to_numeric(agg["n_links"]).astype(int)
    agg["best_score"] = pd.to_numeric(agg["best_score"], errors="coerce")
    # GROUP_CONCAT has no ordering guarantee — restore the sorted a/b form.
    agg["sources"] = agg["sources"].map(
        lambda s: "/".join(sorted(set(str(s).split(",")))) if s else "")
    agg = agg.sort_values("n_links", ascending=False)
    return agg.head(limit).reset_index(drop=True)


def _pandas_disease_search(keyword: str, limit: int) -> pd.DataFrame:
    """Pandas fallback for get_diseases (CSV deployments / tests)."""
    kw = keyword.lower()
    df = _load_table("jnuherb_protein_disease_links")
    mask = (df["disease_name"].str.lower().str.contains(kw, na=False, regex=False)
            | df["disease_id"].str.lower().str.startswith(kw))
    sub = df[mask]
    if sub.empty:
        return pd.DataFrame(columns=_DISEASE_EMPTY_COLUMNS)
    agg = (sub.groupby(["disease_id", "disease_name"], as_index=False)
              .agg(n_links=("Ensembl_ID", "nunique"),
                   best_score=("score", "max"),
                   sources=("source", lambda s: "/".join(sorted(set(s))))))
    agg["best_score"] = pd.to_numeric(agg["best_score"], errors="coerce")
    return agg.sort_values("n_links", ascending=False).head(limit)


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
    params: list = []
    df = _sql_select("jnuherb_proteins",
                     _in_clause(_check_col(by), items_list, params),
                     tuple(params))
    if df is None:
        df = _load_table("jnuherb_proteins")
        df = df[df[by].isin(items_list)].copy()
    df = df.drop_duplicates(subset=["Ensembl_ID"])
    df = df.reset_index(drop=True)
    return df
