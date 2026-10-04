"""HerbiV scoring and combination optimization.

Provides:

* ``herbiv_score`` — probability-union-based scoring for formulas,
  TCMs, and chemicals against target proteins.
* ``component`` — knapsack-based optimizer that finds the best
  TCM or formula combinations.
* ``knapsack`` — internal DP solver with probability-union recurrence.
"""

from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# herbiv_score
# ---------------------------------------------------------------------------


def herbiv_score(
    tcm: pd.DataFrame,
    tcm_chem_links: pd.DataFrame,
    chem: pd.DataFrame,
    chem_protein_links: pd.DataFrame,
    formula: Optional[pd.DataFrame] = None,
    formula_tcm_links: Optional[pd.DataFrame] = None,
    weights: Optional[dict[str, float]] = None,
) -> dict[str, Optional[pd.DataFrame]]:
    """Calculate HerbiV Scores for formulas, TCMs, and chemicals.

    For each target protein, a probability-union score is computed
    at three levels:

    * **Chemical score**: :math:`1 - \\prod (1 - s_i)` over all
      chemical-protein link scores for that protein.
    * **TCM score**: union of chemical scores for all chemicals in
      the TCM.
    * **Formula score**: union of TCM scores for all TCMs in the
      formula.

    An ``Importance Score`` column is then computed as the
    (possibly weighted) mean across all protein score columns.

    Args:
        tcm: TCM information (must have ``HVMID`` column).
        tcm_chem_links: TCM-chemical links (``HVMID``, ``HVCID``).
        chem: Chemical information (must have ``HVCID`` column).
        chem_protein_links: Chemical-protein links with
            ``Combined_score`` in [0, 1] (must have ``HVCID``
            and ``Ensembl_ID`` columns).
        formula: Optional formula information (``HVPID`` column).
        formula_tcm_links: Optional formula-TCM links.
        weights: Optional per-protein weights.  If None equal
            weights are used.

    Returns:
        Dict with keys ``tcm_and_score``, ``chem_and_score``,
        ``formula_and_score`` (the last is None if formula not
        provided).
    """
    formula_and_score: Optional[pd.DataFrame] = None
    if formula is not None:
        formula_and_score = formula.copy()
    tcm_and_score = tcm.copy()
    chem_and_score = chem.copy()

    proteins_id = list(chem_protein_links["Ensembl_ID"].unique())

    # --- pre-group for efficiency ---
    cpl_grouped = chem_protein_links.groupby("HVCID")
    tcl_grouped = tcm_chem_links.groupby("HVMID")["HVCID"].apply(set).to_dict()

    fml_tcm_map: dict[str, set] = {}
    if formula is not None and formula_tcm_links is not None:
        fml_tcm_map = (
            formula_tcm_links.groupby("HVPID")["HVMID"].apply(set).to_dict()
        )

    # --- per-protein HerbiV Score columns (collect, then assign at once) ---
    chem_score_cols: dict[str, list[float]] = {}
    tcm_score_cols: dict[str, list[float]] = {}
    fml_score_cols: dict[str, list[float]] = {}

    for protein in proteins_id:
        score_col = f"{protein} HerbiV Score"

        # Chemical score
        chem_scores = []
        for hvcid in chem_and_score["HVCID"]:
            grp = cpl_grouped.get_group(hvcid) if hvcid in cpl_grouped.groups else None
            if grp is not None:
                subset = grp.loc[grp["Ensembl_ID"] == protein, "Combined_score"]
            else:
                subset = pd.Series([], dtype=float)
            chem_scores.append(_union_score(subset))
        chem_score_cols[score_col] = chem_scores

        # Cache chemical scores for TCM lookup
        chem_score_series = pd.Series(chem_scores, index=chem_and_score["HVCID"])

        # TCM score
        tcm_scores = []
        for hvm in tcm_and_score["HVMID"]:
            cids = tcl_grouped.get(hvm, set())
            subset = chem_score_series.loc[chem_score_series.index.isin(cids)]
            tcm_scores.append(_union_score(subset))
        tcm_score_cols[score_col] = tcm_scores

        # Formula score
        if formula is not None and formula_and_score is not None:
            fml_scores = []
            tcm_score_series = pd.Series(tcm_scores, index=tcm_and_score["HVMID"])
            for hvpid in formula_and_score["HVPID"]:
                mids = fml_tcm_map.get(hvpid, set())
                subset = tcm_score_series.loc[tcm_score_series.index.isin(mids)]
                fml_scores.append(_union_score(subset))
            fml_score_cols[score_col] = fml_scores

    # --- Assign all score columns at once (avoids fragmentation) ---
    if chem_score_cols:
        chem_new = pd.DataFrame(chem_score_cols, index=chem_and_score.index)
        chem_and_score = pd.concat(
            [chem_and_score.drop(columns=list(chem_score_cols.keys()), errors="ignore"),
             chem_new], axis=1
        )
    if tcm_score_cols:
        tcm_new = pd.DataFrame(tcm_score_cols, index=tcm_and_score.index)
        tcm_and_score = pd.concat(
            [tcm_and_score.drop(columns=list(tcm_score_cols.keys()), errors="ignore"),
             tcm_new], axis=1
        )
    if fml_score_cols and formula_and_score is not None:
        fml_new = pd.DataFrame(fml_score_cols, index=formula_and_score.index)
        formula_and_score = pd.concat(
            [formula_and_score.drop(columns=list(fml_score_cols.keys()), errors="ignore"),
             fml_new], axis=1
        )

    # --- weighted Importance Score ---
    weight_cols = [f"{p} HerbiV Score" for p in proteins_id]
    if weights is None:
        w = np.ones(len(weight_cols))
    else:
        w = np.array([weights.get(c, 1.0) for c in weight_cols])

    if formula is not None and formula_and_score is not None:
        formula_and_score["Importance Score"] = (
            formula_and_score[weight_cols].astype(float).mul(w).mean(axis=1)
        )
        formula_and_score = formula_and_score.sort_values(
            "Importance Score", ascending=False
        )
        formula_and_score.reset_index(drop=True, inplace=True)

    tcm_and_score["Importance Score"] = (
        tcm_and_score[weight_cols].astype(float).mul(w).mean(axis=1)
    )
    tcm_and_score = tcm_and_score.sort_values(
        "Importance Score", ascending=False
    )
    tcm_and_score.reset_index(drop=True, inplace=True)

    chem_and_score["Importance Score"] = (
        chem_and_score[weight_cols].astype(float).mul(w).mean(axis=1)
    )
    chem_and_score = chem_and_score.sort_values(
        "Importance Score", ascending=False
    )
    chem_and_score.reset_index(drop=True, inplace=True)

    return {
        "tcm_and_score": tcm_and_score,
        "chem_and_score": chem_and_score,
        "formula_and_score": formula_and_score,
    }


def _union_score(scores: pd.Series) -> float:
    """Probability union of a series of scores: 1 - prod(1 - scores)."""
    s = scores.astype(float).values
    if len(s) == 0:
        return 0.0
    return float(1.0 - np.prod(1.0 - s))


# ---------------------------------------------------------------------------
# component
# ---------------------------------------------------------------------------


def component(
    items_and_score: pd.DataFrame,
    random_state: Optional[int] = None,
    num: int = 1000,
    c: int = 10,
) -> pd.DataFrame:
    """Optimize TCM or formula combinations using knapsack search.

    Repeatedly samples random subsets and applies the knapsack DP
    to find high-scoring combinations.

    Args:
        items_and_score: TCM or formula data with an
            ``Importance Score`` column.
        random_state: Optional random seed for reproducibility.
        num: Number of random restarts (default 1000).
        c: Knapsack capacity parameter (default 10).

    Returns:
        DataFrame with columns ``items``, ``Importance Score``,
        ``Boost``, sorted by Boost descending.
    """
    # Determine ID column
    if "HVPID" in items_and_score.columns:
        by = "HVPID"
    else:
        by = "HVMID"

    rng = np.random.RandomState(random_state)

    n_items = len(items_and_score)
    weights_base = np.ones(n_items, dtype=int)
    names_vec = items_and_score[by].values
    values = items_and_score["Importance Score"].astype(float).values
    # Cap subset size to keep knapsack DP fast even with many candidates
    n_select = min(max(1, int(np.ceil(n_items / 10))), 30)

    dp_scores: list[float] = []
    items_ls: list[list[str]] = []
    forbidden_sets: list[frozenset[str]] = []

    for i in range(num):
        indices = rng.permutation(n_items)[:n_select]
        cur_weights = weights_base[indices]
        cur_names = names_vec[indices]
        cur_values = values[indices]

        result = _knapsack(cur_weights, cur_names, cur_values, forbidden_sets, c)
        if result is not None:
            dp_scores.append(result["dp"])
            items_ls.append(result["items"])
            forbidden_sets.append(frozenset(result["items"]))

    components_df = pd.DataFrame(
        {
            "items": [";".join(it) for it in items_ls],
            "Importance Score": dp_scores,
        }
    )

    # Compute Boost
    def _compute_boost(row: pd.Series) -> float:
        item_names = row["items"].split(";")
        mask = items_and_score[by].isin(item_names)
        if not mask.any():
            return 0.0
        max_individual = items_and_score.loc[mask, "Importance Score"].max()
        if max_individual == 0:
            return 0.0
        return (float(row["Importance Score"]) - max_individual) / max_individual

    components_df["Boost"] = components_df.apply(_compute_boost, axis=1)
    components_df = components_df.sort_values("Boost", ascending=False)
    components_df.reset_index(drop=True, inplace=True)

    return components_df


# ---------------------------------------------------------------------------
# _knapsack  (internal)
# ---------------------------------------------------------------------------


def _knapsack(
    weights: np.ndarray,
    names_vec: np.ndarray,
    values: np.ndarray,
    forbidden_combinations: list[frozenset[str]],
    c: int = 10,
) -> Optional[dict]:
    """Knapsack DP with probability-union recurrence.

    Finds a high-scoring subset of items using a dynamic programming
    approach where the objective is the probability union
    :math:`1 - \\prod(1 - v_i)` rather than the usual sum.

    Args:
        weights: Item weights (all 1 for equal weighting).
        names_vec: Item names/IDs.
        values: Item Importance Scores.
        forbidden_combinations: Previously found combinations to avoid,
            represented as frozensets for fast intersection checks.
        c: Capacity parameter.

    Returns:
        Dict with ``dp`` (best score) and ``items`` (list of selected
        names), or None if no valid solution is found.
    """
    n = len(weights)
    capacity = c + 1  # DP table columns (0..c)

    # Build name -> forbidden sets containing that name for fast lookup
    name_to_forbidden: dict[str, list[frozenset[str]]] = {}
    for forbidden in forbidden_combinations:
        for name in forbidden:
            name_to_forbidden.setdefault(name, []).append(forbidden)

    dp = np.zeros((n + 1, capacity))
    items = np.empty((n + 1, capacity), dtype=object)
    for i in range(n + 1):
        for j in range(capacity):
            items[i, j] = ""

    for i in range(n):
        name_i = names_vec[i]
        relevant_forbidden = name_to_forbidden.get(name_i, [])

        for j in range(1, capacity):
            if weights[i] <= j:
                # Check for conflicts with relevant forbidden combinations
                conflict = False
                if relevant_forbidden:
                    prev_items = items[i, j - weights[i]]
                    if prev_items:
                        prev_set = set(prev_items.split(", "))
                        for forbidden in relevant_forbidden:
                            if prev_set & forbidden:
                                conflict = True
                                break

                if conflict:
                    dp[i + 1, j] = dp[i, j]
                    items[i + 1, j] = items[i, j]
                else:
                    prev_val = dp[i, j - weights[i]]
                    new_val = 1.0 - (1.0 - values[i]) * (1.0 - prev_val)
                    if new_val > dp[i, j]:
                        dp[i + 1, j] = new_val
                        prev_str = items[i, j - weights[i]]
                        if prev_str:
                            items[i + 1, j] = f"{name_i}, {prev_str}"
                        else:
                            items[i + 1, j] = name_i
                    else:
                        dp[i + 1, j] = dp[i, j]
                        items[i + 1, j] = items[i, j]
            else:
                dp[i + 1, j] = dp[i, j]
                items[i + 1, j] = items[i, j]

    # Clean up trailing ", "
    last_row_items = []
    for j in range(capacity):
        s = str(items[n, j]).strip()
        if s.endswith(", "):
            s = s[:-2]
        last_row_items.append(s)

    # MLE-based elbow detection for number of components
    dp_row = dp[n, :]
    row_sum = dp_row.sum()
    if row_sum == 0:
        return None

    score_ratio = np.cumsum(dp_row) / row_sum
    m = len(score_ratio)
    mle_estimates = (score_ratio - 1.0 / m) / np.sqrt(2.0 / m)

    num_components = int(np.argmin(mle_estimates))
    if num_components <= 1:
        num_components = 2

    selected_items_str = last_row_items[num_components]
    selected_items = [s for s in selected_items_str.split(", ") if s]

    return {"dp": float(dp[n, num_components]), "items": selected_items}
