#!/usr/bin/env python
"""JNUHERB — Interactive Network Pharmacology Analysis (Streamlit GUI).

Supports two analysis modes:

* **Forward**: TCMs / formulas / monomers → chemicals → proteins
* **Reverse**: disease / target proteins → chemicals → TCMs → formulas

Launch via::

    python -m jnuherb
    jnuherb
    streamlit run src/jnuherb/data/streamlit_app.py
"""

from __future__ import annotations

import base64
import concurrent.futures
import fcntl
import json
import os
import tempfile
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import re
import streamlit as st

matplotlib.use("Agg")

from jnuherb.analysis import from_proteins, from_tcm_or_formula
from jnuherb.enrichment import extract_genes_from_result, run_enrichment
from jnuherb.network import (
    calculate_degree,
    degree_summary,
    filter_formula_network,
    plot_degree_barchart,
    plot_degree_boxplot,
)

# ---------------------------------------------------------------------------
# Page config
# ---------------------------------------------------------------------------
try:
    from importlib.resources import files as _res_files
    _BRAND_ICON = _res_files("jnuherb.data.assets").joinpath(
        "jnuherb-brand-icon.png").read_bytes()
except Exception:  # noqa: BLE001 - fall back to the emoji favicon
    _BRAND_ICON = "🌿"

st.set_page_config(
    page_title="JNUHERB",
    page_icon=_BRAND_ICON,
    layout="wide",
    initial_sidebar_state="expanded",
)
# Hide Streamlit default chrome (Deploy button, "Made with Streamlit" footer, etc.)
# plus the JNUHERB brand styles (hero banner, entry cards).
# NOTE: do NOT hide `header` because it contains the sidebar toggle button.
_HIDE_STYLES = """
<style>
/* ---- Design tokens (ui-redesign-proposal/DESIGN-DIRECTION.md) ----
   secondary text #65796E · divider #DFE7E0 · cards radius 12px,
   spacing scale 8/12/16/24/32px */
#MainMenu {visibility: hidden;}
footer {visibility: hidden;}
.stDeployButton {display: none !important;}
hr { border: none; border-top: 1px solid #DFE7E0; margin: 1.6rem 0; }

/* ---- Brand: hero banner ---- */
.jnu-hero {
    background: linear-gradient(120deg, #0f4a37 0%, #16634a 55%, #1e7a5a 100%);
    border-radius: 16px;
    padding: 2.4rem 2.6rem 1.6rem;
    color: #ffffff;
    margin-bottom: 1.2rem;
}
.jnu-hero-kicker {
    font-size: 1.0rem;
    font-weight: 700;
    letter-spacing: 0.30em;
    opacity: 0.85;
    text-transform: uppercase;
}
.jnu-kicker-icon {
    height: 1.25em;
    vertical-align: -0.22em;
    margin-right: 0.4em;
}
.jnu-hero-title {
    font-size: 2.1rem;
    font-weight: 800;
    letter-spacing: 0.01em;
    line-height: 1.15;
    margin-top: 0.7rem;
}
.jnu-hero-sub {
    font-size: 1.0rem;
    opacity: 0.88;
    margin-top: 0.6rem;
    line-height: 1.5;
}
.jnu-hero-icons {
    display: flex;
    justify-content: space-between;
    align-items: stretch;
    gap: 0.4rem;
    margin-top: 1.8rem;
}
.jnu-hero-ent {
    text-align: center;
    flex: 1;
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: flex-end;
}
.jnu-hero-ent img {
    height: 104px;
    max-width: 100%;
    object-fit: contain;
    filter: drop-shadow(0 6px 14px rgba(0,0,0,0.28));
}
.jnu-hero-ent-num {
    margin-top: 0.55rem;
    font-size: 1.4rem;
    font-weight: 800;
    letter-spacing: 0.02em;
}
.jnu-hero-ent-lbl {
    font-size: 1.0rem;
    font-weight: 500;
    letter-spacing: 0.04em;
    opacity: 0.85;
    margin-top: 0.1rem;
}

/* ---- Workflow section + cards ---- */
.jnu-workflow-h2 {
    font-size: 1.0rem;
    font-weight: 800;
    color: #173e34;
    margin: 1.2rem 0 0.2rem;
}
.jnu-workflow-sub {
    font-size: 1.0rem;
    color: #65796E;
    margin-bottom: 0.9rem;
}
@media (max-width: 760px) {
    .jnu-card-grid { grid-template-columns: 1fr; }
    .jnu-hero-ent img { height: 60px; }
    .jnu-hero-ent-num { font-size: 1.15rem; }
    .jnu-hero-ent-lbl { font-size: 0.95rem; }
}

/* ---- Verified examples band ---- */

/* ---- Results provenance bar (which submitted run the tabs belong to) ---- */
.jnu-runbar {
    display: flex;
    align-items: center;
    gap: 0.7rem;
    flex-wrap: wrap;
    background: #ffffff;
    border: 1px solid #DFE7E0;
    border-left: 4px solid #16634a;
    border-radius: 10px;
    padding: 0.5rem 0.9rem;
    margin: 0.4rem 0 0.8rem;
    font-size: 1.0rem;
    color: #65796E;
}
.jnu-runbar .jnu-runbar-mode {
    font-weight: 700;
    color: #213B32;
}
.jnu-runbar .jnu-runbar-id {
    font-weight: 700;
    color: #16634a;
    font-family: monospace;
}
.jnu-runbar-stale {
    background: #FFF8E6;
    border-color: #E8C96A;
    border-left-color: #C9962E;
    color: #7A5B00;
}
.jnu-runbar-stale .jnu-runbar-id { color: #7A5B00; }
.jnu-runbar-stale .jnu-runbar-mode { color: #7A5B00; }
/* ---- Flow chip inside native workflow cards ---- */
.jnu-flow-chip {
    display: inline-block;
    font-size: 1.0rem;
    font-weight: 700;
    color: #16634a;
    background: #EEF3EE;
    border: 1px solid #DFE7E0;
    border-radius: 8px;
    padding: 0.15rem 0.6rem;
    margin: 0.1rem 0 0.3rem;
}
</style>
"""
st.markdown(_HIDE_STYLES, unsafe_allow_html=True)
# ---------------------------------------------------------------------------
# Example workflows (verified cases) — prefill machinery
# ---------------------------------------------------------------------------

def _apply_pending_prefill() -> None:
    """Apply a queued example/preset BEFORE any widget is instantiated.

    Streamlit raises if a widget's session-state key is written after the
    widget exists, so pending prefill values are applied here, at the very
    top of the script.
    """
    pend = st.session_state.pop("pending_prefill", None)
    if not pend:
        return
    if pend.get("mode"):
        st.session_state["analysis_mode_sel"] = pend["mode"]
    for k, v in pend.get("widgets", {}).items():
        st.session_state[k] = v
    if "reverse_disease" in pend.get("widgets", {}):
        # A stale disambiguation pick must not survive a new disease query.
        st.session_state.pop("reverse_disease_pick", None)
        st.session_state.pop("reverse_disease_id", None)


_apply_pending_prefill()


@st.dialog("🧪 Try an example", width="small")
def _example_dialog() -> None:
    """Verified example workflows; inputs/params prefilled, never auto-run."""
    st.caption(
        "Verified example workflows. Inputs and parameters are prefilled — "
        "press **🚀 Run Analysis** in the sidebar to reproduce."
    )
    if st.button("🍵 Gui Pi Decoction · Forward (HVP5811, score ≥ 900)",
                 use_container_width=True):
        st.session_state.pending_prefill = {
            "mode": "🚀 Forward",
            "widgets": {"forward_entity_type": "Formula (方剂)",
                        "forward_text_fml": "HVP5811"},
        }
        st.rerun()
    if st.button("🫁 Lung cancer · Reverse (MONDO:0008903, score ≥ 0.1, "
                 "top 30, combined ≥ 900)", use_container_width=True):
        st.session_state.pending_prefill = {
            "mode": "🎯 Reverse",
            "widgets": {"reverse_disease": "MONDO:0008903",
                        "reverse_disease_score": 0.10,
                        "reverse_disease_topn": 30,
                        "reverse_score": 900},
        }
        st.rerun()
    if st.button("😴 Insomnia · AI consultation (症状已填入 💬 AI 页)",
                 use_container_width=True):
        st.session_state.pending_prefill = {
            "widgets": {"ai_symptoms": "失眠多梦，心悸健忘，神疲乏力，舌淡苔薄"},
        }
        st.rerun()

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def parse_lines(text: str) -> list[str] | None:
    """Split multiline text into a cleaned list of non-empty strings."""
    if not text or not text.strip():
        return None
    vals = [v.strip() for v in text.split("\n")]
    vals = [v for v in vals if v]
    return vals if vals else None


def _validate_monomer_names(chem_names: list[str]) -> dict:
    """Check which monomer names exist in the chemical database.

    Performs a case-insensitive lookup and returns the canonical database
    name for each matched input. This lets the UI report clearly which
    monomers are missing and which will be used for the analysis.
    """
    from jnuherb._data import _load_table

    df = _load_table("jnuherb_chemicals")
    name_to_canonical: dict[str, str] = {}
    not_found: list[str] = []

    # Build a lowercase lookup to the canonical (original-cased) name.
    lower_to_canonical = {}
    for canonical in df["Name"].dropna().unique():
        lower_to_canonical[str(canonical).lower()] = str(canonical)

    for name in chem_names:
        canonical = lower_to_canonical.get(name.lower())
        if canonical:
            name_to_canonical[name] = canonical
        else:
            not_found.append(name)

    return {
        "found": name_to_canonical,
        "not_found": not_found,
        "canonical_names": list(name_to_canonical.values()),
    }


# ---------------------------------------------------------------------------
# Chinese font setup for PDF/PNG exports
# ---------------------------------------------------------------------------

_CHINESE_FONT_FAMILY = "Noto Sans CJK SC, sans-serif"


def _setup_chinese_fonts() -> None:
    """Register the bundled Noto Sans CJK SC font for matplotlib and set globals."""
    try:
        from importlib.resources import files
        from matplotlib import font_manager

        font_path = files("jnuherb.data.fonts").joinpath("NotoSansCJKsc-Regular.otf")
        if font_path.exists():
            font_manager.fontManager.addfont(str(font_path))
            plt.rcParams["font.family"] = ["Noto Sans CJK SC", "sans-serif"]
        else:
            plt.rcParams["font.family"] = ["sans-serif"]
    except Exception:
        plt.rcParams["font.family"] = ["sans-serif"]

    plt.rcParams["axes.unicode_minus"] = False


def _apply_chinese_font(fig):
    """Apply Chinese font to a plotly figure layout."""
    fig.update_layout(font=dict(family=_CHINESE_FONT_FAMILY))
    return fig


_setup_chinese_fonts()


# ---------------------------------------------------------------------------
# Matplotlib plotting helpers (publication-style theme_bw bar charts)
# ---------------------------------------------------------------------------


def _theme_bw_axes(ax):
    """Apply a theme_bw-like style to a matplotlib Axes."""
    ax.set_facecolor("white")
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(0.8)
    ax.tick_params(width=0.8, length=4, color="black")
    ax.grid(False)


def _wrap_term(
    term: str,
    max_chars: int = 21,
    max_words: Optional[int] = None,
) -> str:
    """Wrap a long enrichment term into multiple lines for cleaner y-axis labels.

    Rules:
    - Each line is at most ``max_chars`` characters long.
    - If ``max_words`` is set, each line contains at most that many words.
    - Words are not broken in the middle.
    """
    if not isinstance(term, str):
        term = str(term)
    words = term.split()
    if not words:
        return term

    lines: list[str] = []
    current: list[str] = []

    def flush():
        if current:
            lines.append(" ".join(current))
            current.clear()

    for word in words:
        # Check word count constraint
        if max_words is not None and len(current) >= max_words:
            flush()
        # Check character length constraint with the next word
        candidate = " ".join(current + [word])
        if len(candidate) > max_chars and current:
            flush()
        current.append(word)

    flush()
    return "\n".join(lines) if lines else term


def _clean_go_term(term: str) -> str:
    """Remove GO identifiers like '(GO:0001234)' or '（GO:0001234）' from term names."""
    import re

    if not isinstance(term, str):
        term = str(term)
    # Match both half-width and full-width parentheses
    return re.sub(r"[\(（]GO:\d+[\)）]", "", term).strip()


def _plot_top_n_barh(
    df: pd.DataFrame,
    x: str,
    y: str,
    title: str,
    xlabel: str,
    color: str = "#8CC6ED",
    figsize: tuple[float, float] = (6.5, 8),
) -> plt.Figure:
    """Generate a publication-style horizontal bar chart for top-N items.

    Parameters match the style requested in 修改建议.md:
    - single-color bars with black edges
    - value labels on the right
    - theme_bw axis style
    """
    plot_df = df.sort_values(x, ascending=True).copy()
    # Round display values to 2 decimals for cleaner bar labels
    plot_df[x] = plot_df[x].round(2)
    vals = plot_df[x].values
    labels = plot_df[y].values

    fig, ax = plt.subplots(figsize=figsize)
    y_pos = np.arange(len(plot_df))
    bars = ax.barh(
        y_pos,
        vals,
        height=0.55,
        color=color,
        edgecolor="black",
        linewidth=0.5,
        zorder=3,
    )

    # Value labels
    for bar, val in zip(bars, vals):
        ax.text(
            val + max(vals) * 0.015,
            bar.get_y() + bar.get_height() / 2,
            f"{val:.2f}",
            va="center",
            ha="left",
            fontsize=9,
            color="black",
        )

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=10, fontweight="bold")
    ax.set_xlabel(xlabel, fontsize=11, fontweight="bold")
    ax.set_title(title, fontsize=12, fontweight="bold", pad=12)
    ax.set_xlim(0, max(vals) * 1.18)
    _theme_bw_axes(ax)
    plt.tight_layout()
    return fig


def _plot_degree_by_category(
    node_degree: pd.DataFrame,
    title: str = "PPI Network Degree (Top 5 per Category)",
    category_order: Optional[list[str]] = None,
    category_colors: Optional[dict[str, str]] = None,
    top_n: int = 5,
    figsize: tuple[float, float] = (8, 7),
) -> plt.Figure:
    """Generate a grouped horizontal bar chart of node degrees by attribute.

    TCM > Chemicals > Proteins, fixed colors per category, category labels on
    the right.  Matches the style requested in 修改建议.md.
    """
    if category_order is None:
        category_order = ["TCM", "Chemicals", "Proteins"]
    if category_colors is None:
        category_colors = {
            "TCM": "#8CC6ED",
            "Chemicals": "#C41E3A",
            "Proteins": "#2171B5",
        }

    selected = []
    for attr in category_order:
        subset = node_degree[node_degree["Attribute"] == attr].sort_values(
            "Degree", ascending=False
        )
        n = min(top_n, len(subset))
        if n > 0:
            subset = subset.head(n).copy()
            subset["Category"] = attr
            selected.append(subset)

    if not selected:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "No data available", ha="center", va="center")
        _theme_bw_axes(ax)
        return fig

    df_plot = pd.concat(selected, ignore_index=True)
    labels = df_plot["Key"].values
    degrees = df_plot["Degree"].values
    categories = df_plot["Category"].values
    bar_colors = [category_colors.get(cat, "#333333") for cat in categories]

    fig, ax = plt.subplots(figsize=figsize)
    y_pos = np.arange(len(df_plot))
    bars = ax.barh(
        y_pos,
        degrees,
        height=0.55,
        color=bar_colors,
        edgecolor="black",
        linewidth=0.5,
        zorder=3,
    )
    ax.invert_yaxis()

    max_deg = degrees.max()
    for bar, val in zip(bars, degrees):
        ax.text(
            val + max_deg * 0.015,
            bar.get_y() + bar.get_height() / 2,
            str(int(val)),
            va="center",
            ha="left",
            fontsize=10,
            color="black",
        )

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=10)

    # Category separator lines
    cat_boundaries = []
    for attr in category_order:
        idx = np.where(categories == attr)[0]
        if len(idx) > 0:
            cat_boundaries.append((idx.min(), idx.max(), attr))

    for i in range(len(cat_boundaries) - 1):
        _, bottom_max, _ = cat_boundaries[i]
        top_min, _, _ = cat_boundaries[i + 1]
        sep_y = (bottom_max + top_min) / 2 + 0.5
        ax.axhline(
            y=sep_y,
            color="grey",
            linewidth=0.8,
            linestyle="-",
            alpha=0.5,
            zorder=2,
        )

    # Right-side category labels
    for start, end, cat_name in cat_boundaries:
        mid_y = (start + end) / 2
        color = category_colors.get(cat_name, "#333333")
        ax.text(
            max_deg * 1.18,
            mid_y,
            cat_name,
            va="center",
            ha="left",
            fontsize=10,
            fontweight="bold",
            color=color,
            bbox=dict(
                boxstyle="round,pad=0.3",
                facecolor="#F5F5F5",
                edgecolor="grey",
                linewidth=0.5,
            ),
        )

    ax.set_xlabel("Degree", fontsize=12, fontweight="bold")
    ax.set_title(title, fontsize=13, fontweight="bold", pad=15)
    ax.set_xlim(0, max_deg * 1.38)
    _theme_bw_axes(ax)

    # Color y-tick labels by category
    tick_colors = [category_colors.get(cat, "#333333") for cat in categories]
    for tick, color in zip(ax.get_yticklabels(), tick_colors):
        tick.set_color(color)

    plt.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Enrichment plotting helpers (theme_bw matplotlib style)
# ---------------------------------------------------------------------------


def _plot_enrichment_barh(
    df: pd.DataFrame,
    title: str,
    color_low: str = "#8CC6ED",
    color_high: str = "#F0A0A0",
    figsize: tuple[float, float] = (6.5, 7),
) -> plt.Figure:
    """Horizontal bar chart of enrichment terms coloured by -log10(P)."""
    from matplotlib.colors import LinearSegmentedColormap

    plot_df = df.sort_values("Count", ascending=True).copy()
    terms = plot_df["Term"].values
    counts = plot_df["Count"].values
    neg_log_p = plot_df["-log10(P)"].values

    cmap_obj = LinearSegmentedColormap.from_list("enrich_cmap", [color_low, color_high])
    norm = plt.Normalize(neg_log_p.min(), neg_log_p.max())
    bar_colors = [cmap_obj(norm(v)) for v in neg_log_p]

    fig, ax = plt.subplots(figsize=figsize)
    y_pos = np.arange(len(plot_df))
    bars = ax.barh(
        y_pos,
        counts,
        height=0.55,
        color=bar_colors,
        edgecolor="black",
        linewidth=0.5,
        zorder=3,
    )

    max_count = counts.max()
    for bar, val in zip(bars, counts):
        ax.text(
            val + max_count * 0.015,
            bar.get_y() + bar.get_height() / 2,
            str(int(val)),
            va="center",
            ha="left",
            fontsize=9,
            color="black",
        )

    ax.set_yticks(y_pos)
    ax.set_yticklabels(terms, fontsize=10)
    ax.set_xlabel("Gene count", fontsize=11, fontweight="bold")
    ax.set_title(title, fontsize=12, fontweight="bold", pad=12)
    ax.set_xlim(0, max_count * 1.22)
    _theme_bw_axes(ax)

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap_obj, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, orientation="vertical", shrink=0.6, pad=0.02)
    cbar.set_label("-log10(P)", fontsize=10)

    plt.tight_layout()
    return fig


def _plot_enrichment_dot(
    df: pd.DataFrame,
    title: str,
    color_low: str = "#8CC6ED",
    color_high: str = "#F0A0A0",
    figsize: tuple[float, float] = (6.5, 7),
) -> plt.Figure:
    """Dot plot of GeneRatio vs Term, size=Count, colour=-log10(P)."""
    from matplotlib.colors import LinearSegmentedColormap

    plot_df = df.sort_values("Count", ascending=True).copy()
    terms = plot_df["Term"].values
    gene_ratio = plot_df["GeneRatio"].values
    counts = plot_df["Count"].values
    neg_log_p = plot_df["-log10(P)"].values

    cmap_obj = LinearSegmentedColormap.from_list("enrich_cmap", [color_low, color_high])
    norm = plt.Normalize(neg_log_p.min(), neg_log_p.max())
    colors = [cmap_obj(norm(v)) for v in neg_log_p]

    fig, ax = plt.subplots(figsize=figsize)
    y_pos = np.arange(len(plot_df))

    # Scale marker size by Count
    sizes = 50 + (counts - counts.min()) / max(1, counts.max() - counts.min()) * 250
    scatter = ax.scatter(
        gene_ratio,
        y_pos,
        s=sizes,
        c=colors,
        edgecolors="black",
        linewidths=0.5,
        zorder=3,
    )

    ax.set_yticks(y_pos)
    ax.set_yticklabels(terms, fontsize=10)
    ax.set_xlabel("Gene ratio", fontsize=11, fontweight="bold")
    ax.set_title(title, fontsize=12, fontweight="bold", pad=12)
    ax.set_ylim(-0.5, len(plot_df) - 0.5)
    ax.invert_yaxis()
    _theme_bw_axes(ax)

    # Legend for size
    if len(counts) > 0:
        legend_sizes = [counts.min(), int((counts.min() + counts.max()) / 2), counts.max()]
        legend_handles = [
            plt.scatter(
                [], [], s=50 + (s - counts.min()) / max(1, counts.max() - counts.min()) * 250,
                c="grey", edgecolors="black", linewidths=0.5,
            )
            for s in legend_sizes
        ]
        ax.legend(
            legend_handles,
            [str(s) for s in legend_sizes],
            title="Gene count",
            loc="lower right",
            frameon=True,
            fontsize=8,
        )

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap_obj, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, orientation="vertical", shrink=0.6, pad=0.02)
    cbar.set_label("-log10(P)", fontsize=10)

    plt.tight_layout()
    return fig


def _plot_go_barh(
    df: pd.DataFrame,
    title: str,
    category_order: list[str] = None,
    category_colors: dict[str, str] = None,
    figsize: tuple[float, float] = (7, 8),
) -> plt.Figure:
    """GO enrichment bar chart with fixed colours per ontology category."""
    if category_order is None:
        category_order = ["BP", "CC", "MF"]
    if category_colors is None:
        category_colors = {"BP": "#F0A0A0", "CC": "#8CC6ED", "MF": "#b8b4ca"}

    plot_df = df.copy()
    plot_df["Category"] = pd.Categorical(
        plot_df["Category"], categories=category_order, ordered=True
    )
    plot_df = plot_df.sort_values(["Category", "Count"], ascending=[True, True])

    terms = plot_df["Term"].values
    counts = plot_df["Count"].values
    categories = plot_df["Category"].values
    bar_colors = [category_colors.get(cat, "#333333") for cat in categories]

    fig, ax = plt.subplots(figsize=figsize)
    y_pos = np.arange(len(plot_df))
    bars = ax.barh(
        y_pos,
        counts,
        height=0.55,
        color=bar_colors,
        edgecolor="black",
        linewidth=0.5,
        zorder=3,
    )

    max_count = counts.max()
    for bar, val in zip(bars, counts):
        ax.text(
            val + max_count * 0.015,
            bar.get_y() + bar.get_height() / 2,
            str(int(val)),
            va="center",
            ha="left",
            fontsize=9,
            color="black",
        )

    ax.set_yticks(y_pos)
    ax.set_yticklabels(terms, fontsize=10)
    ax.set_xlabel("Gene count", fontsize=11, fontweight="bold")
    ax.set_title(title, fontsize=12, fontweight="bold", pad=12)
    ax.set_xlim(0, max_count * 1.25)
    _theme_bw_axes(ax)

    # Legend
    handles = [
        plt.Rectangle(
            (0, 0), 1, 1,
            facecolor=category_colors[cat],
            edgecolor="black",
            linewidth=0.5,
        )
        for cat in category_order
        if cat in category_colors
    ]
    ax.legend(
        handles,
        category_order,
        title="Ontology",
        loc="lower right",
        frameon=True,
        fontsize=9,
    )

    plt.tight_layout()
    return fig


def _plot_go_dot(
    df: pd.DataFrame,
    title: str,
    category_order: list[str] = None,
    color_low: str = "#8CC6ED",
    color_high: str = "#F0A0A0",
    figsize: tuple[float, float] = (7, 9),
) -> plt.Figure:
    """GO enrichment dot plot faceted by ontology category.

    Points are coloured by -log10(P) using a two-colour gradient.
    """
    from matplotlib.colors import LinearSegmentedColormap

    if category_order is None:
        category_order = ["BP", "CC", "MF"]

    plot_df = df.copy()
    plot_df["Category"] = pd.Categorical(
        plot_df["Category"], categories=category_order, ordered=True
    )
    plot_df = plot_df.sort_values(["Category", "Count"], ascending=[True, True])

    n_cats = len(category_order)
    fig, axes = plt.subplots(
        n_cats, 1, figsize=figsize,
        gridspec_kw={"hspace": 0.25, "top": 0.92}
    )
    if n_cats == 1:
        axes = [axes]

    counts_all = plot_df["Count"].values
    neg_log_p_all = plot_df["-log10(P)"].values
    cmap_obj = LinearSegmentedColormap.from_list("go_dot_cmap", [color_low, color_high])
    norm = plt.Normalize(neg_log_p_all.min(), neg_log_p_all.max())

    for ax, cat in zip(axes, category_order):
        sub = plot_df[plot_df["Category"] == cat]
        if sub.empty:
            ax.set_title(f"{cat} (no terms)", fontsize=11, fontweight="bold")
            _theme_bw_axes(ax)
            ax.set_yticks([])
            ax.set_xlim(0, 1)
            continue

        terms = sub["Term"].values
        gr = sub["GeneRatio"].values
        counts = sub["Count"].values
        neg_log_p = sub["-log10(P)"].values
        sizes = 50 + (counts - counts_all.min()) / max(1, counts_all.max() - counts_all.min()) * 250
        colors = [cmap_obj(norm(v)) for v in neg_log_p]

        ax.scatter(
            gr,
            np.arange(len(sub)),
            s=sizes,
            c=colors,
            edgecolors="black",
            linewidths=0.5,
            zorder=3,
        )
        ax.set_yticks(np.arange(len(sub)))
        ax.set_yticklabels(terms, fontsize=9)
        ax.set_ylim(-0.5, len(sub) - 0.5)
        ax.invert_yaxis()
        ax.set_title(cat, fontsize=11, fontweight="bold", pad=8)
        _theme_bw_axes(ax)

    axes[-1].set_xlabel("Gene ratio", fontsize=11, fontweight="bold")
    fig.suptitle(title, fontsize=12, fontweight="bold", y=1.02)

    # Shared colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap_obj, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, orientation="vertical", shrink=0.5, pad=0.02)
    cbar.set_label("-log10(P)", fontsize=10)

    return fig


def _fig_to_pdf_bytes(fig, engine: str = "plotly") -> bytes:
    """Convert a matplotlib or plotly figure to PDF bytes for download.

    Args:
        fig: A plotly.graph_objects.Figure or matplotlib.figure.Figure.
        engine: "plotly" or "matplotlib".

    Returns:
        PDF bytes. Raises RuntimeError if conversion fails.
    """
    if engine == "plotly":
        try:
            return fig.to_image(format="pdf", engine="kaleido")
        except Exception as e:
            raise RuntimeError(f"Failed to export Plotly figure to PDF: {e}")
    elif engine == "matplotlib":
        import io
        buf = io.BytesIO()
        try:
            fig.savefig(buf, format="pdf", bbox_inches="tight")
        except Exception as e:
            raise RuntimeError(f"Failed to export matplotlib figure to PDF: {e}")
        buf.seek(0)
        return buf.read()
    else:
        raise ValueError(f"Unsupported engine: {engine}")


def _fr_layout(
    edges: list[tuple[str, str]],
    nodes: list[str],
    seed: int = 42,
    iterations: int = 80,
) -> dict[str, tuple[float, float]]:
    """Small deterministic Fruchterman–Reingold spring layout in numpy.

    Dependency-free replacement for networkx.spring_layout (networkx is
    not a declared dependency).  ``nodes`` are the unique node names and
    ``edges`` are (source, target) pairs; positions are scaled to a
    ~100x100 canvas.
    """
    n = len(nodes)
    if n == 0:
        return {}
    idx = {name: i for i, name in enumerate(nodes)}
    rng = np.random.default_rng(seed)
    pos = rng.random((n, 2)) * 100.0
    edge_idx = [(idx[s], idx[t]) for s, t in edges if s in idx and t in idx]
    # Optimal pairwise distance for a pleasant spread.
    k = 100.0 / max(1.0, np.sqrt(n))
    temp = 100.0

    for _ in range(iterations):
        disp = np.zeros((n, 2))
        # Repulsive forces between all pairs.
        diff = pos[:, None, :] - pos[None, :, :]
        dist = np.linalg.norm(diff, axis=-1)
        dist = np.where(dist == 0, 1e-6, dist)
        force = k * k / dist
        with np.errstate(divide="ignore", invalid="ignore"):
            disp += (diff / dist[..., None] * force[..., None]).sum(axis=1)
        # Attractive forces along edges.
        for i, j in edge_idx:
            d = pos[i] - pos[j]
            dist = max(np.linalg.norm(d), 1e-6)
            f = dist * dist / k
            disp[i] -= d / dist * f
            disp[j] += d / dist * f
        length = np.linalg.norm(disp, axis=1, keepdims=True)
        length = np.where(length == 0, 1e-6, length)
        step = np.minimum(length, temp) * disp / length
        pos += step
        temp *= 0.95

    pos -= pos.mean(axis=0)
    return {name: (float(pos[idx[name], 0]), float(pos[idx[name], 1]))
            for name in nodes}


def _three_ring_layout(
    nodes: list[str],
    attr_map: dict[str, str],
) -> dict[str, tuple[float, float]]:
    """Layout with three SEPARATED rings: Herb on the left, Chemicals in
    the center, Proteins on the right — the three types never mix.
    Nodes with an unknown attribute join the Chemicals ring."""
    import math

    groups: dict[str, list[str]] = {
        "Herb": [], "Chemicals": [], "Proteins": [],
    }
    for n in nodes:
        g = attr_map.get(n, "Chemicals")
        g = {"TCM": "Herb"}.get(g, g)  # node Attribute uses "TCM"
        if g not in groups:
            g = "Chemicals"
        groups[g].append(n)

    radii = {g: min(280.0, max(70.0, 14.0 * len(ids) / math.pi))
             for g, ids in groups.items()}
    gap = 80.0
    order = ["Herb", "Chemicals", "Proteins"]
    width = 2 * sum(radii[g] for g in order) + 2 * gap + 120
    height = 2 * max(radii.values()) + 150
    cy = height / 2

    pos: dict[str, tuple[float, float]] = {}
    x = 60 + radii["Herb"]
    for i, g in enumerate(order):
        ids = groups[g]
        r = radii[g]
        for j, n in enumerate(ids):
            ang = -math.pi / 2 + 2 * math.pi * j / max(1, len(ids))
            pos[n] = (x + r * math.cos(ang), cy + r * math.sin(ang))
        if i < len(order) - 1:
            x += r + gap + radii[order[i + 1]]
    return pos


def _network_mpl_fig(
    nodes: list[str],
    pos: dict[str, tuple[float, float]],
    edges: list[tuple[str, str]],
    attr_map: dict[str, str],
    color_map: dict[str, str],
    title: str,
):
    """Matplotlib twin of the on-screen Plotly network, used ONLY for the
    vector-PDF download (no Chromium/kaleido dependency on the server)."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(12, 9))
    for s, t in edges:
        if s in pos and t in pos:
            x0, y0 = pos[s]
            x1, y1 = pos[t]
            ax.plot([x0, x1], [y0, y1], color="#b8c2cc", lw=0.5, zorder=1)
    xs = [pos[m][0] for m in nodes if m in pos]
    ys = [pos[m][1] for m in nodes if m in pos]
    cs = [color_map.get(attr_map.get(m, ""), "#95a5a6") for m in nodes if m in pos]
    ax.scatter(xs, ys, c=cs, s=70, edgecolors="#333333", linewidths=0.6, zorder=2)
    for m in nodes:
        if m in pos:
            ax.annotate(m, (pos[m][0], pos[m][1]),
                        fontsize=6, alpha=0.75, zorder=3)
    # Legend by attribute type.
    for attr, color in color_map.items():
        ax.scatter([], [], c=color, s=70, edgecolors="#333", label=attr)
    ax.legend(loc="upper right", fontsize=9, frameon=False)
    ax.set_title(title, fontsize=13)
    ax.axis("off")
    return fig


def _boost_mpl_fig(df: pd.DataFrame, title: str, color: str):
    """Matplotlib bar chart twin of the Plotly combination-boost charts,
    for vector-PDF download without kaleido."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 5))
    ranks = np.arange(len(df)) + 1
    ax.bar(np.arange(len(df)), df["Boost"].to_numpy(), color=color)
    ax.set_xticks(np.arange(len(df)))
    ax.set_xticklabels([str(r) for r in ranks], rotation=0)
    ax.set_xlabel("Combination Rank")
    ax.set_ylabel("Boost")
    ax.set_title(title)
    fig.tight_layout()
    return fig


def _dl_slug(download_key: Optional[str], title: str) -> str:
    """Short content-based slug for download filenames.

    ``download_key="download_kegg"`` → ``"kegg"``; falls back to a
    slugified plot title when no key is given.
    """
    src = re.sub(r"^download_", "", download_key or "") or (title or "plot")
    slug = re.sub(r"[^A-Za-z0-9]+", "_", src).strip("_").lower()
    return slug or "plot"


def _pdf_download_button(fig, filename: str, label: str, engine: str = "plotly", key: Optional[str] = None):
    """Render a download button that exports a figure as a vector PDF on click."""

    def _make_pdf():
        return _fig_to_pdf_bytes(fig, engine=engine)

    try:
        st.download_button(
            label=label,
            data=_make_pdf,
            file_name=filename,
            mime="application/pdf",
            key=key,
            use_container_width=False,
        )
    except Exception as e:
        st.caption(f"⚠️ PDF export unavailable: {e}")


def _display_enrichment(
    top_df: pd.DataFrame,
    full_df: pd.DataFrame,
    title: str,
    empty_msg: str = "No enrichment terms found.",
    download_key: Optional[str] = None,
    show_dot: bool = True,
    show_relevance: bool = True,
) -> None:
    """Render an enrichment result: table (top N), download CSV (full), bar + dot plots.

    ``show_relevance=False`` (used for KEGG) keeps the topic-relevance
    ordering but hides the Core/Related/Review tier labels and the KEGG
    BRITE class column from both the table and the CSV download.
    """
    st.markdown(f"#### {title}")
    if top_df.empty:
        st.info(empty_msg)
        return

    # Show top N table
    display_cols = [c for c in ["Term", "Category", "Relevance", "KEGG Class", "Count", "GeneRatio", "P-value", "Adjusted P-value", "Genes"] if c in top_df.columns]
    if not show_relevance:
        display_cols = [c for c in display_cols if c not in ("Relevance", "KEGG Class")]
    st.dataframe(top_df[display_cols], use_container_width=True, hide_index=True)

    # CSV download of FULL results
    if not full_df.empty:
        full_display_cols = [c for c in ["Term", "Category", "Relevance", "KEGG Class", "Count", "GeneRatio", "P-value", "Adjusted P-value", "Genes"] if c in full_df.columns]
        if not show_relevance:
            full_display_cols = [c for c in full_display_cols if c not in ("Relevance", "KEGG Class")]
        csv = full_df[full_display_cols].to_csv(index=False).encode("utf-8")
        st.download_button(
            label="📥 Download all results (CSV)",
            data=csv,
            file_name=f"{_dl_slug(download_key, title)}_all.csv",
            mime="text/csv",
            key=download_key or f"download_{title}",
            use_container_width=False,
        )

    if "-log10(P)" not in top_df.columns or "Term" not in top_df.columns:
        return

    # Prepare data: sort by Count for better bar order (top first after flip)
    plot_df = top_df.sort_values("Count", ascending=True).copy()
    # Wrap long KEGG/DO terms: at most 21 chars or 3 words per line
    plot_df["Term"] = plot_df["Term"].apply(
        lambda t: _wrap_term(t, max_chars=21, max_words=3)
    )

    # Bar plot: Term vs Count, fill by -log10(P)
    fig_bar = _plot_enrichment_barh(
        plot_df,
        title=f"{title} - Gene Count",
    )
    st.pyplot(fig_bar)
    _pdf_download_button(
        fig_bar, f"{_dl_slug(download_key, title)}_bar.pdf",
        label="📥 PDF", engine="matplotlib", key=f"pdf_bar_{download_key or title}"
    )

    if show_dot:
        # Dot plot: GeneRatio vs Term, size=Count, color=-log10(P)
        fig_dot = _plot_enrichment_dot(
            plot_df,
            title=f"{title} - Gene Ratio",
        )
        st.pyplot(fig_dot)
        _pdf_download_button(
            fig_dot, f"{_dl_slug(download_key, title)}_dot.pdf",
            label="📥 PDF", engine="matplotlib", key=f"pdf_dot_{download_key or title}"
        )


def _display_go_enrichment(
    top_df: pd.DataFrame,
    full_df: pd.DataFrame,
    title: str = "GO Enrichment",
    download_key: Optional[str] = None,
) -> None:
    """Render GO enrichment with combined table and faceted plots."""
    st.markdown(f"#### {title}")
    if top_df.empty:
        st.info("No GO enrichment terms found for the current gene set.")
        return

    # Clean GO identifiers from term names for display
    top_df = top_df.copy()
    full_df = full_df.copy()
    if "Term" in top_df.columns:
        top_df["Term"] = top_df["Term"].apply(_clean_go_term)
    if "Term" in full_df.columns:
        full_df["Term"] = full_df["Term"].apply(_clean_go_term)

    # Combined table with Category column
    display_cols = [c for c in ["Term", "Category", "Relevance", "Count", "GeneRatio", "P-value", "Adjusted P-value", "Genes"] if c in top_df.columns]
    st.dataframe(top_df[display_cols], use_container_width=True, hide_index=True)

    # CSV download of FULL results
    if not full_df.empty:
        full_display_cols = [c for c in ["Term", "Category", "Relevance", "Count", "GeneRatio", "P-value", "Adjusted P-value", "Genes"] if c in full_df.columns]
        csv = full_df[full_display_cols].to_csv(index=False).encode("utf-8")
        st.download_button(
            label="📥 Download all GO results (CSV)",
            data=csv,
            file_name="go_enrichment.csv",
            mime="text/csv",
            key=download_key or "download_go_all",
            use_container_width=False,
        )

    if "-log10(P)" not in top_df.columns or "Term" not in top_df.columns:
        return

    category_order = ["BP", "CC", "MF"]
    plot_df = top_df.sort_values(["Category", "Count"], ascending=[True, True]).copy()
    plot_df["Category"] = pd.Categorical(plot_df["Category"], categories=category_order, ordered=True)
    plot_df = plot_df.sort_values(["Category", "Count"])
    # Wrap long GO terms: at most 4 words per line
    plot_df["Term"] = plot_df["Term"].apply(
        lambda t: _wrap_term(t, max_chars=30, max_words=4)
    )

    # Category-colored bar chart
    fig_bar = _plot_go_barh(
        plot_df,
        title=f"{title} - Gene Count by Ontology",
        category_order=category_order,
    )
    st.pyplot(fig_bar)
    _pdf_download_button(
        fig_bar, "go_enrichment_bar.pdf",
        label="📥 PDF", engine="matplotlib", key=f"pdf_go_bar_{download_key or 'go'}"
    )

    # Faceted dot plot by Category
    fig_dot = _plot_go_dot(
        plot_df,
        title=f"{title} - Gene Ratio (facet by Ontology)",
        category_order=category_order,
    )
    st.pyplot(fig_dot)
    _pdf_download_button(
        fig_dot, "go_enrichment_dot.pdf",
        label="📥 PDF", engine="matplotlib", key=f"pdf_go_dot_{download_key or 'go'}"
    )


# ---------------------------------------------------------------------------
# Sidebar — Mode & Inputs (rendered on the Network Pharmacology page)
# ---------------------------------------------------------------------------


# ---- Results provenance: sidebar snapshot for the run-binding bar ----
# All sidebar widget keys whose values define "what is about to be run".
# Secrets (e.g. deepseek_api_key) intentionally live under unrelated keys
# and are never part of this snapshot.
_SIDEBAR_STATE_KEYS = (
    "analysis_mode_sel", "forward_entity_type",
    "forward_text_tcm", "forward_text_fml", "forward_text_mono",
    "forward_id_type_tcm", "forward_id_type_fml",
    "forward_proteins", "forward_score",
    "reverse_disease", "reverse_disease_score", "reverse_disease_topn",
    "reverse_disease_pick", "reverse_disease_id",
    "reverse_text", "reverse_id_type",
    "formula_text", "formula_id_type", "reverse_score",
)


def _sidebar_snapshot() -> dict:
    """Capture the current sidebar input state (raw widget values).

    Stored with each run's provenance so the results area can tell whether
    the inputs have changed since the displayed results were produced.
    """
    return {k: st.session_state.get(k) for k in _SIDEBAR_STATE_KEYS}
# ---------------------------------------------------------------------------
# Timeout helper for long-running analysis
# ---------------------------------------------------------------------------

ANALYSIS_TIMEOUT_SECONDS = 360
ANALYSIS_MAX_CONCURRENT = 2
ANALYSIS_LOCK_DIR = "/tmp/jnuherb_analysis_locks"


class _AnalysisInputError(Exception):
    """Raised for input validation errors; the UI message is shown via st.error."""
    pass


class _AnalysisDataError(Exception):
    """Raised when the requested data is not present in the database.

    The message is user-facing and is displayed via ``st.error`` in the
    main run block. Using an exception instead of calling ``st.error``
    directly ensures the message is rendered reliably even when the
    analysis runs inside a cached function.
    """
    pass


def _acquire_analysis_lock() -> int:
    """Acquire one of N file locks to allow up to N concurrent analyses.

    Uses one lock file per slot under ``ANALYSIS_LOCK_DIR``.  A slot is
    available if its file lock can be acquired exclusively without blocking.
    If all slots are busy, raises a user-facing RuntimeError.
    """
    os.makedirs(ANALYSIS_LOCK_DIR, exist_ok=True)
    for i in range(ANALYSIS_MAX_CONCURRENT):
        lock_path = os.path.join(ANALYSIS_LOCK_DIR, f"slot_{i}.lock")
        fd = os.open(lock_path, os.O_CREAT | os.O_RDWR)
        try:
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
            return fd
        except BlockingIOError:
            os.close(fd)
    raise RuntimeError(
        f"Server is busy with {ANALYSIS_MAX_CONCURRENT} analyses. "
        "Please wait a moment and try again."
    )


def _release_analysis_lock(fd: int) -> None:
    """Release a file lock slot; never raise."""
    try:
        try:
            fcntl.flock(fd, fcntl.LOCK_UN)
        except Exception:
            pass
    finally:
        try:
            os.close(fd)
        except Exception:
            pass


def _clear_stale_analysis_lock() -> None:
    """Try to remove stale lock files left by crashed sessions."""
    try:
        if os.path.exists(ANALYSIS_LOCK_DIR):
            for name in os.listdir(ANALYSIS_LOCK_DIR):
                if name.endswith(".lock"):
                    try:
                        os.remove(os.path.join(ANALYSIS_LOCK_DIR, name))
                    except Exception:
                        pass
    except Exception:
        pass


def _run_with_timeout(func, timeout: float, *args, **kwargs):
    """Run a function with a timeout using a background thread.

    This prevents a single long-running analysis from blocking the server
    indefinitely when a user closes the browser before the result returns.
    """
    with concurrent.futures.ThreadPoolExecutor(max_workers=1) as executor:
        future = executor.submit(func, *args, **kwargs)
        try:
            return future.result(timeout=timeout)
        except concurrent.futures.TimeoutError as exc:
            raise TimeoutError(
                f"Analysis exceeded {timeout:.0f} seconds. "
                "Try reducing the number of input genes/herbs/formulas or use a smaller score threshold."
            ) from exc


# ---------------------------------------------------------------------------
# Cached analysis runners
# ---------------------------------------------------------------------------

@st.cache_data(show_spinner="Running reverse analysis...")
def _run_reverse(
    proteins_tuple: tuple[str, ...],
    formula_tuple: tuple[str, ...] | None,
    protein_id_type: str,
    formula_id_type: str,
    score: float,
    num: int,
    random_state: int,
    tcm_comp: bool,
    fml_comp: bool,
    output_dir: str,
) -> dict | None:
    """Reverse: proteins → chemicals → TCMs → formulas."""
    proteins_list = list(proteins_tuple)
    formula_list = list(formula_tuple) if formula_tuple else None

    try:
        result = from_proteins(
            proteins=proteins_list,
            score=score,
            random_state=random_state,
            num=num,
            tcm_component=tcm_comp,
            formula_component=fml_comp,
            out_for_cytoscape=True,
            re=True,
            path=output_dir,
            id_type=protein_id_type,
        )
    except Exception as e:
        raise _AnalysisDataError(f"Reverse analysis failed: {e}") from e

    return _build_network_results(result, formula_list, output_dir, is_forward=False)


# English dosage-form words used in Formula_en_name / Formula_alias.  They
# are dropped before pinyin matching because users type the pinyin of the
# Chinese name ("Xiao Chai Hu Tang") while the database stores an English
# translation ("Xiao Chai Hu Decoction").
_FORM_WORDS_RE = re.compile(
    r"\b(decoctions?|pills?|granules?|capsules?|tablets?|syrups?|powders?|"
    r"injections?|liquids?|sprays?|plasters?|ointments?|drops?|wines?|teas?|"
    r"pastes?|pellets?|lozenges?|mixtures?|solutions?|suspensions?|extracts?|"
    r"tinctures?|soft|oral|effervescent|sliced|decocted)\b",
    re.IGNORECASE,
)


def _formula_pinyin_mask(fml_df: pd.DataFrame, query: str) -> pd.Series:
    """Row mask of formulas whose pinyin matches ``query``.

    Matching ignores spaces, case and English dosage-form words, in both
    directions: "Xiao Chai Hu Tang", "xiao chai hu" and "xiaochaihu" all
    find "Xiao Chai Hu Decoction"; "yougui wan" finds "You Gui Pill".
    """
    q_flat = re.sub(r"[^a-z0-9]", "", str(query).lower())
    mask = pd.Series(False, index=fml_df.index)
    if not q_flat:
        return mask
    for col in (c for c in ("Formula_en_name", "Formula_alias")
                if c in fml_df.columns):
        text = fml_df[col].fillna("").astype(str)
        # Each ;/-separated alias is a separate segment: matching and
        # reverse-containment must not leak across aliases.
        segments = text.str.split(r"[;；/]", expand=True)
        for seg_col in segments.columns:
            seg = segments[seg_col].fillna("")
            flat = (
                seg.str.replace(_FORM_WORDS_RE, " ", regex=True)
                   .str.lower()
                   .str.replace(r"[^a-z0-9]", "", regex=True)
            )
            m = flat.str.contains(re.escape(q_flat), regex=True)
            # Reverse containment (query carries the form word, DB doesn't):
            # require a minimum length so "wan"/"tang" alone can't match all.
            for idx, val in flat.items():
                if len(val) >= 4 and val in q_flat:
                    m.loc[idx] = True
            mask = mask | m
    return mask


@st.cache_data(show_spinner="Running forward analysis...")
def _run_forward(
    entities_tuple: tuple[str, ...],
    id_type: str,
    proteins_tuple: tuple[str, ...] | None,
    score: float,
    num: int,
    random_state: int,
    tcm_comp: bool,
    fml_comp: bool,
    output_dir: str,
) -> dict | None:
    """Forward: TCM/formula/monomer → chemicals → proteins.

    Routes based on id_type:
    - "Name" / "PubChem_id" → monomer query (chemical name or PubChem ID)
    - "Formula Name" / "Pinyin" → search formula by name/pinyin, then run
      forward analysis
    - everything else → standard from_tcm_or_formula
    """
    entities_list = list(entities_tuple)
    proteins_list = list(proteins_tuple) if proteins_tuple else None

    # --- Monomer: chemical name or PubChem ID query ---
    if id_type in ("Name", "PubChem_id"):
        return _run_monomer_forward(
            chem_keys=entities_list,
            chem_id_type="Name" if id_type == "Name" else "PubChem_id",
            proteins_list=proteins_list,
            score=score,
            output_dir=output_dir,
        )

    # --- Formula Name / Pinyin: resolve name → HVPID first ---
    if id_type in ("Formula Name", "Pinyin"):
        from jnuherb._data import _load_table

        fml_all = _load_table("jnuherb_formula")
        all_hvpids = []
        for name in entities_list:
            if id_type == "Formula Name":
                mask = fml_all["name"].str.contains(
                    name, na=False, case=False, regex=False
                )
            else:
                mask = _formula_pinyin_mask(fml_all, name)
            all_hvpids.extend(fml_all.loc[mask, "HVPID"].tolist())

        if not all_hvpids:
            return None  # caller handles the error message

        entities_list = list(set(all_hvpids))
        id_type = "HVMID"  # HVP prefix auto-detected as formula

    try:
        result = from_tcm_or_formula(
            tcm_or_formula_id=entities_list,
            proteins_id=proteins_list,
            score=score,
            out_for_cytoscape=True,
            out_graph=True,
            re=True,
            path=output_dir,
            id_type=id_type,
        )
    except Exception as e:
        raise _AnalysisDataError(f"Forward analysis failed: {e}") from e

    return _build_network_results(result, None, output_dir, is_forward=True)


def _run_monomer_forward(
    chem_keys: list[str],
    chem_id_type: str,
    proteins_list: list[str] | None,
    score: float,
    output_dir: str,
) -> dict | None:
    """Monomer-first forward analysis: chemical (name or PubChem ID) →
    targets → TCMs → formulas."""
    from jnuherb.query import (
        get_chemicals,
        get_tcm_chem_links,
        get_tcm,
        get_formula_tcm_links,
        get_formula,
    )
    from jnuherb.analysis import dfs_filter
    from jnuherb.score import herbiv_score
    from jnuherb.output import out_for_cyto, vis_network

    # Step 1: Find chemicals by name or PubChem ID
    chem = get_chemicals(chem_id_type, chem_keys)
    if len(chem) == 0:
        key_desc = "names" if chem_id_type == "Name" else "PubChem IDs"
        raise _AnalysisDataError(
            f"No chemicals found for the given {key_desc}: "
            f"{', '.join(chem_keys)}. "
            "Please check the spelling or try a different chemical."
        )

    # Step 2: Get TCMs containing these chemicals
    tcm_chem_links = get_tcm_chem_links("PubChem_id", list(chem["PubChem_id"]))
    if len(tcm_chem_links) == 0:
        raise _AnalysisDataError(
            "No TCMs found containing these chemicals in the database."
        )

    tcm = get_tcm("HVMID", list(tcm_chem_links["HVMID"]))

    # Step 3: Get chemical-protein links
    from jnuherb.query import get_chem_protein_links, get_proteins
    chem_prot_links = get_chem_protein_links("PubChem_id", list(chem["PubChem_id"]), score)
    if len(chem_prot_links) == 0:
        raise _AnalysisDataError(
            f"No protein targets found at score >= {score}. Try lowering the threshold."
        )

    # Step 4: Get proteins
    if proteins_list:
        proteins = get_proteins("Ensembl_ID", proteins_list)
    else:
        proteins = get_proteins("Ensembl_ID", list(chem_prot_links["Ensembl_ID"]))

    # Step 5: Get formulas containing these TCMs
    fml_tcm_links = get_formula_tcm_links("HVMID", list(tcm["HVMID"]))
    fml = get_formula("HVPID", list(fml_tcm_links["HVPID"]))

    # Step 6: DFS filter
    filtered = dfs_filter(fml, fml_tcm_links, tcm, tcm_chem_links,
                          chem, chem_prot_links, proteins)
    fml = filtered["formula"]
    fml_tcm_links = filtered["formula_tcm_links"]
    tcm = filtered["tcm"]
    tcm_chem_links = filtered["tcm_chem_links"]
    chem = filtered["chem"]
    chem_prot_links = filtered["chem_protein_links"]
    proteins = filtered["proteins"]

    # Step 7: Score
    score_result = herbiv_score(tcm, tcm_chem_links, chem, chem_prot_links,
                                formula=fml, formula_tcm_links=fml_tcm_links)
    tcm = score_result["tcm_and_score"]
    chem = score_result["chem_and_score"]
    fml = score_result["formula_and_score"]

    # Step 8: Output
    out_for_cyto(tcm, tcm_chem_links, chem, chem_prot_links, proteins, output_dir)
    vis_network(tcm, tcm_chem_links, chem, chem_prot_links, proteins, output_dir)

    result = {
        "formula": fml,
        "formula_tcm_links": fml_tcm_links,
        "tcm": tcm,
        "tcm_chem_links": tcm_chem_links,
        "chem": chem,
        "chem_protein_links": chem_prot_links,
        "proteins": proteins,
    }
    return _build_network_results(result, None, output_dir, is_forward=True)


def _build_network_results(
    api_result: dict,
    formula_filter_names: list[str] | None,
    output_dir: str,
    is_forward: bool,
) -> dict:
    """Read network files from disk and compute degree stats."""
    network_path = Path(output_dir) / "Network.csv"
    node_path = Path(output_dir) / "Type.csv"

    if network_path.exists() and node_path.exists():
        network = pd.read_csv(network_path, dtype=str, keep_default_na=False)
        node = pd.read_csv(node_path, dtype=str, keep_default_na=False)
    else:
        # Build from raw result.  Select edge endpoint columns BY NAME —
        # chem_protein_links is ordered (Ensembl_ID, Combined_score,
        # PubChem_id), so positional iloc would grab the score column and
        # drop every chemical→protein edge.
        tcl = api_result.get("tcm_chem_links")
        cpl = api_result.get("chem_protein_links")
        if tcl is not None and cpl is not None:
            _chem_col = "PubChem_id" if "PubChem_id" in cpl.columns else "HVCID"
            network = pd.DataFrame({
                "SourceNode": list(tcl["HVMID"]) + list(cpl[_chem_col]),
                "TargetNode": list(tcl[_chem_col if _chem_col in tcl.columns else tcl.columns[1]])
                           + list(cpl["Ensembl_ID"]),
            })
        else:
            network = pd.DataFrame(columns=["SourceNode", "TargetNode"])

        tcm_df = api_result.get("tcm", pd.DataFrame())
        chem_df = api_result.get("chem", pd.DataFrame())
        prot_df = api_result.get("proteins", pd.DataFrame())

        node = pd.DataFrame({
            "Key": (
                list(tcm_df.get("cn_name", tcm_df.iloc[:, 0]) if len(tcm_df) > 0 else [])
                + list(chem_df.get("Name", chem_df.iloc[:, 0]) if len(chem_df) > 0 else [])
                + list(prot_df.get("gene_name", prot_df.iloc[:, 0]) if len(prot_df) > 0 else [])
            ),
            "Attribute": (
                ["TCM"] * len(tcm_df)
                + ["Chemicals"] * len(chem_df)
                + ["Proteins"] * len(prot_df)
            ),
        })

    # SCI-friendly display: show herbs by pinyin instead of Chinese names.
    # Applied consistently to edges and node keys (before the formula
    # filter and degree computation), so the on-screen network figure, the
    # degree charts, the matplotlib PDF twins, and the CSV downloads all
    # agree.  Herbs whose pinyin collides with another herb keep their
    # Chinese name to avoid merging two nodes.
    tcm_disp = api_result.get("tcm", pd.DataFrame())
    if len(tcm_disp) > 0 and {"cn_name", "pinyin_name"} <= set(tcm_disp.columns):
        raw_map = {}
        for cn, py in zip(tcm_disp["cn_name"], tcm_disp["pinyin_name"]):
            if cn and py:
                raw_map[cn] = py.strip().title()
        dupes = {p for p in raw_map.values() if list(raw_map.values()).count(p) > 1}
        pin_map = {c: p for c, p in raw_map.items() if c != p and p not in dupes}
        if pin_map:
            network["SourceNode"] = network["SourceNode"].map(
                lambda s: pin_map.get(s, s))
            network["TargetNode"] = network["TargetNode"].map(
                lambda s: pin_map.get(s, s))
            node["Key"] = node["Key"].map(lambda s: pin_map.get(s, s))

    # Formula filter (reverse mode only)
    formula_filter = None
    if formula_filter_names and not is_forward:
        found = [f for f in formula_filter_names if f in node["Key"].values]
        if found:
            formula_filter = filter_formula_network(
                network=network, node=node, formula_names=formula_filter_names,
                include_proteins=True, path=output_dir, prefix="formula",
            )

    # Degree
    node_degree = calculate_degree(network, node)
    deg_summary = degree_summary(node_degree)

    return {
        "protein_result": api_result,
        "network": network,
        "node": node,
        "formula_filter": formula_filter,
        "node_degree": node_degree,
        "degree_summary": deg_summary,
        "output_dir": output_dir,
        "is_forward": is_forward,
    }


# ---------------------------------------------------------------------------
# Molecular docking helpers (cached structure download)
# ---------------------------------------------------------------------------

_AF_CACHE_DIR = Path(tempfile.gettempdir()) / "jnuherb_alphafold"


@st.cache_data(show_spinner="Downloading AlphaFold structure...")
def _fetch_af_pdb_cached(uniprot_id: str) -> bytes:
    from jnuherb.docking import fetch_alphafold_pdb
    return fetch_alphafold_pdb(uniprot_id, _AF_CACHE_DIR).read_bytes()


def _resolve_box(dk, rec_pdb, residues, center_manual, box_size):
    """Resolve the docking box: residues > manual center > receptor centroid."""
    size = (box_size, box_size, box_size)
    if residues is not None:
        return dk.binding_box_from_residues(rec_pdb, residues, size)
    if center_manual is not None:
        return center_manual, size
    return dk.default_box_from_receptor(rec_pdb, size)


# Receptor display styles offered in the docking result panel.
_DOCK_STYLES = ["Cartoon", "Cartoon + Surface", "Surface", "Spheres", "Sticks", "Lines"]


def _show_complex_3d(pdb_text: str, style: str = "Cartoon",
                     height: int = 440) -> None:
    """Render receptor + docked ligand with py3Dmol (graceful if missing).

    The ligand is always shown as yellow sticks; *style* controls the
    receptor representation. The 3Dmol JS is loaded from jsDelivr, which
    is reachable from networks where 3dmol.org is blocked.
    """
    try:
        import py3Dmol
    except ImportError:
        st.caption("Install `py3Dmol` to enable the interactive 3D view.")
        return
    import streamlit.components.v1 as components

    view = py3Dmol.view(
        width="100%", height=height,
        js="https://cdn.jsdelivr.net/npm/3dmol@latest/build/3Dmol-min.js",
    )
    view.addModel(pdb_text, "pdb")
    view.setStyle({"hetflag": True},
                  {"stick": {"colorscheme": "yellowCarbon", "radius": 0.25}})
    if style == "Cartoon":
        view.setStyle({"hetflag": False},
                      {"cartoon": {"color": "spectrum", "thickness": 0.2}})
    elif style == "Cartoon + Surface":
        view.setStyle({"hetflag": False},
                      {"cartoon": {"color": "spectrum", "thickness": 0.2}})
        view.addSurface(py3Dmol.VDW, {"opacity": 0.6, "color": "white"},
                        {"hetflag": False})
    elif style == "Surface":
        view.addSurface(py3Dmol.VDW, {"opacity": 0.75, "color": "spectrum"},
                        {"hetflag": False})
    elif style == "Spheres":
        view.setStyle({"hetflag": False}, {"sphere": {"colorscheme": "spectrum"}})
    elif style == "Sticks":
        view.setStyle({"hetflag": False},
                      {"stick": {"colorscheme": "spectrum", "radius": 0.15}})
    elif style == "Lines":
        view.setStyle({"hetflag": False}, {"line": {"colorscheme": "spectrum"}})
    view.zoomTo()
    components.html(view._make_html(), height=height + 20)



# ---------------------------------------------------------------------------
# Pages
# ---------------------------------------------------------------------------

def _asset_image(fname: str, max_h: int = 260):
    """Load a packaged PNG asset, downscaled for embedding."""
    from importlib.resources import files
    from PIL import Image as _PILImage
    import io as _io

    p = files("jnuherb.data.assets").joinpath(fname)
    im = _PILImage.open(_io.BytesIO(p.read_bytes())).convert("RGBA")
    if im.height > max_h:
        w = int(im.width * max_h / im.height)
        im = im.resize((w, max_h), _PILImage.LANCZOS)
    return im


def _fit_canvas(im, box=(860, 640)):
    """Centre an image on a fixed-aspect transparent canvas.

    All workflow card illustrations get the same display box, so the
    native cards keep identical width and height.
    """
    from PIL import Image as _PILImage

    canvas = _PILImage.new("RGBA", box, (0, 0, 0, 0))
    scale = min(box[0] / im.width, box[1] / im.height)
    w, h = int(im.width * scale), int(im.height * scale)
    im2 = im.resize((w, h), _PILImage.LANCZOS) if (w, h) != im.size else im
    canvas.paste(im2, ((box[0] - w) // 2, (box[1] - h) // 2), im2)
    return canvas


def _brand_icon_b64(size: int = 64) -> str:
    """Base64 of the brand icon, downscaled for inline HTML use."""
    import base64 as _b64
    import io as _io

    im = _asset_image("jnuherb-brand-icon.png", max_h=size)
    buf = _io.BytesIO()
    im.save(buf, "PNG", optimize=True)
    return _b64.b64encode(buf.getvalue()).decode()


def _entity_stats_html() -> str:
    """Six-entity stats row (icon + emphasised count + label) for the hero."""
    from importlib.resources import files
    import base64 as _b64
    import io as _io
    import json as _json

    try:
        man = _json.loads(
            files("jnuherb.data").joinpath("dataset_manifest.json").read_text("utf-8"))
        ent = man["entity_counts"]
    except Exception:  # noqa: BLE001 - never break the page for stats
        ent = {}

    hero_icons = [
        ("tcm-syndrome-meridians.png", "Syndromes",
         ent.get("syndrome_terms", {}).get("unique")),
        ("disease-medical-record.png", "Diseases",
         ent.get("disease_terms", {}).get("unique")),
        ("herbal-formula-transparent.png", "Formulas",
         ent.get("formulas", {}).get("unique")),
        ("medicinal-herb-plant.png", "Herbs",
         ent.get("herbs", {}).get("unique")),
        ("ingredient-molecule.png", "Ingredients",
         ent.get("ingredients", {}).get("unique")),
        ("protein-folded-chain.png", "Proteins",
         ent.get("proteins", {}).get("unique")),
    ]
    cells = ""
    try:
        for fname, lbl, n in hero_icons:
            im = _asset_image(fname, max_h=260)
            buf = _io.BytesIO()
            im.save(buf, "PNG", optimize=True)
            b64 = _b64.b64encode(buf.getvalue()).decode()
            num = (f'<div class="jnu-hero-ent-num">{n:,}</div>'
                   if n is not None else "")
            cells += (f'<div class="jnu-hero-ent"><img alt="{lbl}" '
                      f'src="data:image/png;base64,{b64}"/>'
                      f"{num}"
                      f'<div class="jnu-hero-ent-lbl">{lbl}</div></div>')
    except Exception:  # noqa: BLE001 - icons are decorative
        pass
    return f'<div class="jnu-hero-icons">{cells}</div>' if cells else ""


_HOME_CARDS = [
    ("💬", "AI-assisted consultation",
     "Symptoms ↔ Syndrome ↔ Candidates",
     "Find candidate herbs and formulas from symptom descriptions.",
     "ai-consultation.png", "ai", "Describe symptoms →"),
    ("🔗", "Network pharmacology",
     "TCMs ↔ Chemicals ↔ Proteins · Disease ↔ Targets",
     "Trace herbs, ingredients and target proteins — forward or reverse.",
     "network-pharmacology.png", "analysis", "Explore pathway →"),
    ("🧩", "Combination optimization",
     "Candidates ↔ Combinations",
     "Rank herb and formula combinations from analysis results.",
     "combination-optimization.png", "optimization", "Rank combinations →"),
    ("🧪", "Molecular docking",
     "Ligand + Receptor → Poses",
     "Prepare and dock one ligand against a selected receptor.",
     "molecular-docking.png", "docking", "Prepare docking →"),
]


def _home_card(emoji: str, title: str, flow: str, desc: str,
               img: str, page_key: str, link_label: str) -> None:
    """One workflow card: icon vertically centred, flow chip, native link."""
    with st.container(border=True):
        ci, ct = st.columns([1, 2.4], vertical_alignment="center")
        with ci:
            st.image(_fit_canvas(_asset_image(img, max_h=560)),
                     use_container_width=True)
        with ct:
            st.markdown(f"**{emoji} {title}**")
            st.markdown(f'<div class="jnu-flow-chip">{flow}</div>',
                        unsafe_allow_html=True)
            st.caption(desc)
        st.page_link(_PAGES[page_key], label=link_label,
                     use_container_width=False)


def _page_home() -> None:
    """Research Home: hero, entity scale, workflow cards, verified examples."""
    st.markdown(f"""
    <div class="jnu-hero">
      <div class="jnu-hero-kicker"><img class="jnu-kicker-icon" alt=""
        src="data:image/png;base64,{_brand_icon_b64()}"/> JNUHERB · Research
        Platform</div>
      <div class="jnu-hero-title">An integrated platform for traditional
        Chinese medicine research</div>
      <div class="jnu-hero-sub">AI-assisted syndrome differentiation,
        bidirectional network pharmacology, formula optimization, and
        molecular docking in one platform.</div>
      {_entity_stats_html()}
    </div>
    """, unsafe_allow_html=True)

    st.markdown('<div class="jnu-workflow-h2">Choose your workflow</div>',
                unsafe_allow_html=True)
    st.markdown('<div class="jnu-workflow-sub">Four modules, one data '
                'backbone — pick an entry point below.</div>',
                unsafe_allow_html=True)

    r1c1, r1c2 = st.columns(2)
    with r1c1:
        _home_card(*_HOME_CARDS[0])
    with r1c2:
        _home_card(*_HOME_CARDS[1])
    r2c1, r2c2 = st.columns(2)
    with r2c1:
        _home_card(*_HOME_CARDS[2])
    with r2c2:
        _home_card(*_HOME_CARDS[3])

    st.markdown('<div class="jnu-workflow-h2">🧪 Verified example '
                'workflows</div>', unsafe_allow_html=True)
    st.caption("Inputs and parameters are prefilled — press 🚀 Run Analysis "
               "on the Network Pharmacology page to reproduce.")
    e1, e2, e3 = st.columns(3)
    with e1:
        if st.button("🍵 Gui Pi Decoction", use_container_width=True):
            st.session_state.pending_prefill = {
                "mode": "🚀 Forward",
                "widgets": {"forward_entity_type": "Formula (方剂)",
                            "forward_text_fml": "HVP5811"},
            }
            st.switch_page(_PAGES["analysis"])
    with e2:
        if st.button("🫁 Lung cancer", use_container_width=True):
            st.session_state.pending_prefill = {
                "mode": "🎯 Reverse",
                "widgets": {"reverse_disease": "MONDO:0008903",
                            "reverse_disease_score": 0.10,
                            "reverse_disease_topn": 30,
                            "reverse_score": 900},
            }
            st.switch_page(_PAGES["analysis"])
    with e3:
        if st.button("😴 Insomnia", use_container_width=True):
            st.session_state.pending_prefill = {
                "widgets": {"ai_symptoms": "失眠多梦，心悸健忘，神疲乏力，舌淡苔薄"},
            }
            st.switch_page(_PAGES["ai"])

    with st.expander("📖 Detailed guide & more examples"):
        st.markdown("""
        **Reverse mode** (Disease / Proteins → TCMs):
        1. Open **Network Pharmacology** and select *🎯 Reverse*.
        2. Enter a disease name (English or 中文) or paste target protein
           IDs (Ensembl IDs or gene names).
        3. Optionally filter by specific TCM formula (expand **Filters**).
        4. Click **🚀 Run Analysis**.

        **Forward mode** (TCMs / Formulas / Monomers → Proteins):
        1. Select *🚀 Forward* and choose an entity type.
        2. **TCM**: herb IDs or names (e.g. 柴胡, HVM0367) ·
           **Formula**: HVP1625, 右归丸 · **Monomer**: deoxyshikonin.
        3. Optionally restrict to specific proteins and set a score
           threshold (expand **Filters**).
        4. Click **🚀 Run Analysis**.

        ---
        **Example — Reverse (gene names):** `TNF, IL6, TP53, PDCD1`
        **Example — Forward TCM:** `HVM0367, HVM1695` or `柴胡, 黄芩`
        **Example — Forward Formula:** `HVP1625` or `右归丸`
        **Example — Forward Monomer:** `deoxyshikonin, berberine`

        > ℹ️ The protein database contains ~14,200 proteins with known
        > TCM-related chemical interactions, not the full human proteome.
        > If your Ensembl IDs are not found, try gene names instead.
        """)

def _render_overview() -> None:
    result = st.session_state.get("result")
    is_fwd = result.get("is_forward", False) if result else False
    pr = result["protein_result"]
    ff = result["formula_filter"]

    st.subheader("📋 Analysis Summary")

    c1, c2, c3, c4 = st.columns(4)
    c1.metric("TCMs", len(pr.get("tcm", [])))
    c2.metric("Chemicals", len(pr.get("chem", [])))
    c3.metric("Proteins", len(pr.get("proteins", [])))
    if ff:
        c4.metric("Formula TCMs Found", len(ff.get("tcms_found", [])))
    elif is_fwd and pr.get("formula") is not None:
        c4.metric("Formulas Found", len(pr.get("formula", [])))
    else:
        c4.metric("Network Edges", len(result["network"]))

    st.caption("Next: explore per-entity tables in the **🌿 TCM / ⚗️ Chemical "
               "/ 🧬 Protein** tabs, the **🔗 Network**, **🧬 Enrichment** "
               "and **⚙️ Optimization**.")


def _render_consultation_report(ai: dict) -> None:
    """Render the bilingual consultation report: a full Chinese block
    first, then a full English block, each collapsible at the 🏥
    syndrome title.  Headings keep the same font size in both languages
    because both are expander labels / bold markdown."""
    syndrome_zh = ai.get("syndrome", {}).get("zh", "TCM Analysis")
    syndrome_en = ai.get("syndrome", {}).get("en", "")
    treatment = ai.get("treatment_principle", {})
    analysis_zh = ai.get("analysis", {}).get("zh", "")
    analysis_en = ai.get("analysis", {}).get("en", "")
    explanation_zh = ai.get("explanation", {}).get("zh", "")
    explanation_en = ai.get("explanation", {}).get("en", "")

    # ---- Chinese block (expanded by default) ----
    with st.expander(f"🏥 {syndrome_zh}", expanded=True):
        if treatment.get("zh"):
            st.markdown(f"**治法**\n\n{treatment['zh']}")
        if analysis_zh:
            st.markdown(f"**病机分析**\n\n{analysis_zh}")
        if explanation_zh:
            st.markdown("**用药思路**")
            st.markdown(explanation_zh)

    # ---- English block (collapsed by default) ----
    if (syndrome_en or analysis_en or explanation_en
            or treatment.get("en")):
        with st.expander(
            f"🏥 {syndrome_en or syndrome_zh}", expanded=False
        ):
            if treatment.get("en"):
                st.markdown(
                    f"**Treatment Principle**\n\n{treatment['en']}"
                )
            if analysis_en:
                st.markdown(f"**Pathogenesis Analysis**\n\n{analysis_en}")
            if explanation_en:
                st.markdown("**Treatment Rationale**")
                st.markdown(explanation_en)


def _render_ai() -> None:
    result = st.session_state.get("result")
    is_fwd = result.get("is_forward", False) if result else False
    st.subheader("🤖 AI-Powered TCM Consultation")

    # Check for saved consultation
    if "consultation" not in st.session_state:
        st.session_state.consultation = None

    if st.session_state.consultation is None:
        st.markdown("""
        Describe the patient's symptoms or disease condition below.
        The AI will analyze the presentation using **TCM syndrome
        differentiation** (辨证论治) and recommend herbs and formulas
        from the database.

        请在下方描述患者的症状或疾病情况，AI 将基于中医辨证论治进行分析，
        并从数据库中推荐相关草药与方剂。

        **Examples / 示例：**
        - *头痛发热，恶寒无汗，脉浮紧*
        - *失眠多梦，心悸健忘，神疲乏力，舌淡苔薄*
        - *胃脘胀痛，嗳气频作，食欲不振，大便不畅*
        """)

        col1, col2 = st.columns([3, 1])
        with col1:
            symptoms_input = st.text_area(
                "Symptoms / Disease Description (症状描述)",
                placeholder="请描述症状... / Describe symptoms here...",
                height=100,
                key="ai_symptoms",
            )
        with col2:
            st.caption("💡 **Tips / 提示**")
            st.caption("- Include key signs / 包含关键症状")
            st.caption("- Mention tongue/pulse if known / 舌脉信息有助于辨证")
            st.caption("- More detail = better analysis / 越详细越准确")

        if st.button("🔍 Analyze Symptoms", type="primary", key="ai_analyze_btn"):
            if not symptoms_input.strip():
                st.warning("Please enter symptom description.")
            else:
                try:
                    from jnuherb.chat import consult, check_database

                    with st.spinner("Analyzing by JNUHERB..."):
                        ai_result = consult(
                            symptoms_input.strip(),
                            api_key=st.session_state.get("deepseek_api_key"),
                        )

                    # Daily pool exhausted (server-side refusal, no API call
                    # was made): surface the fixed message and stop here.
                    if ai_result.get("quota_exceeded"):
                        st.error(f"😔 {ai_result['quota_exceeded']}")
                        st.stop()

                    st.session_state.consultation = ai_result

                    # Cross-reference with database (use Chinese names for matching)
                    with st.spinner("Cross-referencing with database..."):
                        db_check = check_database(
                            ai_result.get("herbs", {}).get("zh", []),
                            ai_result.get("formulas", {}).get("zh", []),
                        )
                    st.session_state.consultation["db_check"] = db_check
                    st.rerun()

                except ValueError as e:
                    st.error(
                        f"⚠️ {e}\n\n"
                        "The AI service is not configured on this server."
                    )
                except RuntimeError as e:
                    st.error(f"❌ {e}")
    else:
        # Show consultation results — full Chinese block first, then a
        # full English block, each collapsible at the 🏥 syndrome title.
        ai = st.session_state.consultation

        c1, c2 = st.columns([2, 1])
        with c1:
            # Brand/identity/off-topic reply: single fixed message, no
            # bilingual report blocks.
            if ai.get("identity_reply"):
                st.markdown("### 🏥 JNUHERB")
                st.info(ai["identity_reply"])
            else:
                _render_consultation_report(ai)

            # Validation warnings
            warnings = ai.get("validation_warnings", [])
            if warnings:
                st.markdown("**⚠️ Validation / 辨证校验**")
                for w in warnings:
                    st.warning(w)
            notes = ai.get("validation_notes", [])
            for n in notes:
                st.info(n)
        with c2:
            if st.button("🔄 New Consultation", key="ai_reset"):
                st.session_state.consultation = None
                st.rerun()

        st.markdown("---")

        # Database cross-reference
        db = ai.get("db_check", {})
        herbs_found = db.get("herbs_found", [])
        herbs_missing = db.get("herbs_not_found", [])
        formulas_found = db.get("formulas_found", [])
        formulas_missing = db.get("formulas_not_found", [])

        col_h, col_f = st.columns(2)

        with col_h:
            st.subheader("🌿 Recommended Herbs / 推荐草药")
            ai_herbs_en = ai.get("herbs", {}).get("en", [])
            if herbs_found:
                herb_df_data = []
                for i, h in enumerate(herbs_found):
                    herb_df_data.append({
                        "Herb (中文)": h["name"],
                        "Herb (English)": ai_herbs_en[i] if i < len(ai_herbs_en) else h.get("en_name", ""),
                        "HVMID": h["HVMID"],
                        "Pinyin": h.get("pinyin", ""),
                    })
                st.dataframe(herb_df_data, use_container_width=True)

                # Button to analyze these herbs — prefills Forward mode and
                # jumps straight to the Network Pharmacology page (staying
                # on the AI page would leave the user with no visible
                # feedback, since the sidebar inputs live on that page).
                herb_names = [h["name"] for h in herbs_found]
                if st.button("🔬 Analyze These Herbs / 分析这些草药", key="ai_analyze_herbs"):
                    st.session_state.pending_prefill = {
                        "mode": "🚀 Forward",
                        "widgets": {"forward_entity_type": "Herb (草药)",
                                    "forward_text_tcm": "\n".join(herb_names)},
                    }
                    st.switch_page(_PAGES["analysis"])
            if herbs_missing:
                st.caption(f"⚠️ Not in database / 不在数据库中: {', '.join(herbs_missing)}")

        with col_f:
            st.subheader("📜 Recommended Formulas / 推荐方剂")
            ai_formulas_en = ai.get("formulas", {}).get("en", [])
            if formulas_found:
                fml_df_data = []
                for i, f in enumerate(formulas_found):
                    source_text = f.get("source", "")
                    fml_df_data.append({
                        "Formula (中文)": f["full_name"],
                        "Formula (English)": ai_formulas_en[i] if i < len(ai_formulas_en) else "",
                        "HVPID": f["HVPID"],
                        "Source": (
                            source_text[:60] + "..."
                            if len(source_text) > 60
                            else source_text
                        ),
                    })
                st.dataframe(fml_df_data, use_container_width=True)

                # Prefill by Chinese name (not HVPID): the sidebar Formula
                # ID Type no longer offers an HVPID option.
                fml_names = [f["full_name"] for f in formulas_found]
                if st.button("🔬 Analyze These Formulas / 分析这些方剂", key="ai_analyze_formulas"):
                    st.session_state.pending_prefill = {
                        "mode": "🚀 Forward",
                        "widgets": {"forward_entity_type": "Formula (方剂)",
                                    "forward_text_fml": "\n".join(fml_names)},
                    }
                    st.switch_page(_PAGES["analysis"])
            if formulas_missing:
                st.caption(f"⚠️ Not in database / 不在数据库中: {', '.join(formulas_missing)}")

        st.markdown("---")
        with st.expander("🗂️ Database candidates used by AI / AI 使用的数据库候选"):
            candidates = ai.get("candidates", {})
            herb_candidates = candidates.get("herbs", [])
            formula_candidates = candidates.get("formulas", [])

            if not herb_candidates and not formula_candidates:
                st.info(
                    "No database candidates were retrieved for the current symptoms.\n\n"
                    "当前症状未从数据库中检索到候选草药/方剂。"
                )
            else:
                if herb_candidates:
                    st.markdown("**Candidate herbs / 候选草药**")
                    herb_cand_data = []
                    for h in herb_candidates[:15]:
                        herb_cand_data.append({
                            "Herb (中文)": h.get("cn_name", ""),
                            "Herb (English)": h.get("en_name", ""),
                            "HVMID": h.get("HVMID", ""),
                            "Match score": h.get("match_score", 0),
                        })
                    st.dataframe(herb_cand_data, use_container_width=True)

                if formula_candidates:
                    st.markdown("**Candidate formulas / 候选方剂**")
                    fml_cand_data = []
                    for f in formula_candidates[:10]:
                        fml_cand_data.append({
                            "Formula (中文)": f.get("name", ""),
                            "HVPID": f.get("HVPID", ""),
                            "Match score": f.get("match_score", 0),
                            "In DB": "否 / No" if f.get("_synthetic", False) else "是 / Yes",
                        })
                    st.dataframe(fml_cand_data, use_container_width=True)


def _render_tcm() -> None:
    result = st.session_state.get("result")
    is_fwd = result.get("is_forward", False) if result else False
    if result is None:
        st.info("Run analysis first.")
    else:
        tcm_df = result["protein_result"].get("tcm")
        if tcm_df is not None and len(tcm_df) > 0:
            tcm_df = tcm_df.copy()
            disp_cols = [c for c in ["cn_name", "pinyin_name", "en_name", "Importance Score"]
                        if c in tcm_df.columns]
            tcm_disp = tcm_df[disp_cols].sort_values("Importance Score", ascending=False)

            st.dataframe(tcm_disp, use_container_width=True, height=400)

            y_col = "cn_name" if "cn_name" in tcm_disp.columns else tcm_disp.columns[0]
            fig = _plot_top_n_barh(
                tcm_disp.head(15),
                x="Importance Score",
                y=y_col,
                title="Top 15 Herbs",
                xlabel="Importance Score",
                color="#8CC6ED",
            )
            st.pyplot(fig)
            # SCI-friendly PDF twin: pinyin y labels (screen chart keeps
            # Chinese names for readability alongside the table).
            if "pinyin_name" in tcm_disp.columns:
                pdf_df = tcm_disp.head(15).copy()
                pdf_df[y_col] = (pdf_df["pinyin_name"].astype(str)
                                 .str.strip().str.title())
                fig_pdf = _plot_top_n_barh(
                    pdf_df,
                    x="Importance Score",
                    y=y_col,
                    title="Top 15 Herbs",
                    xlabel="Importance Score",
                    color="#8CC6ED",
                )
            else:
                fig_pdf = fig
            _pdf_download_button(
                fig_pdf, "herbs.pdf",
                label="📥 PDF", engine="matplotlib", key="pdf_tcm_bar"
            )
        else:
            st.info("No TCM results available.")


def _render_chemical() -> None:
    result = st.session_state.get("result")
    is_fwd = result.get("is_forward", False) if result else False
    if result is None:
        st.info("Run analysis first.")
    else:
        chem_df = result["protein_result"].get("chem")
        if chem_df is not None and len(chem_df) > 0:
            chem_df = chem_df.copy()
            # --- OB/DL enrichment (additive; graceful when columns absent) ---
            _DL_COLS = [c for c in ("DL_TCMSP", "DL_HERB", "DL") if c in chem_df.columns]
            if "OB" not in chem_df.columns or not _DL_COLS:
                try:
                    from jnuherb._data import _load_table as _load_chem_tbl
                    _chem_meta = _load_chem_tbl("jnuherb_chemicals")
                    _need = ["PubChem_id", "OB"] + [c for c in ("OB_source", "DL_TCMSP", "DL_HERB", "DL") if c in _chem_meta.columns]
                    if "OB" in _chem_meta.columns and _chem_meta.columns.intersection(["DL_TCMSP", "DL_HERB", "DL"]).any():
                        chem_df = chem_df.merge(_chem_meta[_need].drop_duplicates("PubChem_id"), on="PubChem_id", how="left")
                except Exception:
                    pass
                _DL_COLS = [c for c in ("DL_TCMSP", "DL_HERB", "DL") if c in chem_df.columns]

            disp_cols = [c for c in ["PubChem_id", "Name", "Molecular_Formula", "OB", "OB_source",
                                     "DL_TCMSP", "DL_HERB", "DL", "Importance Score"]
                        if c in chem_df.columns]

            # --- user-adjustable OB/DL filters (default OFF, no fixed threshold) ---
            if "OB" in chem_df.columns or _DL_COLS:
                _fcol1, _fcol2 = st.columns(2)
                with _fcol1:
                    _ob_on = st.checkbox("Filter by OB", value=False, key="chem_ob_on",
                                         help="Enable to keep only chemicals with OB >= your chosen value.")
                    _ob_val = st.number_input("OB (%) >=", min_value=0.0, max_value=100.0,
                                              value=30.0, step=1.0, key="chem_ob_val",
                                              disabled=not _ob_on)
                with _fcol2:
                    _dl_col = None
                    if _DL_COLS:
                        _dl_label = st.selectbox(
                            "DL dataset",
                            _DL_COLS, key="chem_dl_col",
                            help="TCMSP and HERB drug-likeness use different models; pick which to filter by.")
                        _dl_col = _dl_label
                    _dl_on = st.checkbox("Filter by DL", value=False, key="chem_dl_on",
                                         help="Enable to keep only chemicals with DL >= your chosen value.",
                                         disabled=not _DL_COLS)
                    _dl_val = st.number_input("DL >=", min_value=0.0, max_value=1.0,
                                              value=0.18, step=0.01, key="chem_dl_val",
                                              disabled=not (_dl_on and _dl_col))

                _before = len(chem_df)
                _ob_cov = 0
                _dl_cov = 0
                if "OB" in chem_df.columns:
                    _ob_num = pd.to_numeric(chem_df["OB"], errors="coerce")
                    _ob_cov = int(_ob_num.notna().sum())
                    if _ob_on:
                        chem_df = chem_df[_ob_num >= _ob_val]
                if _dl_col:
                    _dl_num = pd.to_numeric(chem_df[_dl_col], errors="coerce")
                    _dl_cov = int(_dl_num.notna().sum())
                    if _dl_on:
                        chem_df = chem_df[_dl_num >= _dl_val]

                _cov_bits = []
                if _ob_cov:
                    _cov_bits.append(f"OB {_ob_cov}/{_before}")
                if _dl_cov:
                    _cov_bits.append(f"DL({_dl_col}) {_dl_cov}/{_before}")
                _cov_note = (" · data: " + ", ".join(_cov_bits)) if _cov_bits else ""
                st.caption(
                    f"Showing **{len(chem_df)}** / {_before} chemicals"
                    + (" after filtering" if (_ob_on or _dl_on) else " (no OB/DL filter applied)")
                    + _cov_note
                )
            # --- end OB/DL filters ---

            chem_disp = chem_df[disp_cols].sort_values("Importance Score", ascending=False)

            st.dataframe(chem_disp, use_container_width=True, height=400)

            y_col = "Name" if "Name" in chem_disp.columns else chem_disp.columns[0]
            fig = _plot_top_n_barh(
                chem_disp.head(15),
                x="Importance Score",
                y=y_col,
                title="Top 15 Chemicals",
                xlabel="Importance Score",
                color="#8CC6ED",
            )
            st.pyplot(fig)
            _pdf_download_button(
                fig, "chemicals.pdf",
                label="📥 PDF", engine="matplotlib", key="pdf_chem_bar"
            )
        else:
            st.info("No chemical results available.")


def _render_protein() -> None:
    result = st.session_state.get("result")
    is_fwd = result.get("is_forward", False) if result else False
    if result is None:
        st.info("Run analysis first.")
    else:
        prot_df = result["protein_result"].get("proteins")
        if prot_df is not None and len(prot_df) > 0:
            st.dataframe(prot_df, use_container_width=True, height=400)
        else:
            st.info("No protein results available.")

        node_deg = result.get("node_degree")
        if node_deg is not None and len(node_deg) > 0 and "Proteins" in node_deg["Attribute"].values:
            prot_deg = node_deg[node_deg["Attribute"] == "Proteins"].sort_values(
                "Degree", ascending=False
            ).head(15)
            fig = _plot_top_n_barh(
                prot_deg,
                x="Degree",
                y="Key",
                title="Top 15 Proteins by Degree",
                xlabel="Degree",
                color="#8CC6ED",
            )
            st.pyplot(fig)
            _pdf_download_button(
                fig, "proteins.pdf",
                label="📥 PDF", engine="matplotlib", key="pdf_protein_deg"
            )
        else:
            st.info("No protein degree data available.")


def _render_network() -> None:
    result = st.session_state.get("result")
    is_fwd = result.get("is_forward", False) if result else False
    if result is None:
        st.info("Run analysis first.")
    else:
        network_df = result.get("network")
        node_df = result.get("node_degree")

        if network_df is None or node_df is None or len(network_df) == 0:
            st.info("No network data available for visualization.")
        elif len(network_df) > 5000:
            st.info(
                f"⚠️ Network is too large ({len(network_df)} edges, "
                f"{len(node_df)} nodes) for in-browser visualization. "
                f"Please export to **Cytoscape** using the CSV files in "
                f"`{result.get('output_dir', 'results')}`."
            )
        else:
            st.subheader("Core Sub-Network (top nodes by degree)")

            # Limit to top nodes for clarity
            top_nodes = node_df.nlargest(100, "Degree")["Key"].tolist()
            sub_edges = network_df[
                network_df["SourceNode"].isin(top_nodes)
                & network_df["TargetNode"].isin(top_nodes)
            ]
            sub_nodes = node_df[node_df["Key"].isin(top_nodes)]

            if len(sub_edges) == 0:
                st.info("No edges in the top-100-node sub-network.")
            else:
                # Three separated rings (Herb / Chemicals / Proteins) so
                # the node types never mix.
                unique_nodes = list(
                    set(sub_edges["SourceNode"]) | set(sub_edges["TargetNode"])
                )
                attr_map = dict(zip(sub_nodes["Key"], sub_nodes["Attribute"]))
                # Display label: TCM -> Herb (the layout groups by these).
                disp_attr = {k: ("Herb" if v == "TCM" else v)
                             for k, v in attr_map.items()}
                pos = _three_ring_layout(unique_nodes, disp_attr)

                color_map = {"Herb": "#1FA9E9", "Chemicals": "#f1c40f",
                             "Proteins": "#e74c3c"}

                edge_x, edge_y = [], []
                for _, row in sub_edges.iterrows():
                    if row["SourceNode"] in pos and row["TargetNode"] in pos:
                        x0, y0 = pos[row["SourceNode"]]
                        x1, y1 = pos[row["TargetNode"]]
                        edge_x.extend([x0, x1, None])
                        edge_y.extend([y0, y1, None])

                edge_trace = go.Scatter(
                    x=edge_x, y=edge_y,
                    line=dict(width=0.5, color="#888"),
                    hoverinfo="none",
                    mode="lines",
                )

                node_x = [pos[n][0] for n in unique_nodes if n in pos]
                node_y = [pos[n][1] for n in unique_nodes if n in pos]
                node_colors = [
                    color_map.get(disp_attr.get(n, ""), "#95a5a6")
                    for n in unique_nodes if n in pos
                ]
                node_text = [
                    f"{n}<br>{disp_attr.get(n, '?')}<br>Degree: "
                    f"{sub_nodes.loc[sub_nodes['Key']==n, 'Degree'].values[0] if n in sub_nodes['Key'].values else '?'}"
                    for n in unique_nodes if n in pos
                ]

                node_trace = go.Scatter(
                    x=node_x, y=node_y,
                    mode="markers+text",
                    hoverinfo="text",
                    hovertext=node_text,
                    marker=dict(
                        size=8,
                        color=node_colors,
                        line=dict(width=1, color="#333"),
                    ),
                    text=[n if len(n) < 15 else n[:14] + "…" for n in unique_nodes if n in pos],
                    textposition="top center",
                    textfont=dict(size=8),
                )

                # Legend-only traces (one per node type).
                legend_traces = [
                    go.Scatter(
                        x=[None], y=[None], mode="markers",
                        marker=dict(size=8, color=c),
                        name=name, showlegend=True,
                    )
                    for name, c in color_map.items()
                ]

                fig = go.Figure(
                    data=[edge_trace, node_trace] + legend_traces,
                    layout=go.Layout(
                        title=f"Sub-network ({len(sub_edges)} edges, {len(set(unique_nodes))} nodes)",
                        showlegend=True,
                        legend=dict(orientation="h", yanchor="bottom", y=1.02,
                                    x=1, xanchor="right"),
                        hovermode="closest",
                        margin=dict(b=20, l=5, r=5, t=40),
                        xaxis=dict(showgrid=False, zeroline=False, visible=False),
                        yaxis=dict(showgrid=False, zeroline=False, visible=False,
                                   scaleanchor="x"),
                        height=640,
                    ),
                )
                _apply_chinese_font(fig)
                st.plotly_chart(fig, use_container_width=True)
                # Vector PDF via the matplotlib twin: kaleido's bundled
                # Chromium cannot launch inside the Docker container, so
                # the plotly-engine export fails there.
                _pdf_download_button(
                    _network_mpl_fig(
                        unique_nodes, pos,
                        [(row["SourceNode"], row["TargetNode"])
                         for _, row in sub_edges.iterrows()],
                        disp_attr, color_map,
                        f"Sub-network ({len(sub_edges)} edges, "
                        f"{len(set(unique_nodes))} nodes)",
                    ),
                    "network.pdf",
                    label="📥 PDF", engine="matplotlib", key="pdf_network"
                )

            st.info(
                "💡 For the **full** interactive network, export to "
                "**Cytoscape** using `Network.csv` and `Type.csv` in "
                f"the output directory. "
                f"Nodes: {len(node_df)}, Edges: {len(network_df)}."
            )


def _render_formula_filter() -> None:
    result = st.session_state.get("result")
    is_fwd = result.get("is_forward", False) if result else False
    if result is None:
        st.info("Run analysis first.")
    elif result.get("formula_filter") is None:
        if is_fwd:
            st.info("Formula filtering is available in **Reverse** mode. "
                     "Switch to Reverse mode and specify a formula to filter by.")
        else:
            # Auto-suggest database formulas composed of the result herbs
            # (以药测方): the reverse analysis tells us which herbs are
            # active against the disease; formulas made of those herbs
            # are candidate prescriptions.  Computed lazily on first
            # visit and cached in the result dict.
            sugg = result.get("formula_suggestions")
            if sugg is None:
                with st.spinner("Finding formulas composed of the result herbs..."):
                    try:
                        from jnuherb.analysis import suggest_formulas_from_tcms

                        sugg = suggest_formulas_from_tcms(result.get("tcm"))
                    except Exception as e:  # noqa: BLE001
                        sugg = pd.DataFrame()
                        st.caption(f"(suggestion lookup failed: {e})")
                result["formula_suggestions"] = sugg

            if sugg is not None and len(sugg) > 0:
                st.markdown(
                    "**Formulas composed of the result herbs** "
                    "(以药测方) — click one to filter the network to its herbs:"
                )
                cols = st.columns(4)
                for i, row in sugg.iterrows():
                    with cols[i % 4]:
                        if st.button(
                            f"{row['name']}\n({row['matched']}/{row['total']})",
                            key=f"ff_sugg_{row['HVPID']}",
                            use_container_width=True,
                            help="Herbs matched / formula composition size",
                        ):
                            node_keys = set(result["node"]["Key"].values)
                            names = [n for n in row["matched_names"]
                                     if n in node_keys]
                            result["formula_filter"] = filter_formula_network(
                                network=result["network"],
                                node=result["node"],
                                formula_names=names,
                                include_proteins=True,
                            )
                            result["formula_filter_label"] = (
                                f"{row['name']} ({row['matched']}/{row['total']})"
                            )
                            st.rerun()
                st.caption(
                    "Or enter specific TCM names in the sidebar "
                    "(*🎛️ Formula filter & score* expander) and re-run."
                )
            else:
                st.info("No database formula is composed mostly of the "
                        "result herbs. Enter TCM names in the sidebar "
                        "(*🎛️ Formula filter & score* expander) and re-run "
                        "to filter by specific herbs.")
    else:
        ff = result["formula_filter"]
        label = result.get("formula_filter_label")
        if label:
            c0, c1 = st.columns([3, 1])
            with c0:
                st.markdown(f"**Filtered by formula: {label}**")
            with c1:
                if st.button("✖ Clear filter", key="ff_clear"):
                    result["formula_filter"] = None
                    result.pop("formula_filter_label", None)
                    st.rerun()
        c1, c2, c3 = st.columns(3)
        c1.metric("Found", len(ff.get("tcms_found", [])))
        c2.metric("Missing", len(ff.get("tcms_missing", [])))
        c3.metric("Filtered Edges", len(ff["network"]))

        st.subheader("Filtered Network Edges")
        st.dataframe(ff["network"], use_container_width=True)
        st.subheader("Filtered Network Nodes")
        st.dataframe(ff["node"], use_container_width=True)


def _render_degree() -> None:
    result = st.session_state.get("result")
    is_fwd = result.get("is_forward", False) if result else False
    if result is None:
        st.info("Run analysis first.")
    else:
        node_deg = result.get("node_degree")

        if node_deg is not None and len(node_deg) > 0:
            fig = _plot_degree_by_category(
                node_deg,
                title="PPI Network Degree (Top 5 per Category)",
                category_order=["TCM", "Chemicals", "Proteins"],
                top_n=5,
            )
            st.pyplot(fig)
            _pdf_download_button(
                fig, "degree.pdf",
                label="📥 PDF", engine="matplotlib", key="pdf_degree_category"
            )

            st.subheader("Full Node Degree Table")
            st.dataframe(
                node_deg.sort_values("Degree", ascending=False),
                use_container_width=True, height=400,
            )
        else:
            st.info("Degree statistics not available.")


def _render_enrichment() -> None:
    result = st.session_state.get("result")
    is_fwd = result.get("is_forward", False) if result else False
    st.subheader("🧬 Enrichment")
    st.markdown("Based on the **Protein Results** from the current analysis.")

    if result is None:
        st.info("👈 Run a Forward/Reverse analysis first.")
    else:
        genes = extract_genes_from_result(result)
        if not genes:
            st.info("No protein genes available for enrichment.")
        else:
            # Disease-seeded reverse runs get KEGG/GO enrichment only —
            # a disease-derived gene list fed back into a disease collection
            # is circular.  The Disease collection stays available for
            # manually entered protein lists (and in Forward mode).
            run_info = st.session_state.get("run_info") or {}
            rev_input = run_info.get("input") or {}
            disease_seeded = (not is_fwd) and bool(rev_input.get("disease"))
            options = ["KEGG", "GO"] if disease_seeded else ["Disease", "KEGG", "GO"]
            prev_kind = st.session_state.get("enr_kind")
            if prev_kind not in options:
                st.session_state.enr_kind = options[0]
            enr_kind = st.radio(
                "Gene set collection",
                options,
                horizontal=True,
                key="enr_kind",
            )
            enr_key = f"enr_{enr_kind.lower()}"
            st.caption(f"{len(genes)} unique gene(s) extracted from Protein Results.")
            # Research topic for relevance-tiered ranking: the resolved
            # disease name in disease-seeded reverse runs.
            topic = rev_input.get("disease_name") if disease_seeded else None
            if enr_key not in st.session_state:
                with st.spinner(f"Running {enr_kind} enrichment..."):
                    st.session_state[enr_key] = run_enrichment(
                        genes, which=enr_kind.lower(), topic=topic)

            enr = st.session_state[enr_key]
            if topic:
                st.caption(
                    "Terms most relevant to the research topic "
                    f"**{topic}** are shown first; the complete result "
                    "set is preserved in the CSV download."
                )
            if enr_kind == "GO":
                st.write(
                    f"Gene mapping: {enr['go']['n_mapped']} / {enr['go']['n_input']} "
                    f"genes matched to GO background gene sets."
                )
                _display_go_enrichment(
                    enr["go"]["top"],
                    enr["go"]["full"],
                    title="GO Enrichment (BP/CC/MF Top 3 each)",
                    download_key="download_go",
                )
            else:
                block = "do" if enr_kind == "Disease" else "kegg"
                bg = "disease background gene sets" if enr_kind == "Disease" \
                    else "KEGG background gene sets"
                st.write(
                    f"Gene mapping: {enr[block]['n_mapped']} / {enr[block]['n_input']} "
                    f"genes matched to {bg}."
                )
                _display_enrichment(
                    enr[block]["top"],
                    enr[block]["full"],
                    f"{enr_kind} (Top 10)",
                    empty_msg=f"No {enr_kind} terms found for the current gene set.",
                    download_key=f"download_{block}",
                    show_relevance=(enr_kind != "KEGG"),
                )


def _render_optimization() -> None:
    result = st.session_state.get("result")
    is_fwd = result.get("is_forward", False) if result else False
    if result is None:
        st.info("Run analysis first.")
    elif is_fwd:
        st.info("Combination optimization is available in **Reverse** mode only. "
                 "Switch to Reverse mode to find optimal TCM/formula combinations.")
    else:
        pr = result["protein_result"]
        tcms_opt = pr.get("tcms")
        fmls_opt = pr.get("formulas")

        # ID -> name mappings so combinations show names, not HVMID/HVPID codes
        tcm_tbl = pr.get("tcm")
        fml_tbl = pr.get("formula")

        def _name_maps(tbl, id_col, cn_col, en_col=None):
            if tbl is None or len(tbl) == 0:
                return {}, {}
            cn = dict(zip(tbl[id_col].astype(str), tbl[cn_col].astype(str)))
            en = (dict(zip(tbl[id_col].astype(str), tbl[en_col].astype(str)))
                  if en_col and en_col in tbl.columns else {})
            return cn, en

        tcm_cn, tcm_en = _name_maps(tcm_tbl, "HVMID", "cn_name", "en_name")
        fml_cn, _ = _name_maps(fml_tbl, "HVPID", "name")

        def _expand_items(items_str, cn_map, en_map=None):
            """Split a ';'-joined ID list into display name columns."""
            ids = [i for i in str(items_str).split(";") if i]
            cn = "; ".join(cn_map.get(i, i) for i in ids)
            if en_map:
                en = "; ".join(en_map.get(i, "") for i in ids).strip("; ")
                return cn, en
            return cn, None

        def _display_combo(df, cn_map, en_map, cn_hdr, en_hdr):
            disp = df.copy()
            expanded = [_expand_items(s, cn_map, en_map) for s in disp["items"]]
            disp[cn_hdr] = [a for a, _ in expanded]
            if en_map:
                disp[en_hdr] = [b for _, b in expanded]
            return disp.drop(columns=["items"])

        has_tcms = tcms_opt is not None and len(tcms_opt) > 0
        has_fmls = fmls_opt is not None and len(fmls_opt) > 0

        if not has_tcms and not has_fmls:
            st.info("No optimization results. Enable optimization and re-run.")
        else:
            c1, c2 = st.columns(2)
            with c1:
                st.subheader("TCM Combinations")
                if has_tcms:
                    tcm_disp = _display_combo(tcms_opt, tcm_cn, tcm_en,
                                              "Chinese Name", "English Name")
                    st.dataframe(tcm_disp, use_container_width=True,
                                 hide_index=True)
                    fig = px.bar(
                        tcms_opt.head(20), y="Boost",
                        title="TCM Combination Boost",
                        color_discrete_sequence=["#00a65a"],
                    )
                    fig.update_layout(xaxis_title="Combination Rank")
                    _apply_chinese_font(fig)
                    st.plotly_chart(fig, use_container_width=True)
                    # matplotlib twin for the PDF (kaleido/Chromium fails
                    # inside the container).
                    _pdf_download_button(
                        _boost_mpl_fig(tcms_opt.head(20).reset_index(drop=True),
                                       "TCM Combination Boost", "#00a65a"),
                        "tcm_combination.pdf",
                        label="📥 PDF", engine="matplotlib", key="pdf_opt_tcm"
                    )
                else:
                    st.info("No TCM combination results.")
            with c2:
                st.subheader("Formula Combinations")
                if has_fmls:
                    fml_disp = _display_combo(fmls_opt, fml_cn, None,
                                              "Formula Name", "")
                    st.dataframe(fml_disp, use_container_width=True,
                                 hide_index=True)
                    fig = px.bar(
                        fmls_opt.head(20), y="Boost",
                        title="Formula Combination Boost",
                        color_discrete_sequence=["#3c8dbc"],
                    )
                    fig.update_layout(xaxis_title="Combination Rank")
                    _apply_chinese_font(fig)
                    st.plotly_chart(fig, use_container_width=True)
                    # matplotlib twin for the PDF (kaleido/Chromium fails
                    # inside the container).
                    _pdf_download_button(
                        _boost_mpl_fig(fmls_opt.head(20).reset_index(drop=True),
                                       "Formula Combination Boost", "#3c8dbc"),
                        "formula_combination.pdf",
                        label="📥 PDF", engine="matplotlib", key="pdf_opt_fml"
                    )
                else:
                    st.info("No formula combination results.")


def _render_docking() -> None:
    result = st.session_state.get("result")
    is_fwd = result.get("is_forward", False) if result else False
    st.subheader("🧪 Molecular Docking")
    st.markdown(
        "Single-ligand docking with **AutoDock Vina** (exhaustiveness 8, best "
        "pose). The docking box is defined by binding-site residue numbers "
        "(centroid of their heavy atoms) or an explicit center. The downloaded "
        "file is a complex PDB merging the receptor with the docked ligand pose."
    )

    dock_mode = st.radio(
        "Mode",
        ["Auto — analysis targets", "Upload — your own files"],
        horizontal=True,
        help="Auto: pick a target protein from the current analysis and the "
             "receptor structure comes from AlphaFold DB. Upload: provide "
             "your own ligand and receptor files.",
    )

    # ---------- shared box settings ----------
    box_col1, box_col2 = st.columns(2)
    with box_col1:
        residues_text = st.text_input(
            "Binding-site residues (comma- or hyphen-separated)",
            placeholder="e.g. 163,168,194 or 163-165",
            help="Residue numbers whose heavy-atom centroid becomes the box "
                 "center. Leave empty to specify the center manually below.",
        )
    with box_col2:
        box_size = st.number_input(
            "Box size (Å, cubic)", min_value=10, max_value=60, value=20, step=5,
        )

    center_manual = None
    if not residues_text.strip():
        st.caption("No residues given — the box will be centered on the "
                   "receptor's C-alpha centroid (blind docking).")
        if st.checkbox("Specify box center manually instead"):
            c1, c2, c3 = st.columns(3)
            center_manual = (
                c1.number_input("Center X", value=0.0, format="%.3f"),
                c2.number_input("Center Y", value=0.0, format="%.3f"),
                c3.number_input("Center Z", value=0.0, format="%.3f"),
            )

    # ---------- mode-specific inputs ----------
    dock_input = None  # dict gathered for the run

    if dock_mode.startswith("Auto"):
        if result is None:
            st.info("👈 Run a Forward/Reverse analysis first — the protein "
                    "list and compound scores come from its results.")
        else:
            pr = result["protein_result"]
            prot_df = pr.get("proteins")
            chem_df = pr.get("chem")
            if prot_df is None or len(prot_df) == 0:
                st.info("No proteins in the current analysis result.")
            elif chem_df is None or len(chem_df) == 0:
                st.info("No chemicals in the current analysis result.")
            else:
                prot_df = prot_df.copy()
                prot_df["_label"] = (
                    prot_df.get("gene_name", prot_df["Ensembl_ID"]).astype(str)
                    + " (" + prot_df["Ensembl_ID"].astype(str) + ")"
                )
                sel_label = st.selectbox("Target protein (from analysis)", prot_df["_label"])
                sel_prot = prot_df.loc[prot_df["_label"] == sel_label, "Ensembl_ID"].iloc[0]

                # Rank chemicals by this protein's HerbiV score
                chem_sel = chem_df.copy()
                score_col = f"{sel_prot} HerbiV Score"
                if score_col in chem_sel.columns:
                    chem_sel["_score"] = pd.to_numeric(chem_sel[score_col], errors="coerce")
                else:
                    fallback = [c for c in chem_sel.columns if c.endswith("HerbiV Score")]
                    if fallback:
                        chem_sel["_score"] = pd.to_numeric(chem_sel[fallback[0]], errors="coerce")
                    else:
                        chem_sel["_score"] = 0.0
                if "SMILES" not in chem_sel.columns:
                    st.info("Chemical SMILES not available in the current result.")
                    chem_sel = None
                else:
                    chem_sel = chem_sel[chem_sel["SMILES"].astype(str).str.len() > 3]
                    if len(chem_sel) == 0:
                        st.info("No chemicals with usable SMILES in the current result.")
                        chem_sel = None
                    else:
                        chem_sel = chem_sel.sort_values("_score", ascending=False)

                if chem_sel is not None:
                    def _chem_label(row):
                        score = row["_score"]
                        score_txt = f"{score:.3f}" if pd.notna(score) else "n/a"
                        return f"{row['Name']} [{score_txt}]"

                    chem_sel["_label"] = chem_sel.apply(_chem_label, axis=1)
                    chem_pick = st.selectbox(
                        "Ligand (ranked by HerbiV score)", chem_sel["_label"],
                        help="Defaults to the highest-scored monomer against the "
                             "selected protein.",
                    )
                    lig_row = chem_sel.loc[chem_sel["_label"] == chem_pick].iloc[0]
                    dock_input = {
                        "kind": "auto",
                        "ensembl_id": sel_prot,
                        "protein_label": sel_label,
                        "ligand_name": str(lig_row["Name"]),
                        "smiles": str(lig_row["SMILES"]),
                    }
                    st.caption(f"Ligand SMILES: `{dock_input['smiles']}`")

    else:  # Upload mode
        up_col1, up_col2 = st.columns(2)
        with up_col1:
            lig_file = st.file_uploader(
                "Ligand file (.sdf / .mol / .pdb / SMILES.txt)", type=["sdf", "mol", "pdb", "txt"],
            )
            smiles_text = st.text_input("…or paste a SMILES string", placeholder="CCO")
        with up_col2:
            prot_file = st.file_uploader("Receptor protein (.pdb)", type=["pdb"])

        ligand_name_up = st.text_input("Ligand name (optional)", value="LIG")

        if prot_file is not None:
            tmp_rec = Path(tempfile.mkdtemp(prefix="jnuherb_up_")) / "receptor.pdb"
            tmp_rec.write_bytes(prot_file.getvalue())
            st.caption(f"Receptor: {prot_file.name} ({len(prot_file.getvalue())/1024:.0f} KB)")
            dock_input = {"kind": "upload", "receptor_pdb": tmp_rec,
                          "ligand_name": ligand_name_up.strip() or "LIG",
                          "lig_file": lig_file, "smiles_text": smiles_text}

    # ---------- run ----------
    if dock_input is not None:
        run_dock = st.button("🚀 Run Docking", type="primary")
        if run_dock:
            try:
                from jnuherb import docking as dk

                with st.spinner("Preparing structures and docking "
                                "(typically 1–5 min on this server)…"):
                    # receptor
                    if dock_input["kind"] == "auto":
                        uniprot_id = dk.ensembl_to_uniprot(dock_input["ensembl_id"])
                        rec_pdb = Path(tempfile.mkdtemp(prefix="jnuherb_dock_")) / "receptor.pdb"
                        rec_pdb.write_bytes(_fetch_af_pdb_cached(uniprot_id))
                        st.caption(f"Receptor: AlphaFold model {uniprot_id}")
                        smiles = dock_input["smiles"]
                        lig_name = dock_input["ligand_name"]
                    else:
                        rec_pdb = dock_input["receptor_pdb"]
                        if dock_input["lig_file"] is not None:
                            up = dock_input["lig_file"]
                            lig_pdbqt = dk.ligand_pdbqt_from_bytes(
                                up.getvalue(), up.name, dock_input["ligand_name"])
                            # need SMILES only for the wrapper; do it manually
                            tmp_dir = rec_pdb.parent
                            receptor_pdbqt = dk.prepare_receptor(rec_pdb, tmp_dir)
                            residues = (dk.parse_residue_list(residues_text)
                                        if residues_text.strip() else None)
                            center, size = _resolve_box(
                                dk, rec_pdb, residues, center_manual, box_size)
                            affinity, pose = dk.run_docking(
                                receptor_pdbqt, lig_pdbqt, center, size,
                                exhaustiveness=8, num_poses=1, workdir=tmp_dir)
                            complex_pdb = dk.make_complex_pdb(
                                rec_pdb, pose, dock_input["ligand_name"], affinity)
                            st.session_state.dock_result = {
                                "affinity": affinity, "complex_pdb": complex_pdb,
                                "center": center, "size": size,
                                "label": dock_input["ligand_name"],
                            }
                            st.success("Docking finished!")
                            st.rerun()
                        smiles = dock_input["smiles_text"].strip()
                        lig_name = dock_input["ligand_name"]
                        if not smiles:
                            raise dk.DockingError("Provide a ligand file or a SMILES string.")

                    # common path (auto mode, or upload-by-SMILES)
                    residues = (dk.parse_residue_list(residues_text)
                                if residues_text.strip() else None)
                    center, size = _resolve_box(
                        dk, rec_pdb, residues, center_manual, box_size)
                    out = dk.dock_smiles_to_protein(
                        smiles=smiles,
                        ligand_name=lig_name,
                        receptor_pdb=rec_pdb,
                        residues=None,  # box already resolved above
                        center=center,
                        box_size=size,
                        exhaustiveness=8,
                    )
                    st.session_state.dock_result = {
                        "affinity": out["affinity"], "complex_pdb": out["complex_pdb"],
                        "center": out["center"], "size": out["box_size"],
                        "label": lig_name,
                    }
                    st.success("Docking finished!")
                    st.rerun()
            except Exception as e:  # noqa: BLE001 - show any failure in the UI
                st.error(f"Docking failed: {e}")

    # ---------- results ----------
    dock_res = st.session_state.get("dock_result")
    if dock_res is not None:
        st.markdown("---")
        st.subheader("Docking Result")
        m1, m2, m3 = st.columns(3)
        m1.metric("Vina score (kcal/mol)", f"{dock_res['affinity']:.2f}")
        m2.metric("Box center", "({:.1f}, {:.1f}, {:.1f})".format(*dock_res["center"]))
        m3.metric("Box size (Å)", "{:.0f} × {:.0f} × {:.0f}".format(*dock_res["size"]))

        st.markdown("**3D view** — receptor (style selectable), docked ligand as yellow sticks")
        dock_style = st.selectbox("Receptor style", _DOCK_STYLES, index=0,
                                  label_visibility="collapsed")
        _show_complex_3d(dock_res["complex_pdb"], style=dock_style)

        st.download_button(
            "📥 Download complex PDB",
            data=dock_res["complex_pdb"],
            file_name=f"{dock_res['label']}_docked_complex.pdb",
            mime="chemical/x-pdb",
            use_container_width=True,
        )


def _page_analysis() -> None:
    """Network Pharmacology: inputs, run, and result exploration."""
    st.sidebar.markdown("---")

    # ---- Analysis Mode ----
    st.sidebar.subheader("Analysis Mode")
    # Short pill labels with a new widget key; direction explained by the caption.
    analysis_mode = st.sidebar.radio(
        "Choose direction",
        ["🚀 Forward", "🎯 Reverse"],
        key="analysis_mode_sel",
        horizontal=True,
        label_visibility="collapsed",
    )
    is_forward = analysis_mode.startswith("🚀")
    if is_forward:
        st.sidebar.caption("🚀 Forward: Herbs / Formulas / Monomers → Proteins")
    else:
        st.sidebar.caption("🎯 Reverse: Disease / Target Proteins → Chemicals → TCMs → Formulas")

    st.sidebar.markdown("---")

    # ---- Inputs (conditional on mode) ----

    # Defined at module level, before the widgets/run branch that call it:
    # Streamlit reruns this script top-to-bottom on every interaction, so a
    # def placed below its call site is not yet bound when that code runs
    # (NameError: name '_resolve_disease' is not defined).
    @st.cache_data(show_spinner="Resolving disease targets...")
    def _resolve_disease(disease: str, min_score: float, max_targets: int,
                         disease_id: Optional[str] = None,
                         api_key: Optional[str] = None) -> dict | None:
        """Disease name/id → seed Ensembl protein IDs (CTD / Open Targets).

        ``disease_id`` pins the lookup to one vocabulary entry when the name
        matches several diseases (see the disambiguation selectbox below).
        """
        from jnuherb.analysis import resolve_disease_to_targets
        return resolve_disease_to_targets(disease, min_score=min_score,
                                          max_targets=max_targets,
                                          disease_id=disease_id, api_key=api_key)


    # Per-session output isolation (P0): every browser session gets its own
    # run directory (session_uuid/run_id), so concurrent users can never read
    # or overwrite each other's downloads, figures, or docking outputs.
    if "_session_uuid" not in st.session_state:
        st.session_state._session_uuid = uuid.uuid4().hex[:12]


    if is_forward:
        # ===== FORWARD MODE =====
        st.sidebar.subheader("Query Entity Type")
        forward_entity_type = st.sidebar.selectbox(
            "What are you searching for?",
            ["Herb (草药)", "Formula (方剂)", "Monomer (单体)"],
            key="forward_entity_type",
        )

        if forward_entity_type == "Herb (草药)":
            forward_text = st.sidebar.text_area(
                "Herb names (one per line)",
                placeholder="柴胡\n黄芩\nhuangqin\nskullcap",
                height=110,
                key="forward_text_tcm",
            )
            forward_id_type = st.sidebar.selectbox(
                "ID Type",
                ["cn_name", "pinyin_name", "en_name", "latin_name"],
                help="cn_name=Chinese name, pinyin_name=Pinyin, "
                "en_name=English name, latin_name=Latin name",
                key="forward_id_type_tcm",
            )
        elif forward_entity_type == "Formula (方剂)":
            forward_text = st.sidebar.text_area(
                "Formula names (one per line)",
                placeholder="右归丸\n小柴胡汤\nXiao Chai Hu Tang",
                height=110,
                key="forward_text_fml",
            )
            forward_id_type = st.sidebar.selectbox(
                "ID Type",
                ["Formula Name", "Pinyin"],
                help="Formula Name=Chinese/English name (e.g. 小柴胡汤), "
                     "Pinyin=space-separated pinyin (e.g. Xiao Chai Hu Tang)",
                key="forward_id_type_fml",
            )
        else:  # Monomer
            forward_text = st.sidebar.text_area(
                "Chemical / Monomer names or PubChem IDs (one per line)",
                placeholder="deoxyshikonin\nberberine\n5280345",
                height=110,
                key="forward_text_mono",
            )
            forward_id_type = st.sidebar.selectbox(
                "ID Type",
                ["Name", "PubChem_id"],
                help="Name=chemical name (e.g. berberine), "
                     "PubChem_id=PubChem compound ID (e.g. 5280345)",
                key="forward_id_type_mono",
            )

        with st.sidebar.expander("🎛️ Filters (optional)", expanded=False):
            forward_proteins_text = st.text_area(
                "Restrict to specific proteins (one per line)",
                placeholder="ENSP00000381588\nENSP00000252519",
                height=80,
                key="forward_proteins",
            )

            # Reset the score threshold when the user switches to/from Monomer.
            # The default score for monomers is 0 (inclusive) because many monomers
            # have only low-confidence targets; for TCM/formula it defaults to 900.
            # Widget state persists across reruns, so we explicitly update it when
            # the entity type changes.
            _prev_entity = st.session_state.get("_prev_forward_entity_type")
            _is_monomer = forward_entity_type == "Monomer (单体)"
            if _prev_entity != forward_entity_type or "forward_score" not in st.session_state:
                st.session_state.forward_score = 0 if _is_monomer else 900
            st.session_state._prev_forward_entity_type = forward_entity_type

            forward_score = st.number_input(
                "Combined Score Threshold",
                min_value=0, max_value=1000,
                help="Lower = more results, Higher = more confident links",
                key="forward_score",
            )
    else:
        # ===== REVERSE MODE =====
        st.sidebar.subheader("Disease or Target Proteins")

        # Disease entry: resolves a disease name/id to seed target proteins via
        # the protein–disease link table (CTD / Open Targets).  When filled it
        # takes precedence over the manual protein list below.
        reverse_disease = st.sidebar.text_input(
            "Disease name lookup",
            placeholder="e.g. lung cancer · diabetes · 高血压",
            key="reverse_disease",
            help="Type a disease name (English or 中文) or a standard disease "
                 "ID; its top associated targets become the seed proteins. "
                 "Chinese names are translated automatically.",
        )
        if reverse_disease.strip():
            _d_col1, _d_col2 = st.sidebar.columns(2)
            reverse_disease_score = _d_col1.number_input(
                "Min evidence score", min_value=0.0, max_value=1.0, value=0.10,
                step=0.05, key="reverse_disease_score",
                help="Evidence score floor for scored disease–target "
                     "associations; curated evidence rows (no numeric score) "
                     "always pass.",
            )
            reverse_disease_topn = _d_col2.number_input(
                "Max targets", min_value=5, max_value=200, value=30, step=5,
                key="reverse_disease_topn",
            )
            try:
                # P0-4 disease disambiguation: a name can match several disease
                # IDs across vocabularies (e.g. "stroke" in MONDO and MeSH).
                # List the candidates and let the user pin one exact ID before
                # any target mixing can happen.
                from jnuherb.query import get_diseases

                _cands = get_diseases(reverse_disease.strip())
                _picked_id = None
                if len(_cands) > 1:
                    _opts = [
                        f"{r.disease_name} ({r.disease_id} · {int(r.n_links)} targets)"
                        for r in _cands.itertuples()
                    ]
                    _sel = st.sidebar.selectbox(
                        "Multiple diseases match — select one:",
                        _opts, key="reverse_disease_pick",
                        help="The same name can refer to different vocabulary "
                             "entries. Targets are resolved only for the "
                             "selected ID.",
                    )
                    _picked_id = _cands.iloc[_opts.index(_sel)]["disease_id"]
                elif len(_cands) == 1:
                    _picked_id = _cands.iloc[0]["disease_id"]
                st.session_state["reverse_disease_id"] = _picked_id

                _resolved = _resolve_disease(reverse_disease.strip(),
                                             float(reverse_disease_score),
                                             int(reverse_disease_topn),
                                             disease_id=_picked_id,
                                             api_key=st.session_state.get("deepseek_api_key"))
            except Exception as e:  # noqa: BLE001
                _resolved = None
                st.sidebar.warning(f"Disease lookup unavailable: {e}")
            if _resolved and _resolved["targets"]:
                _tr = ""
                if _resolved.get("translated_from"):
                    _tr = (f"🔄 {_resolved['translated_from']} → "
                           f"{_resolved['disease_name']} · ")
                st.sidebar.caption(
                    f"{_tr}✅ {_resolved['disease_name']} "
                    f"({_resolved['disease_id']}) → {len(_resolved['targets'])} "
                    f"seed targets "
                    f"(of {_resolved['n_total']} associated; "
                    f"top: {', '.join(_resolved.get('top_genes', [])[:5])})"
                )
            elif _resolved is not None:
                st.sidebar.caption("⚠️ No targets above threshold for this disease.")
        else:
            st.sidebar.caption("…or paste target proteins directly:")

        reverse_text = st.sidebar.text_area(
            "Protein IDs (one per line)",
            placeholder="ENSP00000381588\nENSP00000252519",
            height=120,
            key="reverse_text",
        )
        reverse_id_type = st.sidebar.selectbox(
            "Protein ID Type",
            ["Ensembl_ID", "gene_name"],
            key="reverse_id_type",
        )

        with st.sidebar.expander("🎛️ Formula filter & score (optional)", expanded=False):
            formula_text = st.text_area(
                "TCM names to filter by (one per line)",
                placeholder="熟地黄\n山药\n山茱萸",
                height=80,
                key="formula_text",
            )
            formula_id_type = st.selectbox(
                "Formula ID Type",
                ["cn_name", "HVMID", "pinyin_name", "en_name", "latin_name"],
                key="formula_id_type",
            )

            reverse_score = st.number_input(
                "Combined Score Threshold",
                value=0, min_value=0, max_value=1000,
                key="reverse_score",
            )

    # Optimization defaults (previously exposed in sidebar; now hidden)
    opt_num = 100
    opt_seed = 138192
    tcm_opt_enabled = True
    fml_opt_enabled = True

    st.sidebar.markdown("---")
    run_btn = st.sidebar.button(
        "🚀 Run Analysis", type="primary", use_container_width=True,
        disabled=st.session_state.get("running", False),
        help="Click once and wait. Do not refresh or close the browser during analysis."
    )

    # Show a prominent warning while an analysis is running
    if st.session_state.get("running"):
        st.warning(
            "⏳ Analysis is running on the server. Please keep this tab open. "
            "Closing the browser will not cancel the task immediately, but it will time out after 360 seconds."
        )


    if run_btn:
        if st.session_state.get("running"):
            st.warning("⏳ An analysis is already running. Please wait for it to finish.")
        else:
            st.session_state.running = True
            st.session_state.pop("result", None)
            st.session_state.pop("enr_results", None)  # legacy key (pre-merge)
            for _k in ("enr_disease", "enr_kegg", "enr_go"):
                st.session_state.pop(_k, None)

            run_id = uuid.uuid4().hex[:12]
            output_dir = (Path(tempfile.gettempdir()) / "jnuherb_runs"
                          / st.session_state._session_uuid / f"run_{run_id}")
            output_dir.mkdir(parents=True, exist_ok=True)
            # Run provenance: inputs/params are filled in per branch below and
            # dumped to run_manifest.json in the `finally` block either way.
            run_meta: dict = {
                "run_id": run_id,
                "mode": "forward" if is_forward else "reverse",
                "started_at": datetime.now(timezone.utc).isoformat(),
                "output_dir": str(output_dir),
                "dataset_manifest": "jnuherb/data/dataset_manifest.json",
                "input": None, "parameters": None, "status": "failed",
                "sidebar_snapshot": _sidebar_snapshot(),
            }

            lock_fd = None
            try:
                lock_fd = _acquire_analysis_lock()

                if is_forward:
                    entities = parse_lines(forward_text)
                    if not entities:
                        st.error("Please enter at least one TCM/Formula/Monomer ID or name!")
                        raise _AnalysisInputError()
                    if len(entities) > 30:
                        st.error("Please enter at most 30 herbs/formulas/monomers in Forward mode.")
                        raise _AnalysisInputError()

                    # --- Formula Name: resolve before cache ---
                    entities_for_cache = entities
                    id_type_for_cache = forward_id_type
                    if id_type_for_cache == "HVPID":
                        # from_tcm_or_formula only accepts HVMID-family id_types;
                        # an "HVP" prefix is auto-detected as a formula there.
                        id_type_for_cache = "HVMID"

                    if forward_id_type in ("Formula Name", "Pinyin"):
                        from jnuherb._data import _load_table

                        fml_df = _load_table("jnuherb_formula")
                        all_hvpids = []
                        found, not_found = [], []
                        for name in entities:
                            if forward_id_type == "Formula Name":
                                mask = fml_df["name"].str.contains(
                                    name, na=False, case=False, regex=False
                                )
                            else:
                                mask = _formula_pinyin_mask(fml_df, name)
                            matches = fml_df[mask]
                            if len(matches) > 0:
                                all_hvpids.extend(matches["HVPID"].tolist())
                                found.append(name)
                            else:
                                not_found.append(name)

                        if not_found:
                            st.warning(
                                f"Formula name(s) not found: {', '.join(not_found)}. "
                                f"Try a different spelling."
                            )
                        if not all_hvpids:
                            st.error("No matching formulas found. Check the names.")
                            raise _AnalysisInputError()

                        entities_for_cache = list(set(all_hvpids))
                        if len(entities_for_cache) > 30:
                            st.error("Resolved to too many formulas (>30). Please use more specific names.")
                            raise _AnalysisInputError()
                        id_type_for_cache = "HVMID"
                        st.info(
                            f"Resolved to {len(entities_for_cache)} formula(s): "
                            f"{', '.join(entities_for_cache[:10])}"
                            f"{'...' if len(entities_for_cache) > 10 else ''}"
                        )

                    # --- Monomer Name: validate existence before analysis ---
                    if forward_id_type == "Name":
                        validation = _validate_monomer_names(entities)
                        not_found = validation["not_found"]
                        found = validation["found"]

                        if not_found:
                            if len(not_found) == len(entities):
                                st.error(
                                    "The following monomer(s) do not exist in this database: "
                                    f"{', '.join(not_found)}. "
                                    "Please check the spelling or try a different chemical name."
                                )
                                raise _AnalysisInputError()
                            else:
                                st.warning(
                                    "The following monomer(s) were not found in this database and will be skipped: "
                                    f"{', '.join(not_found)}. "
                                    f"Proceeding with: {', '.join(found.values())}."
                                )

                        if not found:
                            st.error("No valid monomers found. Please check the names and try again.")
                            raise _AnalysisInputError()

                        entities_for_cache = list(found.values())
                        id_type_for_cache = "Name"
                        st.info(
                            f"Found {len(entities_for_cache)} monomer(s) in the database: "
                            f"{', '.join(entities_for_cache[:10])}"
                            f"{'...' if len(entities_for_cache) > 10 else ''}"
                        )

                    # --- Monomer PubChem ID: validate existence before analysis ---
                    if forward_id_type == "PubChem_id":
                        from jnuherb.query import get_chemicals

                        wanted = [str(e).strip() for e in entities]
                        chem_df = get_chemicals("PubChem_id", wanted)
                        found_ids = {str(c) for c in chem_df["PubChem_id"]}
                        ok = [e for e in wanted if e in found_ids]
                        not_found = [e for e in wanted if e not in found_ids]

                        if not_found:
                            if len(not_found) == len(wanted):
                                st.error(
                                    "The following PubChem ID(s) do not exist in this "
                                    f"database: {', '.join(not_found)}. "
                                    "Please check the IDs."
                                )
                                raise _AnalysisInputError()
                            st.warning(
                                "The following PubChem ID(s) were not found and will "
                                f"be skipped: {', '.join(not_found)}. "
                                f"Proceeding with: {', '.join(ok)}."
                            )
                        if not ok:
                            st.error("No valid PubChem IDs found.")
                            raise _AnalysisInputError()

                        entities_for_cache = ok
                        id_type_for_cache = "PubChem_id"
                        st.info(
                            f"Found {len(entities_for_cache)} monomer(s) in the database: "
                            f"{', '.join(entities_for_cache[:10])}"
                            f"{'...' if len(entities_for_cache) > 10 else ''}"
                        )

                    prot_filter = parse_lines(forward_proteins_text)
                    if prot_filter and len(prot_filter) > 50:
                        st.error("Please enter at most 50 filter proteins.")
                        raise _AnalysisInputError()

                    run_meta["input"] = list(entities_for_cache)
                    run_meta["parameters"] = {
                        "id_type": id_type_for_cache,
                        "score": float(forward_score),
                        "num": int(opt_num),
                        "random_state": int(opt_seed),
                        "filter_proteins": prot_filter or None,
                    }

                    with st.spinner("Running forward analysis (timeout 360s)..."):
                        st.session_state.result = _run_with_timeout(
                            _run_forward,
                            ANALYSIS_TIMEOUT_SECONDS,
                            entities_tuple=tuple(entities_for_cache),
                            id_type=id_type_for_cache,
                            proteins_tuple=tuple(prot_filter) if prot_filter else None,
                            score=float(forward_score),
                            num=int(opt_num),
                            random_state=int(opt_seed),
                            tcm_comp=False,
                            fml_comp=False,
                            output_dir=str(output_dir),
                        )
                    run_meta["status"] = "success"
                else:
                    # Disease-seeded run takes precedence over manual proteins.
                    disease_query = st.session_state.get("reverse_disease", "").strip()
                    resolved = None
                    if disease_query:
                        try:
                            resolved = _resolve_disease(
                                disease_query,
                                float(st.session_state.get("reverse_disease_score", 0.10)),
                                int(st.session_state.get("reverse_disease_topn", 30)),
                                disease_id=st.session_state.get("reverse_disease_id"),
                                api_key=st.session_state.get("deepseek_api_key"),
                            )
                        except Exception as e:  # noqa: BLE001
                            st.error(f"Disease lookup failed: {e}")
                            raise _AnalysisInputError()
                        if not resolved or not resolved["targets"]:
                            st.error(
                                f"Disease “{disease_query}” resolved to no targets — "
                                "check the name/ID, or lower the evidence score."
                            )
                            raise _AnalysisInputError()

                    if resolved and resolved["targets"]:
                        proteins = list(resolved["targets"])
                        _tr_note = ""
                        if resolved.get("translated_from"):
                            _tr_note = (f"(translated from "
                                        f"**{resolved['translated_from']}**) ")
                        st.info(
                            f"🎯 Disease **{resolved['disease_name']}** "
                            f"({resolved['disease_id']}) {_tr_note}"
                            f"→ seeding Reverse with "
                            f"{len(proteins)} targets "
                            f"(of {resolved['n_total']} associated; score ≥ "
                            f"{resolved['min_score']})"
                        )
                    else:
                        proteins = parse_lines(reverse_text)
                        if not proteins:
                            st.error("Please enter at least one target protein ID "
                                     "or a disease name!")
                            raise _AnalysisInputError()
                        if len(proteins) > 30:
                            st.error("Please enter at most 30 target proteins in "
                                     "Reverse mode.")
                            raise _AnalysisInputError()

                    fml_list = parse_lines(formula_text)
                    if fml_list and len(fml_list) > 30:
                        st.error("Please enter at most 30 filter formulas.")
                        raise _AnalysisInputError()

                    # P0-3: disease-seeded targets are always Ensembl protein
                    # IDs from the link table — never the user's selectbox
                    # choice, which only applies to the manual-protein branch.
                    protein_id_type_for_run = (
                        "Ensembl_ID" if (resolved and resolved["targets"])
                        else reverse_id_type
                    )

                    run_meta["input"] = {
                        "disease": disease_query or None,
                        "disease_name": (resolved or {}).get("disease_name") or None,
                        "disease_id": (resolved or {}).get("disease_id") or None,
                        "proteins": None if (resolved and resolved["targets"])
                                    else list(proteins),
                    }
                    run_meta["parameters"] = {
                        "protein_id_type": protein_id_type_for_run,
                        "formula_id_type": formula_id_type,
                        "score": float(reverse_score),
                        "num": int(opt_num),
                        "random_state": int(opt_seed),
                        "tcm_comp": bool(tcm_opt_enabled),
                        "fml_comp": bool(fml_opt_enabled),
                        "filter_formulas": fml_list or None,
                    }

                    with st.spinner("Running reverse analysis (timeout 360s)..."):
                        st.session_state.result = _run_with_timeout(
                            _run_reverse,
                            ANALYSIS_TIMEOUT_SECONDS,
                            proteins_tuple=tuple(proteins),
                            formula_tuple=tuple(fml_list) if fml_list else None,
                            protein_id_type=protein_id_type_for_run,
                            formula_id_type=formula_id_type,
                            score=float(reverse_score),
                            num=int(opt_num),
                            random_state=int(opt_seed),
                            tcm_comp=tcm_opt_enabled,
                            fml_comp=fml_opt_enabled,
                            output_dir=str(output_dir),
                        )
                    run_meta["status"] = "success"
            except _AnalysisInputError:
                pass  # Error message already displayed above
            except _AnalysisDataError as e:
                st.error(f"{e}")
                run_meta["error"] = str(e)
            except TimeoutError as e:
                st.error(f"⏱️ {e}")
                run_meta["error"] = str(e)
            except RuntimeError as e:
                st.error(f"🚦 {e}")
                run_meta["error"] = str(e)
            except Exception as e:
                st.error(f"❌ Analysis failed: {e}")
                run_meta["error"] = f"{type(e).__name__}: {e}"
            finally:
                if lock_fd is not None:
                    _release_analysis_lock(lock_fd)
                # Persist run provenance (inputs, parameters, status) next to
                # the run's own outputs. Never log API keys.
                try:
                    with open(output_dir / "run_manifest.json", "w") as _mf:
                        json.dump(run_meta, _mf, ensure_ascii=False, indent=2)
                except Exception:  # noqa: BLE001
                    pass
                if run_meta.get("status") == "success":
                    # Keep the summary for the Overview tab; cleared on the
                    # next run by the pops above.
                    st.session_state.run_info = run_meta
                st.session_state.running = False
                st.rerun()


    result = st.session_state.get("result")
    is_fwd = result.get("is_forward", False) if result else False

    # ---------------------------------------------------------------------------
    # Download results
    # ---------------------------------------------------------------------------
    if result is not None:
        st.sidebar.markdown("---")
        out_dir = Path(result["output_dir"])
        # Collect CSV files for download
        csv_files = list(out_dir.glob("*.csv"))
        if csv_files:
            st.sidebar.download_button(
                label="📥 Download Network.csv",
                data=(out_dir / "Network.csv").read_bytes() if (out_dir / "Network.csv").exists() else b"",
                file_name="network.csv",
                mime="text/csv",
                use_container_width=True,
            )

    if result is None:
        st.info("👈 Set up the analysis in the sidebar and press **🚀 Run "
                "Analysis** — or load a 🧪 example from the **Research "
                "Home** page.")
        return

    tabs = st.tabs([
        "📋 Overview",
        "🌿 TCM",
        "⚗️ Chemical",
        "🧬 Protein",
        "🔗 Network",
        "🔍 Formula Filter",
        "📊 Degree",
        "🧬 Enrichment",
    ])
    with tabs[0]:
        _render_overview()
    with tabs[1]:
        _render_tcm()
    with tabs[2]:
        _render_chemical()
    with tabs[3]:
        _render_protein()
    with tabs[4]:
        _render_network()
    with tabs[5]:
        _render_formula_filter()
    with tabs[6]:
        _render_degree()
    with tabs[7]:
        _render_enrichment()



# ---------------------------------------------------------------------------
# Target Intersection — forward entities ∩ disease targets (Venn + subnetwork)
# ---------------------------------------------------------------------------


def _venn2_fig(
    a_only: int, both: int, b_only: int,
    label_a: str, label_b: str,
    color_a: str = "#00a65a", color_b: str = "#3c8dbc",
):
    """Two-circle Venn (entity targets vs disease targets) as a
    matplotlib figure — no matplotlib-venn dependency."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle

    fig, ax = plt.subplots(figsize=(9, 6.2), dpi=150)
    ax.set_aspect("equal")
    ax.axis("off")
    r, d = 1.0, 1.15
    ax.add_patch(Circle((-d / 2, 0), r, facecolor=color_a, alpha=0.45,
                        edgecolor=color_a, lw=2))
    ax.add_patch(Circle((d / 2, 0), r, facecolor=color_b, alpha=0.45,
                        edgecolor=color_b, lw=2))
    ax.text(-d / 2 - 0.42, 0, str(a_only), ha="center", va="center",
            fontsize=16, fontweight="bold", color=color_a)
    ax.text(0, 0, str(both), ha="center", va="center",
            fontsize=20, fontweight="bold", color="#333333")
    ax.text(d / 2 + 0.42, 0, str(b_only), ha="center", va="center",
            fontsize=16, fontweight="bold", color=color_b)
    ax.text(-d / 2, -r - 0.18, label_a, ha="center", va="top", fontsize=11,
            linespacing=1.4)
    ax.text(d / 2, -r - 0.18, label_b, ha="center", va="top", fontsize=11,
            linespacing=1.4)
    ax.set_xlim(-d / 2 - r - 0.3, d / 2 + r + 0.3)
    ax.set_ylim(-r - 1.25, r + 0.35)
    return fig


def _render_isect_enrichment(result: dict) -> None:
    """Enrichment for the intersection result — GO and KEGGG only (no
    Disease Ontology).  Session keys are isect_-prefixed so results can
    never contaminate or be contaminated by the main analysis page."""
    st.subheader("🧬 Enrichment")
    st.markdown("Based on the **intersected protein targets**.")
    genes = extract_genes_from_result(result)
    if not genes:
        st.info("No protein genes available for enrichment.")
        return
    enr_kind = st.radio(
        "Gene set collection", ["KEGG", "GO"], horizontal=True,
        key="isect_enr_kind",
    )
    enr_key = f"isect_enr_{enr_kind.lower()}"
    st.caption(f"{len(genes)} unique gene(s) from the intersected targets.")
    # Research topic for relevance-tiered ranking: the Overlap disease.
    _isect = st.session_state.get("isect") or {}
    topic = ((_isect.get("resolved") or {}).get("disease_name"))
    if enr_key not in st.session_state:
        with st.spinner(f"Running {enr_kind} enrichment..."):
            st.session_state[enr_key] = run_enrichment(
                genes, which=enr_kind.lower(), topic=topic)
    enr = st.session_state[enr_key]
    if enr_kind == "GO":
        st.write(
            f"Gene mapping: {enr['go']['n_mapped']} / {enr['go']['n_input']} "
            f"genes matched to GO background gene sets."
        )
        _display_go_enrichment(
            enr["go"]["top"], enr["go"]["full"],
            title="GO Enrichment (BP/CC/MF Top 3 each)",
            download_key="isect_download_go",
        )
    else:
        st.write(
            f"Gene mapping: {enr['kegg']['n_mapped']} / {enr['kegg']['n_input']} "
            f"genes matched to KEGG background gene sets."
        )
        _display_enrichment(
            enr["kegg"]["top"], enr["kegg"]["full"],
            "KEGG (Top 10)",
            empty_msg="No KEGG terms found for the current gene set.",
            download_key="isect_download_kegg",
            show_relevance=False,
        )


def _render_intersection() -> None:
    """Forward entities ∩ disease targets: Venn diagram + filtered
    network with the standard result tabs (no DO enrichment)."""
    st.subheader("🎯 Target Intersection / 靶点交集")
    st.markdown(
        "Enter herbs / formulas / monomers **and** a disease.  The forward "
        "targets are intersected with the disease-associated targets; the "
        "overlapping sub-network is explored with the standard tabs below."
    )

    # ---- sidebar inputs ----
    st.sidebar.subheader("🎯 Overlap")
    ent_type = st.sidebar.selectbox(
        "Entity type",
        ["Herb (草药)", "Formula (方剂)", "Monomer (单体)"],
        key="isect_ent_type",
    )
    if ent_type == "Herb (草药)":
        ent_text = st.sidebar.text_area(
            "Herb names (one per line)", placeholder="柴胡\n黄芩",
            height=100, key="isect_text",
        )
        id_type = st.sidebar.selectbox(
            "ID Type", ["cn_name", "pinyin_name", "en_name", "latin_name"],
            key="isect_id_type",
        )
    elif ent_type == "Formula (方剂)":
        ent_text = st.sidebar.text_area(
            "Formula names (one per line)",
            placeholder="右归丸\nXiao Chai Hu Tang",
            height=100, key="isect_text",
        )
        id_type = st.sidebar.selectbox(
            "ID Type", ["Formula Name", "Pinyin"], key="isect_id_type",
        )
    else:
        ent_text = st.sidebar.text_area(
            "Monomer names / PubChem IDs (one per line)",
            placeholder="berberine\n2353",
            height=100, key="isect_text",
        )
        id_type = st.sidebar.selectbox(
            "ID Type", ["Name", "PubChem_id"], key="isect_id_type",
        )

    disease = st.sidebar.text_input(
        "Disease / 疾病", placeholder="e.g. osteoporosis · 高血压",
        key="isect_disease",
    )
    # Monomers usually need a low score threshold, mirroring Forward mode.
    _is_mono = ent_type == "Monomer (单体)"
    _prev = st.session_state.get("_isect_prev_type")
    if _prev != ent_type or "isect_score" not in st.session_state:
        st.session_state.isect_score = 0 if _is_mono else 900
    st.session_state._isect_prev_type = ent_type
    score = st.sidebar.number_input(
        "Combined Score Threshold", min_value=0, max_value=1000,
        key="isect_score",
    )
    run = st.sidebar.button(
        "🚀 Run Intersection", type="primary", use_container_width=True,
        key="isect_run",
    )

    # ---- run ----
    if run:
        entities = parse_lines(ent_text)
        disease_q = (disease or "").strip()
        if not entities:
            st.sidebar.error("Please enter at least one herb/formula/monomer.")
        elif not disease_q:
            st.sidebar.error("Please enter a disease name.")
        else:
            from jnuherb.analysis import resolve_disease_to_targets

            try:
                resolved = resolve_disease_to_targets(disease_q)
            except Exception as e:  # noqa: BLE001
                resolved = None
                st.error(f"Disease lookup failed: {e}")
            if resolved is not None and not resolved.get("targets"):
                st.error(
                    f"Could not resolve disease **{disease_q}** to any "
                    "associated targets. Check the spelling or try English."
                )
                resolved = None
            if resolved is not None:
                out_dir = Path(tempfile.mkdtemp(prefix="jnuherb_isect_"))
                try:
                    with st.spinner("Running forward analysis..."):
                        fwd = _run_forward(
                            entities_tuple=tuple(entities),
                            id_type=id_type,
                            proteins_tuple=None,
                            score=float(score),
                            num=0,
                            random_state=0,
                            tcm_comp=False,
                            fml_comp=False,
                            output_dir=str(out_dir),
                        )
                except Exception as e:  # noqa: BLE001
                    fwd = None
                    st.error(f"Forward analysis failed: {e}")
                if fwd is not None:
                    with st.spinner("Intersecting with disease targets..."):
                        _isect_run_intersection(fwd, resolved, entities,
                                                ent_type, str(out_dir))
                    st.rerun()

    # ---- results ----
    isect = st.session_state.get("isect")
    if not isect:
        st.info("👈 Set up the inputs in the sidebar and press "
                "**🚀 Run Intersection**.")
        return

    resolved = isect["resolved"]
    fwd_t: set = isect["fwd_targets"]
    dis_t: set = isect["disease_targets"]
    inter: set = isect["overlap"]
    _tr = ""
    if resolved.get("translated_from"):
        _tr = f" (translated from **{resolved['translated_from']}**)"
    st.markdown(
        f"**{resolved['disease_name']}** ({resolved['disease_id']}){_tr} "
        f"→ {len(dis_t)} disease targets · **{len(inter)} overlapping**"
    )

    # Venn
    vfig = _venn2_fig(
        len(fwd_t - dis_t), len(inter), len(dis_t - fwd_t),
        label_a=f"Herb\n(n={len(fwd_t)})",
        label_b=f"Disease\n(n={len(dis_t)})",
    )
    st.pyplot(vfig)
    _pdf_download_button(
        vfig, "intersection_venn.pdf", label="📥 PDF",
        engine="matplotlib", key="pdf_isect_venn",
    )

    # CSV download: one row per intersected protein with its source
    # monomers and the disease.
    st.download_button(
        "📥 Download intersected targets (CSV)",
        data=isect["overlap_df"].to_csv(index=False).encode("utf-8"),
        file_name="intersection_targets.csv",
        mime="text/csv",
        key="dl_isect_csv",
    )
    with st.expander("📋 Intersected targets table"):
        st.dataframe(isect["overlap_df"], use_container_width=True)

    if len(inter) == 0:
        st.warning(
            "No overlapping targets — try a lower Combined Score "
            "Threshold or different inputs."
        )
        return

    st.markdown("---")
    st.markdown(f"**Sub-network through the {len(inter)} intersected targets**")

    # Standard result tabs, rendered against the intersection result by
    # temporarily swapping st.session_state.result.
    res = isect["result"]
    saved = st.session_state.get("result")
    st.session_state.result = res
    try:
        tabs = st.tabs([
            "🌿 TCM", "⚗️ Chemical", "🧬 Protein",
            "🔗 Network", "📊 Degree", "🧬 Enrichment",
        ])
        with tabs[0]:
            _render_tcm()
        with tabs[1]:
            _render_chemical()
        with tabs[2]:
            _render_protein()
        with tabs[3]:
            _render_network()
        with tabs[4]:
            _render_degree()
        with tabs[5]:
            _render_isect_enrichment(res)
    finally:
        if saved is None:
            st.session_state.pop("result", None)
        else:
            st.session_state.result = saved


def _isect_run_intersection(
    fwd: dict, resolved: dict, entities: list[str],
    ent_type: str, out_dir: str,
) -> None:
    """Restrict the forward result to disease targets, rebuild the
    filtered sub-network and store everything in session state."""
    from jnuherb.analysis import dfs_filter
    from jnuherb.output import out_for_cyto, vis_network
    from jnuherb.score import herbiv_score

    pr = fwd["protein_result"]
    dis_targets = set(resolved["targets"])
    cpl_full = pr["chem_protein_links"]
    fwd_targets = set(cpl_full["Ensembl_ID"].astype(str))
    dis_targets = {str(t) for t in dis_targets}
    inter = fwd_targets & dis_targets

    cpl_i = cpl_full[cpl_full["Ensembl_ID"].astype(str).isin(dis_targets)].copy()

    if len(cpl_i) > 0:
        filt = dfs_filter(
            pr.get("formula"), pr.get("formula_tcm_links"),
            pr["tcm"], pr["tcm_chem_links"], pr["chem"], cpl_i,
            pr["proteins"],
        )
        sr = herbiv_score(
            filt["tcm"], filt["tcm_chem_links"], filt["chem"],
            filt["chem_protein_links"],
            formula=filt["formula"],
            formula_tcm_links=filt["formula_tcm_links"],
        )
        trimmed = {
            "formula": sr["formula_and_score"],
            "formula_tcm_links": filt["formula_tcm_links"],
            "tcm": sr["tcm_and_score"],
            "tcm_chem_links": filt["tcm_chem_links"],
            "chem": sr["chem_and_score"],
            "chem_protein_links": filt["chem_protein_links"],
            "proteins": filt["proteins"],
        }
        if len(trimmed["tcm_chem_links"]) > 0:
            out_for_cyto(trimmed["tcm"], trimmed["tcm_chem_links"],
                         trimmed["chem"], trimmed["chem_protein_links"],
                         trimmed["proteins"], out_dir)
            vis_network(trimmed["tcm"], trimmed["tcm_chem_links"],
                        trimmed["chem"], trimmed["chem_protein_links"],
                        trimmed["proteins"], out_dir)
    else:
        trimmed = {
            k: pr.get(k) for k in (
                "formula", "formula_tcm_links", "tcm", "tcm_chem_links",
                "chem", "chem_protein_links", "proteins")
        }

    new_res = _build_network_results(trimmed, None, out_dir, is_forward=True)

    # Overlap table: one row per intersected protein — source monomers
    # (chemicals) and the disease.
    chem = pr["chem"]
    name_of_chem = dict(zip(chem["PubChem_id"].astype(str),
                            chem["Name"].astype(str)))
    proteins = pr["proteins"]
    gene_of = dict(zip(proteins["Ensembl_ID"].astype(str),
                       proteins.get("gene_name",
                                    pd.Series([""] * len(proteins))).astype(str)))
    cpl_full = cpl_full.assign(
        Ensembl_ID=cpl_full["Ensembl_ID"].astype(str),
        PubChem_id=cpl_full["PubChem_id"].astype(str),
    )
    rows = []
    for ens in sorted(inter):
        mono_ids = cpl_full.loc[cpl_full["Ensembl_ID"] == ens, "PubChem_id"]
        monomers = sorted({name_of_chem.get(m, m) for m in mono_ids})
        rows.append({
            "Ensembl_ID": ens,
            "gene_name": gene_of.get(ens, ""),
            "monomers": "; ".join(monomers),
            "disease_name": resolved["disease_name"],
            "disease_id": resolved["disease_id"],
        })
    overlap_df = pd.DataFrame(rows)

    st.session_state.isect = {
        "result": new_res,
        "resolved": resolved,
        "fwd_targets": fwd_targets,
        "disease_targets": dis_targets,
        "overlap": inter,
        "entities": entities,
        "entity_type": ent_type,
        "overlap_df": overlap_df,
    }


# ---------------------------------------------------------------------------
# Native navigation
# ---------------------------------------------------------------------------
_PAGES = {
    "home": st.Page(_page_home, title="Home", icon="🏠"),
    "ai": st.Page(_render_ai, title="Assistant", icon="💬"),
    "analysis": st.Page(_page_analysis, title="NetPharm", icon="🔗"),
    "optimization": st.Page(_render_optimization, title="Optimization", icon="⚙️"),
    "intersection": st.Page(_render_intersection, title="Overlap", icon="🎯"),
    "docking": st.Page(_render_docking, title="Docking", icon="🧪"),
}

# Brand (icon + wordmark) pinned above the native navigation, on every page.
try:
    st.logo(_res_files("jnuherb.data.assets").joinpath(
        "jnuherb-brand.png").read_bytes())
except Exception:  # noqa: BLE001 - fall back to a plain text brand
    st.sidebar.markdown("## 🌿 JNUHERB")

pg = st.navigation(list(_PAGES.values()))
pg.run()
