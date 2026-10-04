"""Network filtering, degree centrality, and visualization.

Functions for post-analysis network exploration:

* ``filter_formula_network`` — sub-network for a specific formula.
* ``filter_target_network`` — sub-network for specific target proteins.
* ``calculate_degree`` — node degree centrality.
* ``degree_summary`` — aggregate degree statistics by attribute.
* ``plot_degree_boxplot`` — degree distribution boxplot (matplotlib).
* ``plot_degree_barchart`` — degree bar chart for key nodes (matplotlib).
* ``network_analysis_pipeline`` — run all post-analysis steps at once.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Union

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# Use non-interactive backend for headless environments
matplotlib.use("Agg")

# Register bundled Chinese font so PDF/PNG exports render CJK characters correctly.
try:
    from importlib.resources import files
    from matplotlib import font_manager

    _font_path = files("jnuherb.data.fonts").joinpath("NotoSansCJKsc-Regular.otf")
    if _font_path.exists():
        font_manager.fontManager.addfont(str(_font_path))
        plt.rcParams["font.family"] = ["Noto Sans CJK SC", "sans-serif"]
    else:
        plt.rcParams["font.family"] = ["sans-serif"]
except Exception:
    plt.rcParams["font.family"] = ["sans-serif"]
plt.rcParams["axes.unicode_minus"] = False


# ---------------------------------------------------------------------------
# filter_formula_network
# ---------------------------------------------------------------------------


def filter_formula_network(
    network: pd.DataFrame,
    node: pd.DataFrame,
    formula_names: list[str],
    include_proteins: bool = True,
    path: Optional[str] = None,
    prefix: str = "formula",
) -> dict:
    """Filter network to a specific formula's TCM components.

    Args:
        network: Edge list with ``SourceNode`` and ``TargetNode``.
        node: Node attributes with ``Key`` and ``Attribute``
            (TCM / Chemicals / Proteins).
        formula_names: TCM names (typically Chinese names) to filter by.
        include_proteins: If True, also include downstream protein
            targets connected via chemicals.
        path: Optional directory to export filtered CSV files.
        prefix: Filename prefix for exported CSVs.

    Returns:
        Dict with keys ``network``, ``node``, ``tcms_found``,
        ``tcms_missing``.
    """
    tcms_found = [n for n in formula_names if n in node["Key"].values]
    tcms_missing = [n for n in formula_names if n not in node["Key"].values]

    if not tcms_found:
        import warnings

        warnings.warn("None of the specified TCM names were found in the network.")
        return {
            "network": network.iloc[0:0].copy(),
            "node": node.iloc[0:0].copy(),
            "tcms_found": [],
            "tcms_missing": formula_names,
        }

    # Edges: source in formula TCM list
    network_tcm = network[network["SourceNode"].isin(tcms_found)]

    if include_proteins and len(network_tcm) > 0:
        compounds = network_tcm["TargetNode"].unique()
        network_proteins = network[network["SourceNode"].isin(compounds)]
        final_network = pd.concat(
            [network_tcm, network_proteins], ignore_index=True
        )
    else:
        final_network = network_tcm

    # Extract affected nodes
    nodes_in_network = set(final_network["SourceNode"]) | set(
        final_network["TargetNode"]
    )
    final_node = node[node["Key"].isin(nodes_in_network)]

    # Export
    if path is not None:
        out_dir = Path(path)
        out_dir.mkdir(parents=True, exist_ok=True)
        final_network.to_csv(
            out_dir / f"{prefix}_Network_Edges.csv",
            index=False, encoding="utf-8",
        )
        final_node.to_csv(
            out_dir / f"{prefix}_Network_Nodes.csv",
            index=False, encoding="utf-8",
        )

    return {
        "network": final_network,
        "node": final_node,
        "tcms_found": tcms_found,
        "tcms_missing": tcms_missing,
    }


# ---------------------------------------------------------------------------
# filter_target_network
# ---------------------------------------------------------------------------


def filter_target_network(
    network: pd.DataFrame,
    node: pd.DataFrame,
    target_proteins: list[str],
    upstream: bool = True,
    path: Optional[str] = None,
    prefix: str = "target",
) -> dict:
    """Filter network to specific target proteins.

    Args:
        network: Edge list with ``SourceNode`` and ``TargetNode``.
        node: Node attributes with ``Key`` and ``Attribute``.
        target_proteins: Protein names (gene_name) or Ensembl IDs.
        upstream: If True, also include upstream chemicals.
        path: Optional directory to export filtered CSV files.
        prefix: Filename prefix for exported CSVs.

    Returns:
        Dict with keys ``network``, ``node``, ``proteins_found``,
        ``proteins_missing``.
    """
    proteins_found = [
        p for p in target_proteins if p in node["Key"].values
    ]
    proteins_missing = [
        p for p in target_proteins if p not in node["Key"].values
    ]

    if not proteins_found:
        import warnings

        warnings.warn("None of the specified proteins were found in the network.")
        return {
            "network": network.iloc[0:0].copy(),
            "node": node.iloc[0:0].copy(),
            "proteins_found": [],
            "proteins_missing": target_proteins,
        }

    # Edges: target in protein list (proteins are TargetNodes)
    network_protein = network[network["TargetNode"].isin(proteins_found)]

    if upstream and len(network_protein) > 0:
        chemicals = network_protein["SourceNode"].unique()
        network_chemicals = network[network["TargetNode"].isin(chemicals)]
        final_network = pd.concat(
            [network_chemicals, network_protein], ignore_index=True
        )
    else:
        final_network = network_protein

    nodes_in_network = set(final_network["SourceNode"]) | set(
        final_network["TargetNode"]
    )
    final_node = node[node["Key"].isin(nodes_in_network)]

    if path is not None:
        out_dir = Path(path)
        out_dir.mkdir(parents=True, exist_ok=True)
        final_network.to_csv(
            out_dir / f"{prefix}_Network_Edges.csv",
            index=False, encoding="utf-8",
        )
        final_node.to_csv(
            out_dir / f"{prefix}_Network_Nodes.csv",
            index=False, encoding="utf-8",
        )

    return {
        "network": final_network,
        "node": final_node,
        "proteins_found": proteins_found,
        "proteins_missing": proteins_missing,
    }


# ---------------------------------------------------------------------------
# calculate_degree
# ---------------------------------------------------------------------------


def calculate_degree(
    network: pd.DataFrame,
    node: pd.DataFrame,
    export_path: Optional[str] = None,
    prefix: str = "network",
) -> pd.DataFrame:
    """Calculate degree centrality for every node in the network.

    Args:
        network: Edge list with ``SourceNode`` and ``TargetNode``.
        node: Node attributes with ``Key`` and ``Attribute``.
        export_path: Optional directory to save degree CSV.
        prefix: Filename prefix.

    Returns:
        Node DataFrame with an added ``Degree`` column (0 for isolated
        nodes).
    """
    all_nodes = pd.concat(
        [network["SourceNode"], network["TargetNode"]]
    )
    degree_counts = all_nodes.value_counts().reset_index()
    degree_counts.columns = ["Key", "Degree"]
    degree_counts["Degree"] = degree_counts["Degree"].astype(int)

    node_degree = node.merge(degree_counts, on="Key", how="left")
    node_degree["Degree"] = node_degree["Degree"].fillna(0).astype(int)

    if export_path is not None:
        out_dir = Path(export_path)
        out_dir.mkdir(parents=True, exist_ok=True)
        node_degree.to_csv(
            out_dir / f"{prefix}_degree.csv", index=False, encoding="utf-8",
        )

    return node_degree


# ---------------------------------------------------------------------------
# degree_summary
# ---------------------------------------------------------------------------


def degree_summary(node_degree: pd.DataFrame) -> pd.DataFrame:
    """Summarize degree statistics grouped by node attribute.

    Args:
        node_degree: Output of ``calculate_degree`` — must have
            ``Attribute`` and ``Degree`` columns.

    Returns:
        DataFrame with columns ``Attribute``, ``Mean_Degree``,
        ``Max_Degree``, ``Node_Count``.
    """
    if "Attribute" not in node_degree.columns:
        raise ValueError("node_degree must contain an 'Attribute' column.")
    if "Degree" not in node_degree.columns:
        raise ValueError(
            "node_degree must contain a 'Degree' column. "
            "Run calculate_degree() first."
        )

    summary = (
        node_degree.groupby("Attribute")["Degree"]
        .agg(Mean_Degree="mean", Max_Degree="max", Node_Count="count")
        .reset_index()
    )
    summary["Mean_Degree"] = summary["Mean_Degree"].round(2)
    summary["Max_Degree"] = summary["Max_Degree"].astype(int)
    summary["Node_Count"] = summary["Node_Count"].astype(int)
    return summary


# ---------------------------------------------------------------------------
# plot_degree_boxplot
# ---------------------------------------------------------------------------


def plot_degree_boxplot(
    node_degree: pd.DataFrame,
    colors: Optional[list[str]] = None,
    main: str = "Degree Distribution by Attribute",
    export_path: Optional[str] = None,
    prefix: str = "degree_boxplot",
) -> None:
    """Generate a degree distribution boxplot grouped by attribute.

    Args:
        node_degree: Output of ``calculate_degree``.
        colors: Bar colors per attribute.  Defaults to
            ``["gold", "lightgreen", "coral"]``.
        main: Plot title.
        export_path: Optional directory to save the plot as PDF.
        prefix: Filename prefix.
    """
    if colors is None:
        colors = ["gold", "lightgreen", "coral"]

    n_attrs = node_degree["Attribute"].nunique()
    if n_attrs > len(colors):
        colors = sns.color_palette("husl", n_attrs)
    else:
        colors = colors[:n_attrs]

    fig, ax = plt.subplots(figsize=(8, 6))
    sns.boxplot(
        data=node_degree,
        x="Attribute",
        y="Degree",
        palette=colors,
        ax=ax,
    )
    ax.set_title(main)
    ax.set_xlabel("Node Attribute")
    ax.set_ylabel("Degree")

    if export_path is not None:
        out_dir = Path(export_path)
        out_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(
            out_dir / f"{prefix}.pdf", bbox_inches="tight",
        )
    plt.close(fig)


# ---------------------------------------------------------------------------
# plot_degree_barchart
# ---------------------------------------------------------------------------


def plot_degree_barchart(
    node_degree: pd.DataFrame,
    filter_list: Optional[list[str]] = None,
    top_n: int = 10,
    fill_color: str = "#BBBDE6",
    export_path: Optional[str] = None,
    prefix: str = "Bar",
) -> Optional[dict[str, plt.Figure]]:
    """Generate horizontal bar charts of node degrees by attribute.

    Args:
        node_degree: Output of ``calculate_degree``.
        filter_list: Specific node names to include.  If None, the
            top ``top_n`` nodes per attribute are used.
        top_n: Number of top nodes per attribute when ``filter_list``
            is None.
        fill_color: Bar fill colour.
        export_path: Optional directory to save plots as PDFs.
        prefix: Filename prefix.

    Returns:
        Dict mapping attribute names to matplotlib Figures, or None
        if no nodes match.
    """
    if filter_list is not None:
        plot_data = node_degree[node_degree["Key"].isin(filter_list)]
    else:
        top_parts = []
        for _, grp in node_degree.groupby("Attribute"):
            top_parts.append(grp.nlargest(top_n, "Degree"))
        plot_data = pd.concat(top_parts, ignore_index=True) if top_parts else pd.DataFrame()

    if len(plot_data) == 0:
        import warnings

        warnings.warn("No nodes matched the criteria. Plot not generated.")
        return None

    plot_list: dict[str, plt.Figure] = {}

    for attr, sub_data in plot_data.groupby("Attribute"):
        sub_data = sub_data.sort_values("Degree")  # ascending for horizontal bar

        fig, ax = plt.subplots(figsize=(5, max(3, len(sub_data) * 0.3)))
        bars = ax.barh(sub_data["Key"], sub_data["Degree"], color=fill_color, alpha=0.8)

        # Add value labels
        for bar, val in zip(bars, sub_data["Degree"]):
            ax.text(
                bar.get_width() + bar.get_width() * 0.02,
                bar.get_y() + bar.get_height() / 2,
                str(val),
                va="center",
                fontsize=8,
            )

        ax.set_title(f"Top Candidates - {attr}", fontweight="bold")
        ax.set_xlabel("Degree")
        ax.margins(x=0.15)
        fig.tight_layout()

        if export_path is not None:
            out_dir = Path(export_path)
            out_dir.mkdir(parents=True, exist_ok=True)
            fig.savefig(
                out_dir / f"{prefix}_{attr}.pdf", bbox_inches="tight",
            )

        plot_list[attr] = fig
        plt.close(fig)

    return plot_list


# ---------------------------------------------------------------------------
# network_analysis_pipeline
# ---------------------------------------------------------------------------


def network_analysis_pipeline(
    network: pd.DataFrame,
    node: pd.DataFrame,
    filter_list: Optional[list[str]] = None,
    export_path: str = "results",
    prefix: str = "analysis",
) -> dict:
    """Run a complete post-analysis pipeline.

    Executes ``calculate_degree`` → ``degree_summary`` →
    ``plot_degree_boxplot`` → ``plot_degree_barchart`` in sequence.

    Args:
        network: Edge list (``SourceNode``, ``TargetNode``).
        node: Node attributes (``Key``, ``Attribute``).
        filter_list: Specific nodes to highlight in bar charts.
        export_path: Output directory.
        prefix: Filename prefix.

    Returns:
        Dict with keys ``node_degree`` and ``summary``.
    """
    node_degree = calculate_degree(network, node, export_path, prefix)

    summary_stats = degree_summary(node_degree)
    print("\n=== Degree Summary by Attribute ===")
    print(summary_stats.to_string(index=False))

    plot_degree_boxplot(node_degree, export_path=export_path, prefix=prefix)

    plot_degree_barchart(
        node_degree,
        filter_list=filter_list,
        export_path=export_path,
        prefix=prefix,
    )

    return {"node_degree": node_degree, "summary": summary_stats}
