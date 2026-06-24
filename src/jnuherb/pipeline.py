"""Integrated analysis pipeline — all-in-one workflow.

``herb_analysis`` combines reverse network pharmacology analysis,
formula filtering, degree centrality, and visualization into a
single function call.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import pandas as pd

from jnuherb.analysis import from_proteins
from jnuherb.network import (
    calculate_degree,
    degree_summary,
    filter_formula_network,
    plot_degree_barchart,
    plot_degree_boxplot,
)


def herb_analysis(
    proteins: list[str],
    formula: Optional[list[str]] = None,
    protein_id_type: str = "Ensembl_ID",
    formula_id_type: str = "cn_name",
    score: float = 0,
    path: str = "results",
    out_for_cytoscape: bool = True,
    out_graph: bool = True,
    tcm_component: bool = True,
    formula_component: bool = True,
    random_state: Optional[int] = None,
    num: int = 1000,
    degree_boxplot: bool = True,
    key_nodes: Optional[list[str]] = None,
) -> dict:
    """Integrated herb-target network pharmacology analysis.

    Runs reverse analysis from target proteins, optionally filters by
    a formula, calculates degree centrality, and generates
    visualizations — all in one call.

    Args:
        proteins: Target protein IDs (Ensembl_ID or gene_name).
        formula: Optional TCM names (Chinese names by default) to
            filter the network.
        protein_id_type: ``"Ensembl_ID"`` (default) or ``"gene_name"``.
        formula_id_type: For ``formula``: ``"cn_name"`` (default),
            ``"HVMID"``, ``"pinyin_name"``, ``"en_name"``,
            or ``"latin_name"``.
        score: Minimum ``Combined_score`` threshold (default 0).
        path: Output directory.
        out_for_cytoscape: Whether to write Cytoscape CSV files.
        out_graph: Whether to generate interactive HTML graph.
        tcm_component: Whether to run TCM combination optimization.
        formula_component: Whether to run formula combination optimization.
        random_state: Random seed for optimization.
        num: Number of optimization solution sets.
        degree_boxplot: Whether to generate a degree distribution boxplot.
        key_nodes: Node names to highlight in bar charts.  If None,
            auto-selects from formula components and top targets.

    Returns:
        Dict with keys ``protein_result``, ``network``, ``node``,
        ``formula_filter``, ``node_degree``, ``degree_summary``.
    """
    # ---- Step 1: Reverse analysis from targets ----
    print("Step 1: Running reverse network pharmacology analysis from targets...")

    protein_result = from_proteins(
        proteins=proteins,
        score=score,
        random_state=random_state,
        num=num,
        tcm_component=tcm_component,
        formula_component=formula_component,
        out_for_cytoscape=out_for_cytoscape,
        re=True,
        path=path,
        id_type=protein_id_type,
    )

    print(f"  - TCMs found: {len(protein_result['tcm'])}")
    print(f"  - Chemicals found: {len(protein_result['chem'])}")
    print(f"  - Proteins found: {len(protein_result['proteins'])}")

    # ---- Step 2: Read or build network ----
    network_path = Path(path) / "Network.csv"
    node_path = Path(path) / "Type.csv"

    if network_path.exists() and node_path.exists():
        network = pd.read_csv(str(network_path), dtype=str, keep_default_na=False)
        node = pd.read_csv(str(node_path), dtype=str, keep_default_na=False)
    else:
        if out_for_cytoscape:
            raise FileNotFoundError(
                "Cytoscape output files not found. Set out_for_cytoscape=True."
            )
        # Build from result data
        network = pd.DataFrame(
            {
                "SourceNode": list(protein_result["tcm_chem_links"].iloc[:, 0])
                + list(protein_result["chem_protein_links"].iloc[:, 0]),
                "TargetNode": list(protein_result["tcm_chem_links"].iloc[:, 1])
                + list(protein_result["chem_protein_links"].iloc[:, 1]),
            }
        )
        node = pd.DataFrame(
            {
                "Key": (
                    list(protein_result["tcm"]["cn_name"])
                    + list(protein_result["chem"]["Name"])
                    + list(protein_result["proteins"]["gene_name"])
                ),
                "Attribute": (
                    ["TCM"] * len(protein_result["tcm"])
                    + ["Chemicals"] * len(protein_result["chem"])
                    + ["Proteins"] * len(protein_result["proteins"])
                ),
            }
        )

    # ---- Step 3: Formula filter (optional) ----
    formula_filter = None

    if formula is not None:
        print("Step 2: Filtering network for specified formula...")

        found = [f for f in formula if f in node["Key"].values]
        missing = [f for f in formula if f not in node["Key"].values]

        print(f"  - Formula components found in network: {len(found)}/{len(formula)}")
        if missing:
            print(f"  - Not found: {', '.join(missing)}")

        formula_filter = filter_formula_network(
            network=network,
            node=node,
            formula_names=formula,
            include_proteins=True,
            path=path,
            prefix="formula",
        )

        print(f"  - Filtered network edges: {len(formula_filter['network'])}")
        print(f"  - Filtered network nodes: {len(formula_filter['node'])}")

    # ---- Step 4: Degree centrality ----
    print("Step 3: Calculating degree centrality...")

    node_degree = calculate_degree(
        network=network,
        node=node,
        export_path=path,
        prefix="network",
    )

    deg_summary = degree_summary(node_degree)

    # ---- Step 5: Visualizations ----
    if degree_boxplot:
        print("Step 4: Generating degree distribution boxplot...")
        plot_degree_boxplot(
            node_degree=node_degree,
            main="Degree Distribution by Attribute",
            export_path=path,
            prefix="degree_boxplot",
        )

    # Auto-select key_nodes
    if formula is not None and key_nodes is None:
        formula_in_network = [f for f in formula if f in node["Key"].values]
        protein_nodes = node_degree[node_degree["Attribute"] == "Proteins"]
        top_proteins = (
            protein_nodes.nlargest(10, "Degree")["Key"].tolist()
            if len(protein_nodes) > 0
            else []
        )
        chem_nodes = node_degree[node_degree["Attribute"] == "Chemicals"]
        top_chems = (
            chem_nodes.nlargest(5, "Degree")["Key"].tolist()
            if len(chem_nodes) > 0
            else []
        )
        key_nodes = formula_in_network + top_chems + top_proteins
        key_nodes = [k for k in key_nodes if k in node["Key"].values]

    if key_nodes:
        print("Step 5: Generating bar charts for key nodes...")
        try:
            plot_degree_barchart(
                node_degree=node_degree,
                filter_list=key_nodes,
                export_path=path,
                prefix="key_nodes",
            )
        except Exception as e:
            print(f"  Bar chart skipped: {e}")

    # ---- Step 6: Summary ----
    print("=== Analysis Complete ===")
    print(f"Results saved to: {path}")

    result = {
        "protein_result": protein_result,
        "network": network,
        "node": node,
        "formula_filter": formula_filter,
        "node_degree": node_degree,
        "degree_summary": deg_summary,
    }

    return result
