"""Tests for jnuherb.network."""

import pandas as pd
import pytest
from jnuherb.network import (
    calculate_degree,
    degree_summary,
    filter_formula_network,
    filter_target_network,
    network_analysis_pipeline,
    plot_degree_barchart,
    plot_degree_boxplot,
)


@pytest.fixture
def sample_network():
    """Create a small sample network for testing."""
    edges = pd.DataFrame(
        {
            "SourceNode": [
                "TCM1", "TCM1", "TCM2",
                "Chem1", "Chem2", "Chem2",
            ],
            "TargetNode": [
                "Chem1", "Chem2", "Chem2",
                "Prot1", "Prot1", "Prot2",
            ],
        }
    )
    nodes = pd.DataFrame(
        {
            "Key": [
                "TCM1", "TCM2", "Chem1", "Chem2", "Prot1", "Prot2",
            ],
            "Attribute": [
                "TCM", "TCM", "Chemicals", "Chemicals", "Proteins", "Proteins",
            ],
        }
    )
    return edges, nodes


def test_calculate_degree(sample_network):
    edges, nodes = sample_network
    result = calculate_degree(edges, nodes)
    assert "Degree" in result.columns
    assert len(result) == len(nodes)
    # TCM1 has 2 edges (to Chem1, Chem2)
    tcm1_deg = result.loc[result["Key"] == "TCM1", "Degree"].values[0]
    assert tcm1_deg == 2


def test_calculate_degree_isolated_nodes(sample_network):
    edges, nodes = sample_network
    # Add an isolated node not in edges
    isolated = pd.DataFrame(
        {"Key": ["Isolated"], "Attribute": ["TCM"]}
    )
    nodes_ext = pd.concat([nodes, isolated], ignore_index=True)
    result = calculate_degree(edges, nodes_ext)
    iso_deg = result.loc[result["Key"] == "Isolated", "Degree"].values[0]
    assert iso_deg == 0


def test_degree_summary(sample_network):
    edges, nodes = sample_network
    node_deg = calculate_degree(edges, nodes)
    summary = degree_summary(node_deg)
    assert "Attribute" in summary.columns
    assert "Mean_Degree" in summary.columns
    assert "Max_Degree" in summary.columns
    assert "Node_Count" in summary.columns
    assert len(summary) == 3  # TCM, Chemicals, Proteins


def test_degree_summary_missing_columns():
    """Should raise with helpful error messages."""
    bad_df = pd.DataFrame({"wrong": [1, 2]})
    with pytest.raises(ValueError, match="Attribute"):
        degree_summary(bad_df)

    bad_df2 = pd.DataFrame({"Attribute": ["TCM"], "wrong": [1]})
    with pytest.raises(ValueError, match="Degree"):
        degree_summary(bad_df2)


def test_filter_formula_network(sample_network):
    edges, nodes = sample_network
    result = filter_formula_network(
        edges, nodes, ["TCM1", "TCM2"],
        include_proteins=True,
    )
    assert len(result["tcms_found"]) == 2
    assert len(result["tcms_missing"]) == 0
    assert len(result["network"]) > 0
    assert len(result["node"]) > 0


def test_filter_formula_network_missing():
    """Should warn when no TCMs found."""
    edges = pd.DataFrame(
        {"SourceNode": ["A"], "TargetNode": ["B"]}
    )
    nodes = pd.DataFrame(
        {"Key": ["A", "B"], "Attribute": ["TCM", "Chemicals"]}
    )
    with pytest.warns(UserWarning, match="None of the specified"):
        result = filter_formula_network(edges, nodes, ["NonExistent"])
    assert len(result["tcms_found"]) == 0


def test_filter_target_network(sample_network):
    edges, nodes = sample_network
    result = filter_target_network(
        edges, nodes, ["Prot1", "Prot2"],
        upstream=True,
    )
    assert len(result["proteins_found"]) == 2
    assert len(result["network"]) > 0


def test_plot_degree_boxplot(sample_network):
    edges, nodes = sample_network
    node_deg = calculate_degree(edges, nodes)
    # Should not raise
    plot_degree_boxplot(node_deg, export_path="/tmp/jnuherb_test")


def test_plot_degree_barchart(sample_network):
    edges, nodes = sample_network
    node_deg = calculate_degree(edges, nodes)
    result = plot_degree_barchart(node_deg, top_n=5, export_path="/tmp/jnuherb_test")
    assert result is not None
    assert len(result) > 0


def test_network_analysis_pipeline(sample_network):
    edges, nodes = sample_network
    result = network_analysis_pipeline(
        edges, nodes,
        export_path="/tmp/jnuherb_test",
        prefix="test",
    )
    assert "node_degree" in result
    assert "summary" in result
