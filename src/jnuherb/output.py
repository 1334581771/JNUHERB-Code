"""Cytoscape export and D3.js HTML network visualization.

* ``out_for_cyto`` — writes ``Network.csv`` (edges) and ``Type.csv``
  (node attributes) for import into Cytoscape.
* ``vis_network`` — generates a self-contained interactive HTML
  force-directed network graph using D3.js v7.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

import pandas as pd


# ---------------------------------------------------------------------------
# _re_name — internal ID→name mapping
# ---------------------------------------------------------------------------


def _re_name(
    tcm: pd.DataFrame,
    tcm_chem_links: pd.DataFrame,
    chem: pd.DataFrame,
    chem_protein_links: pd.DataFrame,
    protein: pd.DataFrame,
) -> dict[str, pd.DataFrame]:
    """Rename internal IDs to human-readable names for network output.

    Returns a dict with keys ``tcm``, ``tcm_chem``, ``chem``,
    ``chem_protein``, ``gene``.
    """
    tcm_c = tcm.copy()
    tcm_chem_links_c = tcm_chem_links.copy()
    chem_c = chem.copy()
    chem_protein_links_c = chem_protein_links.copy()
    protein_c = protein.copy()

    # Build lookup dicts for fast mapping
    chem_name_map = dict(zip(chem_c["HVCID"], chem_c["Name"]))
    prot_name_map = dict(zip(protein_c["Ensembl_ID"], protein_c["gene_name"]))
    tcm_name_map = dict(zip(tcm_c["HVMID"], tcm_c["cn_name"]))

    # chemical → protein edges
    out_chem_prot = chem_protein_links_c.iloc[:, :2].copy()
    out_chem_prot.columns = ["SourceNode", "TargetNode"]
    out_chem_prot["SourceNode"] = out_chem_prot["SourceNode"].map(chem_name_map)
    out_chem_prot = out_chem_prot.dropna(subset=["SourceNode"])
    out_chem_prot["TargetNode"] = out_chem_prot["TargetNode"].map(prot_name_map)
    out_chem_prot = out_chem_prot.dropna(subset=["TargetNode"])

    # TCM → chemical edges
    out_tcm_chem = tcm_chem_links_c.iloc[:, :2].copy()
    out_tcm_chem.columns = ["SourceNode", "TargetNode"]
    out_tcm_chem["SourceNode"] = out_tcm_chem["SourceNode"].map(tcm_name_map)
    out_tcm_chem = out_tcm_chem.dropna(subset=["SourceNode"])
    out_tcm_chem["TargetNode"] = out_tcm_chem["TargetNode"].map(chem_name_map)
    out_tcm_chem = out_tcm_chem.dropna(subset=["TargetNode"])

    # Node type DataFrames
    out_tcm = pd.DataFrame(
        {"Key": tcm_c["cn_name"], "Attribute": "TCM"}
    )
    out_chem = pd.DataFrame(
        {"Key": chem_c["Name"], "Attribute": "Chemicals"}
    )
    out_gene = pd.DataFrame(
        {"Key": protein_c["gene_name"], "Attribute": "Proteins"}
    )

    return {
        "tcm": out_tcm,
        "tcm_chem": out_tcm_chem,
        "chem": out_chem,
        "chem_protein": out_chem_prot,
        "gene": out_gene,
    }


# ---------------------------------------------------------------------------
# out_for_cyto
# ---------------------------------------------------------------------------


def out_for_cyto(
    tcm: pd.DataFrame,
    tcm_chem_links: pd.DataFrame,
    chem: pd.DataFrame,
    chem_protein_links: pd.DataFrame,
    protein: pd.DataFrame,
    path: str = "result",
) -> None:
    """Export Cytoscape-compatible CSV files.

    Writes two files into ``path``:

    * ``Network.csv`` — edge list with ``SourceNode``, ``TargetNode``.
    * ``Type.csv`` — node attributes with ``Key``, ``Attribute``
      (one of TCM / Chemicals / Proteins).

    Args:
        tcm: TCM information.
        tcm_chem_links: TCM-chemical links.
        chem: Chemical information.
        chem_protein_links: Chemical-protein links.
        protein: Protein information.
        path: Output directory (created if it does not exist).
    """
    out_dir = Path(path)
    out_dir.mkdir(parents=True, exist_ok=True)

    renamed = _re_name(tcm, tcm_chem_links, chem, chem_protein_links, protein)

    # Edges
    network = pd.concat(
        [renamed["chem_protein"], renamed["tcm_chem"]], ignore_index=True
    )
    network.to_csv(out_dir / "Network.csv", index=False, encoding="utf-8")

    # Nodes
    type_df = pd.concat(
        [renamed["tcm"], renamed["chem"], renamed["gene"]], ignore_index=True
    )
    type_df.to_csv(out_dir / "Type.csv", index=False, encoding="utf-8")


# ---------------------------------------------------------------------------
# vis_network  (interactive D3.js HTML)
# ---------------------------------------------------------------------------


def vis_network(
    tcm: pd.DataFrame,
    tcm_chem_links: pd.DataFrame,
    chem: pd.DataFrame,
    chem_protein_links: pd.DataFrame,
    protein: pd.DataFrame,
    path: str = "result",
) -> None:
    """Generate an interactive D3.js force-directed network as HTML.

    Writes ``Graph.html`` into ``path``.  The HTML is self-contained
    (loads D3.js v7 from CDN) and supports drag, pan, and zoom.

    Args:
        tcm: TCM information.
        tcm_chem_links: TCM-chemical links.
        chem: Chemical information.
        chem_protein_links: Chemical-protein links.
        protein: Protein information.
        path: Output directory (created if it does not exist).
    """
    out_dir = Path(path)
    out_dir.mkdir(parents=True, exist_ok=True)

    renamed = _re_name(tcm, tcm_chem_links, chem, chem_protein_links, protein)

    # Edges
    edges_df = pd.concat(
        [renamed["tcm_chem"], renamed["chem_protein"]], ignore_index=True
    )
    edges_df.columns = ["from", "to"]

    # Nodes
    nodes_tcm = pd.DataFrame(
        {
            "name": renamed["tcm"]["Key"],
            "type": "TCM",
            "color": "#1FA9E9",
            "size": 20,
        }
    )
    nodes_chem = pd.DataFrame(
        {
            "name": renamed["chem"]["Key"],
            "type": "Chemical",
            "color": "#FFFF00",
            "size": 20,
        }
    )
    nodes_prot = pd.DataFrame(
        {
            "name": renamed["gene"]["Key"],
            "type": "Protein",
            "color": "#000000",
            "size": 20,
        }
    )

    nodes = pd.concat(
        [nodes_tcm, nodes_chem, nodes_prot], ignore_index=True
    )
    nodes = nodes.drop_duplicates(subset=["name"])

    if len(edges_df) > 0:
        _create_network_html(nodes, edges_df, str(out_dir / "Graph.html"))


# ---------------------------------------------------------------------------
# _create_network_html  (internal)
# ---------------------------------------------------------------------------


def _create_network_html(
    nodes: pd.DataFrame,
    edges: pd.DataFrame,
    output_path: str,
) -> None:
    """Write a self-contained D3.js network visualization HTML file.

    Args:
        nodes: DataFrame with columns ``name``, ``type``, ``color``, ``size``.
        edges: DataFrame with columns ``from``, ``to``.
        output_path: Destination file path.
    """
    # Build JSON arrays
    node_records = []
    for _, row in nodes.iterrows():
        node_records.append(
            {
                "name": row["name"],
                "group": row["type"],
                "color": row["color"],
                "size": int(row["size"]),
            }
        )
    nodes_json = json.dumps(node_records, ensure_ascii=False)

    edge_records = []
    for _, row in edges.iterrows():
        edge_records.append({"source": row["from"], "target": row["to"]})
    edges_json = json.dumps(edge_records, ensure_ascii=False)

    html = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>HerbiV Network Visualization</title>
<script src="https://d3js.org/d3.v7.min.js"></script>
<style>
  body {{ margin: 0; padding: 0; font-family: Arial, sans-serif; }}
  #network {{ width: 100%; height: 100vh; }}
  .node {{ stroke: #fff; stroke-width: 1.5px; }}
  .link {{ stroke: #999; stroke-opacity: 0.6; }}
  text {{ font-family: Arial, sans-serif; font-size: 10px; }}
</style>
</head>
<body>
<div id="network"></div>
<script>
var nodes = {nodes_json};
var links = {edges_json};

var width = window.innerWidth;
var height = window.innerHeight;

var svg = d3.select("#network").append("svg")
    .attr("width", width)
    .attr("height", height);

var colorMap = {{
  "TCM": "#1FA9E9",
  "Chemical": "#FFFF00",
  "Protein": "#000000"
}};

var simulation = d3.forceSimulation(nodes)
    .force("link", d3.forceLink(links).id(function(d) {{ return d.name; }}).distance(100))
    .force("charge", d3.forceManyBody().strength(-8000))
    .force("center", d3.forceCenter(width / 2, height / 2));

var link = svg.append("g")
    .selectAll("line")
    .data(links)
    .enter().append("line")
    .attr("class", "link")
    .attr("stroke-width", 1);

var node = svg.append("g")
    .selectAll("g")
    .data(nodes)
    .enter().append("g")
    .call(d3.drag()
        .on("start", dragstarted)
        .on("drag", dragged)
        .on("end", dragended));

node.append("circle")
    .attr("class", "node")
    .attr("r", function(d) {{ return d.size / 2; }})
    .attr("fill", function(d) {{ return d.color; }});

node.append("text")
    .attr("dx", 12)
    .attr("dy", ".35em")
    .text(function(d) {{ return d.name; }});

simulation.on("tick", function() {{
    link
        .attr("x1", function(d) {{ return d.source.x; }})
        .attr("y1", function(d) {{ return d.source.y; }})
        .attr("x2", function(d) {{ return d.target.x; }})
        .attr("y2", function(d) {{ return d.target.y; }});

    node.attr("transform", function(d) {{
        return "translate(" + d.x + "," + d.y + ")";
    }});
}});

function dragstarted(event, d) {{
    if (!event.active) simulation.alphaTarget(0.3).restart();
    d.fx = d.x;
    d.fy = d.y;
}}

function dragged(event, d) {{
    d.fx = event.x;
    d.fy = event.y;
}}

function dragended(event, d) {{
    if (!event.active) simulation.alphaTarget(0);
    d.fx = null;
    d.fy = null;
}}
</script>
</body>
</html>"""

    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html)
