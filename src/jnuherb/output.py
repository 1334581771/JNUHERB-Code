"""Cytoscape export and D3.js HTML network visualization.

* ``out_for_cyto`` — writes ``Network.csv`` (edges) and ``Type.csv``
  (node attributes) for import into Cytoscape.
* ``vis_network`` — generates a self-contained interactive HTML
  network graph using D3.js v7 (cluster layout: TCM / Chemicals rings,
  Proteins rectangle; node color and size encode degree).
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
    chem_name_map = dict(zip(chem_c["PubChem_id"], chem_c["Name"]))
    prot_name_map = dict(zip(protein_c["Ensembl_ID"], protein_c["gene_name"]))
    tcm_name_map = dict(zip(tcm_c["HVMID"], tcm_c["cn_name"]))

    # Select edge endpoint columns BY NAME — the chemical-protein link table
    # is ordered (Ensembl_ID, Combined_score, PubChem_id), so positional
    # slicing would grab the score column and silently drop every edge.
    chem_id_col = "PubChem_id" if "PubChem_id" in chem_protein_links_c.columns else "HVCID"
    out_chem_prot = pd.DataFrame({
        "SourceNode": chem_protein_links_c[chem_id_col],
        "TargetNode": chem_protein_links_c["Ensembl_ID"],
    })
    out_chem_prot["SourceNode"] = out_chem_prot["SourceNode"].map(chem_name_map)
    out_chem_prot = out_chem_prot.dropna(subset=["SourceNode"])
    out_chem_prot["TargetNode"] = out_chem_prot["TargetNode"].map(prot_name_map)
    out_chem_prot = out_chem_prot.dropna(subset=["TargetNode"])

    # TCM → chemical edges
    tcl_chem_col = "PubChem_id" if "PubChem_id" in tcm_chem_links_c.columns else "HVCID"
    out_tcm_chem = pd.DataFrame({
        "SourceNode": tcm_chem_links_c["HVMID"],
        "TargetNode": tcm_chem_links_c[tcl_chem_col],
    })
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


def _cluster_layout(
    node_records: list[dict],
) -> dict:
    """Assign fixed cluster positions (in place): TCM nodes form a ring
    on the left, Chemicals a ring in the center, and Proteins a
    rectangle grid on the right.  Returns the dashed-outline frame
    geometry for the three clusters (single source of truth — the HTML
    template draws frames from this dict)."""
    import math

    def _ring(records, cx, cy):
        n = len(records)
        if n == 0:
            return 0.0
        # Radius grows with node count so labels stay readable.
        r = max(90.0, 26.0 * n / math.pi)
        for i, rec in enumerate(records):
            ang = -math.pi / 2 + 2 * math.pi * i / n
            rec["x"] = cx + r * math.cos(ang)
            rec["y"] = cy + r * math.sin(ang)
        return r

    def _rect(records, cx, cy):
        n = len(records)
        if n == 0:
            return {"cols": 1, "rows": 0, "cx": cx, "cy": cy}
        cols = max(1, math.ceil(math.sqrt(n * 1.6)))
        rows = math.ceil(n / cols)
        gw, gh = 118.0, 56.0
        for i, rec in enumerate(records):
            col, row = i % cols, i // cols
            rec["x"] = cx + (col - (cols - 1) / 2) * gw
            rec["y"] = cy + (row - (rows - 1) / 2) * gh
        return {"cols": cols, "rows": rows, "cx": cx, "cy": cy,
                "gw": gw, "gh": gh}

    tcm = [r for r in node_records if r["group"] == "TCM"]
    chem = [r for r in node_records if r["group"] == "Chemical"]
    prot = [r for r in node_records if r["group"] == "Protein"]

    tcm_r = _ring(tcm, 250, 470)
    chem_r = _ring(chem, 810, 470)
    # Keep the protein rectangle clear of the chemical ring: its left
    # edge must sit well right of the ring's dashed frame.
    prot_cols = max(1, math.ceil(math.sqrt(len(prot) * 1.6))) if prot else 1
    prot_cx = max(1360.0, 810 + chem_r + 104 + (prot_cols * 118) / 2)
    pr = _rect(prot, prot_cx, 470)

    return {
        "tcm": {"cx": 250, "cy": 470, "r": tcm_r + 34},
        "chem": {"cx": 810, "cy": 470, "r": chem_r + 34},
        "prot": {
            "x": pr["cx"] - (pr["cols"] * pr["gw"]) / 2 - 30,
            "y": pr["cy"] - (pr["rows"] * pr["gh"]) / 2 - 30,
            "w": pr["cols"] * pr["gw"] + 60,
            "h": max(pr["rows"], 1) * pr["gh"] + 60,
        },
    }


def _create_network_html(
    nodes: pd.DataFrame,
    edges: pd.DataFrame,
    output_path: str,
) -> None:
    """Write a self-contained D3.js network visualization HTML file.

    Layout: TCM nodes form a circle on the left, Chemicals a circle in
    the center, and Proteins a rectangle grid on the right.  Node color
    (light → dark red, YlOrRd) and size both encode the node's degree,
    so hub herbs / ingredients / targets stand out at a glance.  Nodes
    are draggable; positions are deterministic, not force-directed.
    """
    # Build JSON arrays
    node_records = []
    for _, row in nodes.iterrows():
        node_records.append(
            {
                "name": row["name"],
                "group": row["type"],
            }
        )

    # Degree from the edge list (undirected count across both endpoints)
    degree: dict[str, int] = {}
    for _, row in edges.iterrows():
        degree[row["from"]] = degree.get(row["from"], 0) + 1
        degree[row["to"]] = degree.get(row["to"], 0) + 1
    for rec in node_records:
        rec["degree"] = degree.get(rec["name"], 0)

    # Fixed cluster geometry (canvas size derived from content below)
    frames = _cluster_layout(node_records)
    frames_json = json.dumps(frames)

    nodes_json = json.dumps(node_records, ensure_ascii=False)

    edge_records = []
    for _, row in edges.iterrows():
        edge_records.append({"source": row["from"], "target": row["to"]})
    edges_json = json.dumps(edge_records, ensure_ascii=False)

    # Canvas must fit every cluster frame plus label breathing room;
    # large formulas grow the chemical ring, so derive W/H from content.
    max_x = max(
        frames["tcm"]["cx"] + frames["tcm"]["r"],
        frames["chem"]["cx"] + frames["chem"]["r"],
        frames["prot"]["x"] + frames["prot"]["w"],
    )
    min_y = min(
        frames["tcm"]["cy"] - frames["tcm"]["r"],
        frames["chem"]["cy"] - frames["chem"]["r"],
        frames["prot"]["y"],
    )
    max_y = max(
        frames["tcm"]["cy"] + frames["tcm"]["r"],
        frames["chem"]["cy"] + frames["chem"]["r"],
        frames["prot"]["y"] + frames["prot"]["h"],
    )
    canvas_w = int(max(1680, max_x + 280))   # right margin for name labels
    canvas_h = int(max(940, max_y + 90, 940 - min_y + 120))

    html = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>JNUHERB Network</title>
<script src="https://d3js.org/d3.v7.min.js"></script>
<style>
  body {{ margin: 0; padding: 0; font-family: "Segoe UI", Arial, sans-serif;
         background: #fafbfc; }}
  #network {{ width: 100vw; height: 100vh; }}
  svg {{ display: block; margin: 0 auto; }}
  .node {{ stroke: #fff; stroke-width: 1.2px; cursor: grab; }}
  .node:active {{ cursor: grabbing; }}
  .link {{ stroke: #b8c2cc; stroke-opacity: 0.55; fill: none; }}
  .frame {{ fill: none; stroke: #c9d4de; stroke-dasharray: 6 5; }}
  .gtitle {{ font-size: 17px; font-weight: 700; fill: #46525e; }}
  text.nlabel {{ font-family: "Segoe UI", Arial, sans-serif; font-size: 10px;
                 fill: #33404c; pointer-events: none; }}
</style>
</head>
<body>
<div id="network"></div>
<script>
var nodes = {nodes_json};
var links = {edges_json};
var frames = {frames_json};

var W = {canvas_w}, H = {canvas_h};
var svg = d3.select("#network").append("svg")
    .attr("viewBox", "0 0 " + W + " " + H)
    .attr("preserveAspectRatio", "xMidYMid meet");

var maxDeg = d3.max(nodes, function(d) {{ return d.degree; }}) || 1;
var color = d3.scaleSequential(d3.interpolateYlOrRd)
    .domain([0, maxDeg]);
var radius = d3.scaleSqrt()
    .domain([0, maxDeg]).range([7, 20]);

// ---- cluster frames and titles ----
svg.append("circle").attr("class", "frame")
    .attr("cx", frames.tcm.cx).attr("cy", frames.tcm.cy)
    .attr("r", frames.tcm.r);
svg.append("circle").attr("class", "frame")
    .attr("cx", frames.chem.cx).attr("cy", frames.chem.cy)
    .attr("r", frames.chem.r);
svg.append("rect").attr("class", "frame")
    .attr("x", frames.prot.x).attr("y", frames.prot.y)
    .attr("width", frames.prot.w).attr("height", frames.prot.h)
    .attr("rx", 14);
svg.append("text").attr("class", "gtitle")
    .attr("x", frames.tcm.cx).attr("y", frames.tcm.cy - frames.tcm.r - 16)
    .attr("text-anchor", "middle").text("TCM");
svg.append("text").attr("class", "gtitle")
    .attr("x", frames.chem.cx).attr("y", frames.chem.cy - frames.chem.r - 16)
    .attr("text-anchor", "middle").text("Chemicals");
svg.append("text").attr("class", "gtitle")
    .attr("x", frames.prot.x + frames.prot.w / 2)
    .attr("y", frames.prot.y - 16)
    .attr("text-anchor", "middle").text("Proteins");

// ---- degree legend (gradient bar) ----
var lg = svg.append("g")
    .attr("transform", "translate(60," + (H - 70) + ")");
var grad = lg.append("defs").append("linearGradient")
    .attr("id", "deggrad").attr("x1", "0%").attr("x2", "100%");
grad.append("stop").attr("offset", "0%")
    .attr("stop-color", color(0));
grad.append("stop").attr("offset", "100%")
    .attr("stop-color", color(maxDeg));
lg.append("rect").attr("width", 180).attr("height", 12)
    .attr("fill", "url(#deggrad)").attr("rx", 6);
lg.append("text").attr("class", "nlabel").attr("y", 28)
    .text("Degree: 0");
lg.append("text").attr("class", "nlabel").attr("x", 180).attr("y", 28)
    .attr("text-anchor", "end").text("Degree: " + maxDeg);

// ---- edges (curved) ----
function linkPath(d) {{
    var sx = d.source.x, sy = d.source.y, tx = d.target.x, ty = d.target.y;
    var mx = (sx + tx) / 2, my = (sy + ty) / 2;
    var dx = tx - sx, dy = ty - sy;
    // perpendicular offset for a gentle curve
    var off = Math.min(60, Math.hypot(dx, dy) / 6);
    return "M" + sx + "," + sy + " Q" + (mx - dy / 8) + "," + (my + dx / 8)
         + " " + tx + "," + ty;
}}

var link = svg.append("g").selectAll("path").data(links)
    .enter().append("path").attr("class", "link");

// ---- nodes ----
var node = svg.append("g").selectAll("g").data(nodes)
    .enter().append("g")
    .attr("transform", function(d) {{ return "translate(" + d.x + "," + d.y + ")"; }})
    .call(d3.drag()
        .on("drag", function(event, d) {{
            d.x += event.dx; d.y += event.dy;
            d3.select(this).attr("transform", "translate(" + d.x + "," + d.y + ")");
            link.attr("d", linkPath);
        }}));

node.append("circle")
    .attr("class", "node")
    .attr("r", function(d) {{ return radius(d.degree); }})
    .attr("fill", function(d) {{ return color(d.degree); }})
    .append("title")
    .text(function(d) {{ return d.name + "  |  degree: " + d.degree; }});

node.append("text")
    .attr("class", "nlabel")
    .attr("dx", function(d) {{ return radius(d.degree) + 3; }})
    .attr("dy", ".35em")
    .text(function(d) {{ return d.name; }});

link.attr("d", linkPath);
</script>
</body>
</html>"""

    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html)
