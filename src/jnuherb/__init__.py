"""JNUHERB — A Network Pharmacology Analysis Tool for Traditional Chinese Medicine.

Provides forward and reverse network pharmacology analysis integrating
over 578,000 records of formulas, TCMs, chemicals, and protein targets
with HerbiV scoring and combination optimization.

Quick start
-----------
>>> import jnuherb

>>> # Reverse analysis: proteins → chemicals → TCMs → formulas
>>> result = jnuherb.from_proteins(
...     ["ENSP00000381588", "ENSP00000252519"],
...     score=0, random_state=138192, num=100,
... )
>>> print(result["tcm"][["cn_name", "Importance Score"]].head())

>>> # Forward analysis: TCM → chemicals → proteins
>>> result = jnuherb.from_tcm_or_formula(
...     ["HVM0367", "HVM1695"],
...     score=900,
... )

>>> # One-step integrated analysis
>>> result = jnuherb.herb_analysis(
...     proteins=["PDCD1", "TP53"],
...     protein_id_type="gene_name",
...     formula=["柴胡", "黄芩"],
...     path="results",
... )

>>> # Launch the web GUI
>>> jnuherb.launch_app()
"""

__version__ = "0.1.0"

# Query functions
from jnuherb.query import (
    get_chem_protein_links,
    get_chemicals,
    get_formula,
    get_formula_tcm_links,
    get_proteins,
    get_tcm,
    get_tcm_chem_links,
)

# Analysis functions
from jnuherb.analysis import (
    dfs_filter,
    from_proteins,
    from_tcm_or_formula,
)

# Scoring functions
from jnuherb.score import (
    component,
    herbiv_score,
)

# Network functions
from jnuherb.network import (
    calculate_degree,
    degree_summary,
    filter_formula_network,
    filter_target_network,
    network_analysis_pipeline,
    plot_degree_barchart,
    plot_degree_boxplot,
)

# Output functions
from jnuherb.output import (
    out_for_cyto,
    vis_network,
)

# Pipeline
from jnuherb.pipeline import (
    herb_analysis,
)

# App
from jnuherb.app import (
    launch_app,
)

__all__ = [
    # query
    "get_formula",
    "get_formula_tcm_links",
    "get_tcm",
    "get_tcm_chem_links",
    "get_chemicals",
    "get_chem_protein_links",
    "get_proteins",
    # analysis
    "from_tcm_or_formula",
    "from_proteins",
    "dfs_filter",
    # score
    "herbiv_score",
    "component",
    # network
    "filter_formula_network",
    "filter_target_network",
    "calculate_degree",
    "degree_summary",
    "plot_degree_boxplot",
    "plot_degree_barchart",
    "network_analysis_pipeline",
    # output
    "out_for_cyto",
    "vis_network",
    # pipeline
    "herb_analysis",
    # app
    "launch_app",
]
