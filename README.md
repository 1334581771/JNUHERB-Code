# JNUHERB

**A Network Pharmacology Analysis Tool for Traditional Chinese Medicine**

JNUHERB is both a database and a powerful data analysis platform. It integrates over 578,000 records of formulas, traditional Chinese medicines (TCMs), ingredients, and protein targets, along with validated evaluation models for assessing the potential effects of Chinese medicines and their combinations on disease targets.

This Python package is a conversion of the R package `herbtar` (based on the [HerbiV](https://github.com/MLi-lab-Bioinformatics-NJUCM/HerbiV) project), completely rebuilt following Python conventions.

## Installation

```bash
cd /path/to/jnuherb
pip install -e .
```

> **This repository contains source code only.** The data tables are not
> distributed here — see [Data](#data) before running an analysis.

### Dependencies

- pandas >= 1.5
- numpy >= 1.24
- matplotlib >= 3.7
- seaborn >= 0.12
- plotly >= 5.14
- streamlit >= 1.28 (for web GUI)

## Quick Start

```python
import jnuherb

# Reverse analysis: proteins -> chemicals -> TCMs -> formulas
result = jnuherb.from_proteins(
    ["ENSP00000381588", "ENSP00000252519"],
    score=0, random_state=138192, num=100,
)
print(result["tcm"][["cn_name", "Importance Score"]].head())

# Forward analysis: TCM -> chemicals -> proteins
result = jnuherb.from_tcm_or_formula(
    ["HVM0367", "HVM1695"],
    score=900,
)

# One-step integrated analysis
result = jnuherb.herb_analysis(
    proteins=["PDCD1", "TP53"],
    protein_id_type="gene_name",
    formula=["柴胡", "黄芩"],
    path="results",
)
```

## Web GUI (Streamlit)

```python
import jnuherb
jnuherb.launch_app()  # Opens http://127.0.0.1:8501
```

The Streamlit app provides:
- **Overview**: Analysis summary and metrics
- **AI Consultation**: DeepSeek-powered TCM syndrome differentiation with knowledge-graph validation
- **TCM Results**: Interactive table and score distribution plots
- **Chemical Results**: Chemical ingredients with Importance Scores
- **Protein Results**: Target proteins and degree statistics
- **Network Viz**: Interactive drag-and-drop network visualization
- **Formula Filter**: Filter results to your formula of interest
- **Degree Stats**: Degree distribution boxplots and bar charts
- **DO Enrichment**: Disease Ontology enrichment of protein targets
- **KEGG Enrichment**: KEGG pathway enrichment of protein targets
- **GO Enrichment**: GO BP/CC/MF enrichment of protein targets
- **Optimization**: Optimal TCM/formula combination suggestions

## Features

- **Forward Network Pharmacology Analysis**: From TCMs/formulas to chemicals to protein targets
- **Reverse Network Pharmacology Analysis**: From target proteins to chemicals to TCMs to formulas
- **HerbiV Score Calculation**: Probability-union scoring for evaluating TCM/formula effects
- **Combination Optimization**: Knapsack-based algorithm for optimal TCM/formula combinations
- **Functional Enrichment Analysis**: Offline KEGG/GO/DO enrichment of protein targets with bundled gene sets
- **Network Visualization**: Cytoscape-compatible CSV export and interactive D3.js HTML graphs
- **Built-in Database**: Over 578,000 records across 7 CSV datasets
- **Monomer Query Support**: Query targets for individual chemical compounds (e.g., deoxyshikonin)

## API Reference

### Query Functions

| Function | Description |
|----------|-------------|
| `get_formula(by, items)` | Query formula information |
| `get_formula_tcm_links(by, items)` | Query formula-TCM links |
| `get_tcm(by, items)` | Query TCM information |
| `get_tcm_chem_links(by, items)` | Query TCM-chemical links |
| `get_chemicals(by, items)` | Query chemical information |
| `get_chem_protein_links(by, items, score)` | Query chemical-protein links |
| `get_proteins(by, items)` | Query protein information |

### Analysis Functions

| Function | Description |
|----------|-------------|
| `from_tcm_or_formula(ids, ...)` | Forward analysis: TCM/formula to proteins |
| `from_proteins(proteins, ...)` | Reverse analysis: proteins to TCMs/formulas |
| `dfs_filter(...)` | DFS filter for complete pathways |
| `herb_analysis(proteins, ...)` | Integrated one-step analysis |

### Scoring Functions

| Function | Description |
|----------|-------------|
| `herbiv_score(tcm, ...)` | Calculate HerbiV scores |
| `component(items, ...)` | Optimize TCM/formula combinations |

### Network Functions

| Function | Description |
|----------|-------------|
| `filter_formula_network(...)` | Filter by formula composition |
| `filter_target_network(...)` | Filter by target proteins |
| `calculate_degree(network, node)` | Node degree centrality |
| `degree_summary(node_degree)` | Degree statistics by attribute |
| `plot_degree_boxplot(...)` | Degree distribution boxplot |
| `plot_degree_barchart(...)` | Key node bar chart |
| `network_analysis_pipeline(...)` | Complete post-analysis pipeline |

### Output Functions

| Function | Description |
|----------|-------------|
| `out_for_cyto(...)` | Export Cytoscape CSV files |
| `vis_network(...)` | Generate interactive D3.js HTML |

### GUI

| Function | Description |
|----------|-------------|
| `launch_app(port, host)` | Launch Streamlit web app |

## Data

This repository ships **source code only**. The relational tables, the
enrichment libraries, the syndrome profiles and the bundled CJK font are not
part of the source tree, and the MIT license below does **not** cover any
dataset.

| Dataset | Records |
|---------|---------|
| jnuherb_formula | ~6,200 |
| jnuherb_tcm | ~4,700 |
| jnuherb_chemicals | ~13,100 |
| jnuherb_proteins | ~14,200 |
| jnuherb_formula_tcm_links | ~28,300 |
| jnuherb_tcm_chemical_links | ~60,600 |
| jnuherb_chemical_protein_links | ~450,800 |

To run an analysis, place the tables in `src/jnuherb/data/` — or point
`JNUHERB_DATA_DIR` at any directory that holds them:

```bash
export JNUHERB_DATA_DIR=/path/to/jnuherb-data
```

Expected file names (see `_TABLE_FILES` in `src/jnuherb/_data.py`):

```
jnuherb_formula.csv              jnuherb_formula_tcm_links.csv
jnuherb_tcm.csv                  jnuherb_tcm_chemical_links.csv
jnuherb_chemicals.csv            jnuherb_chemical_protein_links.csv
jnuherb_proteins.csv             jnuherb_protein_disease_links.csv.gz
```

Optional extras: `data/enrichment/*.gmt` (+ `kegg_brite.tsv`) for enrichment
analysis, `data/fonts/NotoSansCJKsc-Regular.otf` for CJK output in exported
figures, and `data/*.json` for syndrome profiles.

Once the CSV tables are in place you can build the indexed SQLite backend
(faster queries, optional):

```bash
python -m jnuherb.build_db
```

## License

MIT — for the source code. Data is not distributed with this repository and
carries its own terms.
