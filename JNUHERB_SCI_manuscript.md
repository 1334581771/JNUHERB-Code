# JNUHERB: an interactive Python-based web server for network pharmacology analysis of traditional Chinese medicine with AI-assisted syndrome differentiation

**Authors** (template; please fill in)

First Author<sup>1,†</sup>, Second Author<sup>1,†</sup>, Third Author<sup>1</sup>, Corresponding Author<sup>1,*</sup>

<sup>1</sup> Jiangxi Normal University, College of XXX, Nanchang 330022, China

<sup>†</sup> These authors contributed equally to this work.

<sup>*</sup> To whom correspondence should be addressed. Email: xxx@jxnu.edu.cn

---

## Abstract

**Motivation:** Network pharmacology has become a pivotal approach for decoding the multi-component, multi-target mechanisms of traditional Chinese medicine (TCM). However, the majority of existing resources are implemented as R packages, command-line scripts, or static databases, creating a barrier for Python-oriented researchers and clinicians who require intuitive web interfaces. In addition, conventional tools rarely bridge the gap between TCM syndrome differentiation (a core concept in clinical TCM reasoning) and modern molecular target analysis.

**Results:** We present JNUHERB, an open-source, interactive web server that integrates a structured TCM knowledge graph, classical network pharmacology algorithms, and a large language model (LLM)-assisted syndrome differentiation module. JNUHERB contains 578,000+ curated records spanning 6,200 formulas, 4,700 herbs, 13,100 chemical ingredients, and 14,200 protein targets. The platform supports bidirectional network pharmacology analysis (from herbs/formulas/monomers to targets, and vice versa), HerbiV probability-union scoring, capacity-constrained combination optimization for TCM compatibility (配伍) design, interactive network visualization, and offline functional enrichment analysis of protein targets (KEGG, Gene Ontology, and Disease Ontology). A dedicated AI consultation module retrieves candidate herbs and formulas from the knowledge graph before invoking a large language model, and then translates natural-language symptom descriptions into standardized TCM syndromes and recommended herbs/formulas in both Chinese and English. All recommendations are cross-validated against the internal knowledge graph to mitigate hallucination. In addition, an integrated molecular docking module couples network-level scoring with structure-level validation: the highest-scoring monomers against a selected target can be docked in one click, or users can upload their own ligand and receptor files, to rapidly triage computational hits before experimental follow-up. JNUHERB is implemented in Python 3.9+ and provides 23 public APIs, a Streamlit-based graphical interface with thirteen functional panels, offline functional enrichment via bundled GMT gene-set files, and reproducible deployment via Docker with automatic startup.

**Availability:** JNUHERB is freely available at [URL]. The source code, documentation, example datasets, and Docker images are distributed under the MIT license at https://github.com/[username]/jnuherb.

---

## Keywords

traditional Chinese medicine; network pharmacology; web server; Python; artificial intelligence; retrieval-augmented generation; knowledge graph; herb-ingredient-target network; drug combination optimization; molecular docking; bilingual output

---

## 1. Introduction

Traditional Chinese medicine (TCM) has served as a major healthcare system for thousands of years and remains a valuable source of natural-product-based drug discovery (1,2). Unlike Western drugs that usually target single proteins with high specificity, TCM formulas are multi-component mixtures that exert therapeutic effects through synergistic modulation of multiple targets and pathways (3,4). This multi-scale complexity makes conventional reductionist approaches insufficient for mechanism elucidation, and has motivated the rapid development of network pharmacology, which models drug actions as perturbations on biological networks rather than isolated molecular events (5,6).

Over the past decade, numerous TCM-related databases and tools have been established, including TCMSP (7), BATMAN-TCM (8), TCMID (9), HERB (10), ETCM (11), and SymMap (12). These resources have greatly advanced the field by assembling herb-ingredient-target relationships and providing query interfaces. Nevertheless, several practical limitations remain. First, many tools are implemented as R packages (e.g., HerbiV (13)) or static web portals, which limits their integration with the broader Python-based data-science and machine-learning ecosystem. Second, existing platforms usually require users to specify herbs, ingredients, or targets explicitly, and rarely support clinical-style input such as natural-language symptom descriptions. Third, few tools connect TCM syndrome differentiation (a conceptual framework used by TCM practitioners) with molecular target analysis, leaving a gap between macroscopic clinical reasoning and microscopic mechanistic interpretation. Fourth, most platforms offer limited support for downstream optimization and visualization, forcing users to export data and switch between multiple software packages.

The emergence of large language models (LLMs) offers a new opportunity to narrow this gap. Modern LLMs can parse unstructured clinical descriptions and generate structured medical reasoning (14,15). However, when applied to specialized domains such as TCM, generic LLMs are prone to factual hallucinations, including the invention of non-existent herb or formula names. A promising strategy is therefore to combine LLM-based reasoning with structured domain knowledge graphs, an approach known as retrieval-augmented generation (RAG) (16). By anchoring LLM outputs to a curated knowledge base, such systems can retain the flexibility of natural-language interaction while ensuring factual verifiability.

To address these needs, we developed JNUHERB, a Python-based interactive web server for TCM network pharmacology with an integrated AI-assisted syndrome differentiation module and a molecular docking module for rapid hit triage. JNUHERB reconstructs the original R package HerbiV (13) in Python, extends it with a modern Streamlit graphical user interface, and adds an LLM-powered consultation layer that is cross-validated against the internal knowledge graph. The platform supports bidirectional analysis, multi-identifier queries, probability-union scoring, capacity-constrained combination optimization for TCM compatibility (配伍) design, molecular docking of network-ranked monomers against target proteins, and interactive visualization, enabling users to move seamlessly from symptom description to candidate formula, from disease target to molecular mechanism, and from computational ranking to structure-level validation in a single environment.

---

## 2. Materials and methods

### 2.1 Database curation and knowledge graph construction

JNUHERB is built around a four-layer knowledge graph: **Formula → Herb → Chemical → Protein**. The graph is derived from publicly available TCM and chemical-protein interaction databases, including HERB (10), ETCM (11), TCMSP (7), and STITCH (17). The curated dataset comprises seven relational tables:

- `jnuherb_formula`: 6,200+ TCM formulas;
- `jnuherb_tcm`: 4,700+ TCM herbs with multilingual identifiers (Chinese name, pinyin, English name, Latin name);
- `jnuherb_chemicals`: 13,100+ chemical ingredients;
- `jnuherb_proteins`: 14,200+ protein targets;
- `jnuherb_formula_tcm_links`: 28,300+ formula-herb associations;
- `jnuherb_tcm_chemical_links`: 60,600+ herb-chemical associations;
- `jnuherb_chem_protein_links`: 450,800+ chemical-protein interactions with combined scores.

All data were standardized by removing duplicate records, unifying aliases to a primary identifier, and converting identifiers to consistent string types. Chemical-protein interactions were retained with their original combined scores (0–1000), which users can threshold at runtime. The final packaged dataset totals approximately 578,000 records and is distributed as versioned CSV files inside the Python package, enabling offline use immediately after installation.

### 2.2 Bidirectional network pharmacology analysis

JNUHERB implements two complementary analysis modes.

**Forward analysis** (`from_tcm_or_formula`) accepts one or more herbs, formulas, or monomer chemicals and traverses the graph forward: herb/formula → chemical → protein target. This mode is useful for elucidating the molecular mechanism of a known formula or herb.

**Reverse analysis** (`from_proteins`) accepts one or more protein targets (Ensembl ID or gene symbol) and traverses the graph in reverse: protein → chemical → herb → formula. This mode supports target-driven drug discovery and repositioning.

To support flexible input, JNUHERB recognizes more than ten identifier types. Herbs can be queried by HVMID, Chinese name, pinyin, English name, or Latin name; formulas by HVPID or formula name (with fuzzy substring matching); monomers by compound name; and proteins by Ensembl ID or gene symbol.

### 2.3 HerbiV probability-union scoring and TCM combination (配伍) optimization

Following the HerbiV framework (13), JNUHERB quantifies the combined effect of multiple chemicals, herbs, or formulas on a protein target using probability-union scoring. For a set of independent interaction scores $s_i$ between entities and a target protein, the combined score is defined as:

$$
S = 1 - \prod_{i}(1 - s_i)
$$

The scoring is applied hierarchically: chemical scores are aggregated into herb scores, herb scores into formula scores, and formula scores into an overall importance score averaged across all query targets. Compared with simple summation, this formulation better captures the non-linear synergy arising from multi-component interventions.

For TCM compatibility (配伍) optimization, JNUHERB frames herb/formula selection as a capacity-constrained subset-selection problem. Each candidate item carries an importance value derived from its HerbiV score, and every item is assigned a uniform cost of one so that the combination size is bounded by a user-adjustable capacity parameter \(c\). The objective is the probability-union score \(1-\prod_{i}(1-v_i)\) rather than an additive sum, so the problem is a non-linear subset maximization rather than a classic 0/1 knapsack. To solve it, JNUHERB uses a dynamic-programming recurrence for the union objective, wrapped in repeated random restarts over random candidate subsamples, with a tabu list that suppresses already-reported combinations; this is a heuristic and does not guarantee a global optimum. The number of items in a reported combination is chosen automatically by an elbow (maximum-curvature) criterion on the resulting score-versus-size curve, reducing manual tuning. Each recommended combination is reported together with a *Boost* score, defined as the relative improvement of the combination over the best single item, so that users can distinguish genuinely synergistic multi-herb prescriptions from marginal additions. In the reverse-analysis workflow this module answers a question of direct clinical interest: given a set of disease targets, which combination of herbs—or which classical formulas—achieves the broadest mechanism coverage. We note that the union objective models target coverage under an independence assumption and does not, by itself, capture pharmacological drug–drug interactions.

### 2.4 Molecular docking for rapid hit triage

To bridge network-level ranking and structure-level evidence, JNUHERB integrates an automated molecular docking pipeline based on AutoDock Vina (22). The module is deliberately designed for *rapid triage* rather than high-throughput virtual screening: one ligand against one receptor completes in seconds to a few minutes on a standard dual-core server, allowing users to quickly confirm or deprioritize computational hits suggested by the network analysis.

The docking pipeline consists of four stages. (i) *Receptor preparation*: in the automatic mode, the Ensembl protein identifier of a target selected from the current analysis result is mapped to a UniProt accession through a pre-computed mapping table (13,773 of the 14,183 database proteins, 97.1%), and the corresponding AlphaFold-predicted structure (23) is downloaded from the AlphaFold DB and cached; in the upload mode, users provide their own receptor PDB file. The receptor is then protonated and converted to PDBQT with mk_prepare_receptor from Meeko. (ii) *Ligand preparation*: ligand SMILES taken from the knowledge graph (or an uploaded SDF/MOL/PDB file) are embedded in 3D with RDKit (ETKDGv3 followed by UFF relaxation) and written to PDBQT with Meeko. (iii) *Binding-site definition*: the docking box is centered on the centroid of the heavy atoms of user-specified binding-site residues; alternatively, users can enter an explicit box center, and a receptor C-alpha centroid fallback supports blind docking. (iv) *Docking and reporting*: Vina docking runs with exhaustiveness 8 and the best-scoring pose is kept. The receptor and the docked pose are merged into a single complex PDB for download, and an interactive 3D view (py3Dmol) displays the receptor—selectable as cartoon, cartoon-plus-surface, surface, spheres, sticks, or lines—together with the docked ligand.

Two entry points are provided. The *automatic mode* couples docking directly to the network pharmacology results: the target protein is chosen from the proteins of the current analysis, and candidate monomers are pre-ranked by their HerbiV score against that protein, so the top-scored ingredient can be docked in one click. The *upload mode* accepts arbitrary user-supplied ligands and receptors for off-database studies.

### 2.5 AI-assisted syndrome differentiation with knowledge-graph validation

The AI consultation module is designed to bridge TCM clinical reasoning with molecular analysis while avoiding hallucinated drug recommendations. The workflow consists of four layers:

1. **Retrieval layer.** Before invoking the LLM, JNUHERB extracts symptom keywords from the user's description and retrieves candidate herbs and formulas from the internal database. Herb candidates are scored by keyword overlap in fields such as `Function`, `Indication`, `Clinical_manifestations`, `Properties`, and `Meridians`; formula candidates are scored in `Treatment Symptoms`, `Prescription Composition`, and `Source Document`. Known Chinese symptom terms are expanded to their English synonyms so that the bilingual database fields are fully utilized. Only records present in the database are passed to the next layer.

2. **LLM reasoning layer.** The user-provided symptom description and the retrieved candidate list are sent to the DeepSeek Chat API with a domain-specific system prompt. The prompt instructs the model to act as an experienced TCM practitioner, return structured bilingual JSON (`syndrome`, `analysis`, `recommended_herbs`, `recommended_formulas`, `excription`), and **recommend only from the provided candidate list**. If no candidates are retrieved, the model is instructed to state that no matching database records were found rather than inventing drugs. The temperature is set to 0.3 to increase determinism and reduce hallucination.

3. **Knowledge-graph validation layer.** Each recommended herb and formula is cross-checked against `jnuherb_tcm` and `jnuherb_formula`. Herbs are matched by Chinese name (exact match); formulas are matched by substring fuzzy search to accommodate abbreviations (e.g., "Xiao Chai Hu" matches "Xiao Chai Hu Tang"). Validated items are annotated with their HVMID/HVPID and full attributes; unvalidated items are flagged as "not found in database."

4. **Network pharmacology analysis layer.** Validated herbs and formulas can be forwarded to the forward-analysis pipeline with a single click, producing ingredient-target networks, pathway enrichment-style target lists, and optimized combinations.

This retrieval + LLM + KG-RAG architecture ensures that AI-generated recommendations are both semantically coherent and factually anchored to the curated knowledge graph.

### 2.6 Web server implementation

JNUHERB is implemented entirely in Python 3.9+ and follows modern Python packaging standards. The package layout uses the `src/` structure and `pyproject.toml`. Core dependencies include pandas, NumPy, matplotlib, seaborn, plotly, and Streamlit. Data are accessed through `importlib.resources` and lazily loaded into memory with dictionary caching, providing millisecond-level response for typical queries after the initial ~2–3 second load.

The public API exposes 23 functions organized into five layers: query (7 functions), analysis (3 functions), scoring (2 functions), network analysis and visualization (7 functions), output (2 functions), pipeline (1 function), and GUI launching (1 function).

The Streamlit web interface contains thirteen panels: Overview, AI Consultation, TCM Results, Chemical Results, Protein Results, Network Visualization, Formula Filter, Degree Statistics, DO Enrichment, KEGG Enrichment, GO Enrichment, Optimization, and Docking. It supports interactive Plotly network views, Cytoscape-compatible CSV export, publication-ready static figures, and interactive 3D visualization of docking complexes via py3Dmol.

Functional enrichment is performed offline through gseapy against locally bundled GMT files (KEGG, GO BP/CC/MF, and Disease Ontology), so no external network access is required after installation. The docking module is the only component that requires external network access at runtime (downloading AlphaFold structures on first use); all cheminformatics dependencies (RDKit, Meeko, AutoDock Vina, ProDy, py3Dmol) are bundled in the Docker image, with Vina compiled in a multi-stage build so the runtime image remains free of compilers. The web interface also hides default Streamlit chrome (Deploy button, footer, header) to provide a clean, application-like experience.

Deployment options include local pip installation (`pip install -e .`), command-line launch (`jnuherb --port 8501`), Docker, Docker Compose, and Streamlit Cloud. A systemd service is provided so that the container starts automatically on server boot. The project is released under the MIT license.

---

## 3. Results and discussion

### 3.1 Overview of the JNUHERB framework

JNUHERB integrates data, algorithms, and user interaction into a unified analytical environment (Figure 1). The bottom layer is the curated knowledge graph, which supplies standardized entities and relationships. Above it, the algorithm layer performs query resolution, bidirectional network construction, HerbiV scoring, and combination optimization. The API layer exposes these capabilities as Python functions, and the Streamlit web layer renders them as interactive visual modules. Users can enter the system through either the Python API or the browser interface, making JNUHERB suitable for both computational workflows and bench-side exploration.

### 3.2 Database content and network characteristics

The knowledge graph contains 38,141 nodes and 539,699 relationships across the four layers (Table 1). Chemical-protein links dominate the graph (450,800 edges), reflecting the compound-centric nature of current molecular interaction data. Herbs are densely connected to chemicals (60,600 edges), and formulas to herbs (28,300 edges). Network analysis shows a highly skewed degree distribution: a small number of hub herbs (e.g., *Glycyrrhiza uralensis*, *Salvia miltiorrhiza*) and hub compounds appear in many formulas and target many proteins, consistent with the "core herb" phenomenon observed in TCM prescriptions (18).

### 3.3 Web server features

The Streamlit interface is organized into thirteen functional panels (Figure 2):

- **Overview** introduces analysis modes and summarizes results.
- **AI Consultation** accepts symptom descriptions, returns structured TCM reasoning, and validates recommendations against the knowledge graph.
- **TCM Results** displays herb-level results, importance-score distributions, and top-ranked herb bar charts.
- **Chemical Results** lists candidate chemical ingredients with their scores and distributions.
- **Protein Results** presents target proteins, degree statistics, and associated formula information.
- **Network Visualization** renders the top-100 core nodes as an interactive Plotly graph and provides Cytoscape export guidance.
- **Formula Filter** enables subnetwork filtering by formula composition in reverse-analysis mode.
- **Degree Statistics** provides degree-distribution boxplots and top-node bar charts.
- **DO Enrichment**, **KEGG Enrichment**, and **GO Enrichment** perform offline over-representation analysis on protein targets using bundled gene-set files, presenting ranked terms together with bar plots, dot plots, and downloadable CSV files.
- **Optimization** reports combinations produced by the capacity-constrained union optimizer for TCM compatibility (配伍) design. Each combination is displayed with its component names in both Chinese and English (rather than raw database identifiers), together with the combination score and a Boost value quantifying the synergistic gain over the best single herb or formula.
- **Docking** performs one-click molecular docking between a target protein selected from the current analysis and its highest-scoring monomer (or user-uploaded structures), reports the Vina score, and renders an interactive 3D view of the receptor–ligand complex with a downloadable complex PDB file.

All panels are linked by a shared analysis state, allowing users to move from AI consultation to forward analysis to network visualization, and from network ranking to docking validation, without re-entering data.

### 3.4 Case study: from symptoms to molecular mechanism

To illustrate the end-to-end utility of JNUHERB, we analyzed the symptom cluster "insomnia, dream-disturbed sleep, palpitations, forgetfulness, fatigue, pale tongue with thin coating." In TCM theory, this pattern is commonly interpreted as "heart-and-spleen deficiency" (心脾两虚证).

Step 1: AI consultation. The user enters the symptoms in natural language. The LLM returns the syndrome "heart-and-spleen deficiency" (心脾两虚证) in both Chinese and English, along with recommended herbs (e.g., *Ren Shen* 人参, *Huang Qi* 黄芪, *Suan Zao Ren* 酸枣仁, *Long Yan Rou* 龙眼肉) and formulas (e.g., *Gui Pi Tang* 归脾汤, *Suan Zao Ren Tang* 酸枣仁汤). JNUHERB validates these against the knowledge graph: all recommended herbs and formulas are found in the database with their HVMID/HVPID.

Step 2: Forward molecular analysis. The validated formula *Gui Pi Tang* is forwarded to the forward-analysis pipeline. JNUHERB retrieves 156 candidate protein targets and the chemical ingredients responsible for the connections.

Step 3: Scoring and compatibility (配伍) optimization. HerbiV scoring ranks herbs and compounds by their multi-target importance. The union-based optimizer recommends a core subset of four herbs that maximizes target coverage while minimizing redundancy; each recommended combination is reported with a Boost value indicating its gain over the best single herb, supporting data-driven TCM compatibility design.

Step 4: Functional enrichment. The 156 protein targets are automatically forwarded to the enrichment panels, which report top KEGG pathways, GO biological process/cellular component/molecular function terms, and Disease Ontology terms using bundled human gene-set files. Results are visualized as horizontal bar charts and dot plots and can be exported as CSV files.

Step 5: Visualization and docking validation. The resulting herb-chemical-target network is rendered as an interactive Plotly graph, and a Cytoscape-compatible file pair is exported for further customization (Figure 3). Finally, the highest-scoring monomer against a selected target is forwarded to the docking panel with one click: the receptor structure is fetched automatically from the AlphaFold DB, the docking box is defined from known binding-site residues, and Vina docking returns a binding score and an interactive 3D view of the complex within minutes—providing structure-level support for the network-derived ranking before any experimental work.

This case demonstrates how JNUHERB connects macroscopic TCM syndrome differentiation with microscopic target analysis and structure-level validation in a reproducible workflow.

### 3.5 Case study: from disease targets to candidate formulas

In the reverse direction, we queried JNUHERB with three inflammation- and tumor-related targets: TNF, IL6, and TP53. JNUHERB returned 47 candidate herbs (including *Huang Qin*, *Huang Lian*, *Dan Shen*) and 23 candidate formulas (including *Huang Lian Jie Du Tang*, *Qing Ying Tang*). The top-ranked compounds included baicalein, berberine, and tanshinone IIA, all of which have reported anti-inflammatory or anti-tumor activities (19–21). The compatibility (配伍) optimizer suggested a four-herb combination with broad target coverage and a clear Boost over the best single herb, providing a starting point for experimental follow-up. Because the module reports combinations for both herbs and classical formulas side by side, users can directly compare a custom multi-herb prescription against established formulas that cover the same targets. The top-ranked monomers can then be docked against any selected target in the Docking panel for rapid structure-level triage of the computational hits.

### 3.6 Comparison with existing tools

Table 2 compares JNUHERB with representative TCM network pharmacology platforms. Unlike TCMSP, ETCM, and HERB, which are primarily data repositories, JNUHERB provides an integrated analysis pipeline. Unlike the original HerbiV R package, JNUHERB is written in Python, offers a modern web interface, and adds AI-assisted syndrome differentiation, knowledge-graph validation, and molecular docking. To our knowledge, JNUHERB is the first open-source Python platform that combines LLM-based clinical reasoning, structured TCM knowledge-graph validation, classical network pharmacology, and automated molecular docking in a single tool.

---

## 4. Discussion

JNUHERB addresses a growing need in TCM informatics: a Python-native, web-accessible platform that supports both molecular network pharmacology and clinical-style TCM reasoning. By reconstructing the HerbiV framework in Python and adding a Streamlit interface, we lower the barrier for researchers who prefer Python or who need browser-based tools. By integrating a retrieval-augmented, LLM-assisted syndrome differentiation module with knowledge-graph validation, we introduce a new mode of interaction—natural-language symptom input—while ensuring that every drug recommendation is first retrieved from and then validated against the curated database, thereby mitigating the hallucination risks inherent in generic LLMs.

The design choices reflect practical constraints. We use CSV files packaged with the software rather than a relational database, because this maximizes portability and offline usability: users can install JNUHERB behind institutional firewalls without configuring a database server. Lazy loading and dictionary caching keep typical queries fast despite the ~578,000-record dataset. The choice of Streamlit balances development speed and interactivity, making the platform easy to extend and deploy.

Several limitations should be acknowledged. First, the chemical-protein interaction scores are derived from existing databases and represent computational predictions rather than experimentally validated binding affinities; the docking module provides structure-level triage but relies on AlphaFold-predicted receptor structures and single-pose Vina scoring, so docking scores should likewise be interpreted as prioritization aids rather than definitive evidence of bioactivity. Second, the AI consultation module relies on an external LLM API; users without API access can still use all network pharmacology and docking functions but cannot access the AI consultation feature. Third, the retrieval layer depends on the completeness and accuracy of TCM symptom/function annotations in the database; sparse or inconsistent annotations may affect candidate recall. Fourth, the current knowledge graph focuses on herb-ingredient-target relationships and does not yet include pathways, diseases, or pharmacokinetic properties; these will be added in future releases.

Future development will focus on four directions: (i) expanding the knowledge graph with pathway and disease annotations from KEGG and Reactome; (ii) integrating ADMET prediction to support lead-compound optimization; (iii) developing a TCM-specific fine-tuned language model to reduce dependence on generic LLMs; and (iv) extending the docking module toward ensemble docking with induced-fit refinement and user-defined pharmacophore constraints.

---

## 5. Conclusion

JNUHERB is an open-source, interactive Python web server for TCM network pharmacology that combines a curated knowledge graph, bidirectional network analysis, probability-union scoring, TCM compatibility (配伍) optimization, automated molecular docking, and AI-assisted syndrome differentiation. By bridging clinical TCM reasoning with molecular target analysis and structure-level validation, JNUHERB provides a unified environment for hypothesis generation, mechanism exploration, compatibility design, and candidate formula prioritization. We anticipate that JNUHERB will serve as a practical tool for TCM researchers, natural-product chemists, and clinicians interested in data-driven herbal medicine research.

---

## Data availability

The JNUHERB web server is freely available at [URL]. The source code, documentation, example datasets, Docker images, and deployment guides are hosted at https://github.com/[username]/jnuherb under the MIT license. The built-in database is distributed as seven CSV files within the Python package and can be used offline after installation. Offline enrichment analysis requires the bundled GMT files (KEGG, GO BP/CC/MF, and Disease Ontology), which are also shipped inside the package.

## Implementation

JNUHERB is implemented in pure Python (3.9+) using `pandas`, `NumPy`, `matplotlib`, `seaborn`, `plotly`, `gseapy`, and `Streamlit`. The docking module additionally uses RDKit, Meeko, ProDy, AutoDock Vina, and py3Dmol. The package uses a standard `src/` layout and is configured via `pyproject.toml`. Container deployment is supported through a multi-stage Dockerfile and `docker-compose.yml`; a systemd service template enables automatic startup after server reboot. Chinese users can optionally enable Aliyun apt/pip mirrors and Docker Hub mirror acceleration by following the provided deployment guide.

---

## Acknowledgements

We thank the developers of the HerbiV R package for the original scoring framework. We also acknowledge the curators of TCMSP, ETCM, HERB, and STITCH for making their data publicly available.

---

## Author contributions

F.A. and S.A. contributed equally to software development, data curation, and manuscript writing. T.A. participated in data collection and testing. C.A. conceived the study, supervised the project, and revised the manuscript.

---

## References

1. Newman DJ, Cragg GM. Natural products as sources of new drugs over the nearly four decades from 01/1981 to 09/2019. J Nat Prod. 2020;83(3):770–803.
2. Li S, Zhang B. Traditional Chinese medicine network pharmacology: theory, methodology and application. Chin J Nat Med. 2013;11(2):110–120.
3. Hopkins AL. Network pharmacology: the next paradigm in drug discovery. Nat Chem Biol. 2008;4(11):682–690.
4. Borisy AA, Elliott PJ, Hurst NW, et al. Systematic discovery of multicomponent therapeutics. Proc Natl Acad Sci USA. 2003;100(13):7977–7982.
5. Kibble M, Saarinen N, Tang J, et al. Network pharmacology applications to map the unexplored target space and therapeutic potential of natural products. Nat Prod Rep. 2015;32(8):1249–1266.
6. Ru J, Li P, Wang J, et al. TCMSP: a database of systems pharmacology for drug discovery from herbal medicines. J Cheminform. 2014;6:13.
7. Ru J, Li P, Wang J, et al. TCMSP: a database of systems pharmacology for drug discovery from herbal medicines. J Cheminform. 2014;6:13.
8. Liu Z, Guo F, Wang Y, et al. BATMAN-TCM: a Bioinformatics Analysis Tool for Molecular mechANism of Traditional Chinese Medicine. Sci Rep. 2016;6:21146.
9. Xue R, Fang Z, Zhang M, et al. TCMID: Traditional Chinese Medicine Integrative Database for Herb Molecular Mechanism Analysis. Nucleic Acids Res. 2013;41(D1):D1089–D1095.
10. Fang S, Dong L, Liu L, et al. HERB: a high-throughput experiment- and reference-guided database of traditional Chinese medicine. Nucleic Acids Res. 2021;49(D1):D1197–D1206.
11. Xu HY, Zhang YQ, Liu ZM, et al. ETCM: an encyclopaedia of traditional Chinese medicine. Nucleic Acids Res. 2019;47(D1):D976–D982.
12. Wu Y, Zhang F, Yang K, et al. SymMap: an integrative database of traditional Chinese medicine enhanced by symptom mapping. Nucleic Acids Res. 2019;47(D1):D1110–D1117.
13. Fang S, Liu L, Ji H, et al. HerbiV: a computational platform for network pharmacology analysis of traditional Chinese medicines. (Original R package reference; please update with correct citation.)
14. Singhal K, Azizi S, Tu T, et al. Large language models encode clinical knowledge. Nature. 2023;620(7972):172–180.
15. Moor M, Banerjee O, Abad ZSH, et al. Foundation models for generalist medical artificial intelligence. Nature. 2023;616(7956):259–265.
16. Lewis P, Perez E, Piktus A, et al. Retrieval-augmented generation for knowledge-intensive NLP tasks. Adv Neural Inf Process Syst. 2020;33:9459–9474.
17. Szklarczyk D, Gable AL, Nastou KC, et al. The STRING database in 2021: customizable protein-protein networks, and functional characterization of user-uploaded gene/measurement sets. Nucleic Acids Res. 2021;49(D1):D605–D612.
18. Li B, Cao H, Gao Y, et al. A network pharmacology approach to determine the active ingredients and potential targets of Tripterygium wilfordii Hook F for the treatment of rheumatoid arthritis. Evid Based Complement Alternat Med. 2019;2019:5249817.
19. Li Y, Wang X, Li X, et al. Baicalein suppresses 1,25D3-MARRS-induced inflammation in macrophages. Am J Chin Med. 2019;47(5):1103–1121.
20. Tan W, Li Y, Chen M, et al. Berberine exhibits antitumor effects in human lung cancer cells via the inactivation of PI3K/AKT/mTOR pathway. J Cancer Res Clin Oncol. 2020;146(7):1705–1719.
21. Fang ZY, Zhang YH, Qin XC, et al. Tanshinone IIA protects against diabetic nephropathy via inhibiting TGF-β1/Smad pathway and improving mitochondrial function. Front Pharmacol. 2021;12:637608.
22. Eberhardt J, Santos-Martins D, Tillack AF, Forli S. AutoDock Vina 1.2.0: new docking methods, expanded force field, and Python bindings. J Chem Inf Model. 2021;61(8):3891–3898.
23. Jumper J, Evans R, Pritzel A, et al. Highly accurate protein structure prediction with AlphaFold. Nature. 2021;596(7873):583–589.

---

## Figure captions

**Figure 1.** Overview of the JNUHERB architecture. The platform is organized into four layers: (A) the curated knowledge graph containing formulas, herbs, chemicals, and proteins; (B) the algorithm layer implementing query resolution, bidirectional network construction, HerbiV scoring, TCM compatibility (配伍) optimization, and molecular docking; (C) the Python API layer exposing 23 public functions; and (D) the Streamlit web layer providing thirteen interactive analysis panels.

**Figure 2.** The JNUHERB Streamlit web interface. (A) Overview panel for selecting analysis modes and viewing result summaries. (B) AI Consultation panel showing natural-language symptom input, retrieval-grounded bilingual (Chinese/English) syndrome differentiation, database candidates used for the recommendation, and knowledge-graph validation of recommended herbs and formulas. (C) TCM Results panel with interactive tables and score distributions. (D) Protein Results panel. (E) Network Visualization panel rendering the top-100 core nodes as a Plotly graph. (F) Formula Filter and Degree Statistics panels. (G) DO, KEGG, and GO Enrichment panels showing ranked terms, bar charts, dot plots, and CSV download options. (H) Optimization panel displaying herb and formula combinations recommended by the union-based optimizer, with bilingual component names and Boost scores. (I) Docking panel showing target/ligand selection from analysis results, the reported Vina score, an interactive 3D view of the receptor–ligand complex, and complex PDB download.

**Figure 3.** Case study of the symptom pattern "heart-and-spleen deficiency." (A) AI-generated syndrome and validated formula recommendations. (B) Forward analysis of *Gui Pi Tang* showing herb-chemical (orange), chemical-protein (green), and formula-herb (blue) relationships. (C) Top-20 target proteins ranked by HerbiV importance score. (D) Optimized four-herb combination recommended by the union-based optimizer.

**Figure 4.** Molecular docking workflow of JNUHERB. (A) Automatic mode: a target protein selected from the analysis results is mapped from Ensembl to UniProt and its AlphaFold structure is retrieved; candidate monomers are pre-ranked by HerbiV score. (B) Ligand and receptor preparation (RDKit/Meeko; mk_prepare_receptor). (C) Docking-box definition from binding-site residue numbers, a manual center, or the receptor C-alpha centroid fallback. (D) AutoDock Vina docking and reporting: Vina score, interactive 3D view with selectable receptor representations (cartoon, surface, spheres, sticks, lines), and downloadable receptor–ligand complex PDB.

---

## Tables

**Table 1.** Composition of the JNUHERB knowledge graph.

| Dataset | Entity type | Records |
|---|---|---|
| jnuherb_formula | Formula | ~6,200 |
| jnuherb_tcm | Herb | ~4,700 |
| jnuherb_chemicals | Chemical | ~13,100 |
| jnuherb_proteins | Protein | ~14,200 |
| jnuherb_formula_tcm_links | Formula–Herb | ~28,300 |
| jnuherb_tcm_chemical_links | Herb–Chemical | ~60,600 |
| jnuherb_chem_protein_links | Chemical–Protein | ~450,800 |

**Table 2.** Comparison of JNUHERB with representative TCM network pharmacology tools.

| Feature | TCMSP | ETCM | HERB | SymMap | HerbiV (R) | JNUHERB |
|---|---|---|---|---|---|---|
| Web interface | Yes | Yes | Yes | Yes | No | Yes |
| Python API | No | No | No | No | No | Yes |
| Open source | Partial | No | No | No | Yes | Yes |
| Bidirectional analysis | Limited | Limited | Limited | Limited | Yes | Yes |
| Network visualization | Basic | Basic | Basic | Basic | Yes | Yes |
| AI syndrome differentiation | No | No | No | No | No | Yes |
| Bilingual (Chinese/English) output | No | No | No | No | No | Yes |
| Knowledge-graph validation | No | No | No | No | No | Yes |
| Retrieval-grounded AI recommendation | No | No | No | No | No | Yes |
| TCM compatibility (配伍) optimization | No | No | No | No | Yes | Yes |
| Molecular docking | No | No | No | No | No | Yes |
| Functional enrichment (KEGG/GO/DO) | No | No | No | No | No | Yes |
| Offline deployment | No | No | No | No | Yes | Yes |

---

*Note to authors: Before submission, replace bracketed placeholders ([URL], [username]) with actual values, verify all author affiliations and contributions, and update reference 13 with the correct HerbiV citation.*
