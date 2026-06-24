#!/usr/bin/env python
"""JNUHERB — Interactive Network Pharmacology Analysis (Streamlit GUI).

Supports two analysis modes:

* **Reverse**: proteins → chemicals → TCMs → formulas
* **Forward**: TCMs / formulas / monomers → chemicals → proteins

Launch via::

    python -m jnuherb
    jnuherb
    streamlit run src/jnuherb/data/streamlit_app.py
"""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Optional

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st

matplotlib.use("Agg")

from jnuherb.analysis import from_proteins, from_tcm_or_formula
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
st.set_page_config(
    page_title="JNUHERB",
    page_icon="🌿",
    layout="wide",
    initial_sidebar_state="expanded",
)

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


# ---------------------------------------------------------------------------
# Sidebar — Mode & Inputs
# ---------------------------------------------------------------------------
st.sidebar.title("🌿 JNUHERB")
st.sidebar.markdown("---")

# ---- Analysis Mode ----
st.sidebar.subheader("Analysis Mode")
analysis_mode = st.sidebar.radio(
    "Choose direction",
    ["Reverse (Proteins → TCMs)", "Forward (TCMs / Formulas → Proteins)"],
    key="analysis_mode",
)
is_forward = analysis_mode.startswith("Forward")

st.sidebar.markdown("---")

# ---- Inputs (conditional on mode) ----
if is_forward:
    # ===== FORWARD MODE =====
    st.sidebar.subheader("Query Entity Type")
    forward_entity_type = st.sidebar.selectbox(
        "What are you searching for?",
        ["TCM (Herb)", "Formula (方剂)", "Monomer (单体)"],
        key="forward_entity_type",
    )

    if forward_entity_type == "TCM (Herb)":
        forward_text = st.sidebar.text_area(
            "TCM IDs or names (one per line)",
            placeholder="HVM0367\nHVM1695\n柴胡\n黄芩",
            height=110,
            key="forward_text_tcm",
        )
        forward_id_type = st.sidebar.selectbox(
            "ID Type",
            ["HVMID", "cn_name", "pinyin_name", "en_name", "latin_name"],
            help="HVMID=herb ID, cn_name=Chinese name, pinyin_name=Pinyin, "
            "en_name=English name, latin_name=Latin name",
            key="forward_id_type_tcm",
        )
    elif forward_entity_type == "Formula (方剂)":
        forward_text = st.sidebar.text_area(
            "Formula IDs or names (one per line)",
            placeholder="HVP1625\n右归丸\n小柴胡汤",
            height=110,
            key="forward_text_fml",
        )
        forward_id_type = st.sidebar.selectbox(
            "ID Type",
            ["HVPID", "Formula Name"],
            help="HVPID=formula ID (e.g. HVP1625), Formula Name=search by name",
            key="forward_id_type_fml",
        )
    else:  # Monomer
        forward_text = st.sidebar.text_area(
            "Chemical / Monomer names (one per line)",
            placeholder="deoxyshikonin\nberberine\nquercetin",
            height=110,
            key="forward_text_mono",
        )
        forward_id_type = "Name"

    st.sidebar.subheader("Optional Filters")
    forward_proteins_text = st.sidebar.text_area(
        "Restrict to specific proteins (optional, one per line)",
        placeholder="ENSP00000381588\nENSP00000252519",
        height=80,
        key="forward_proteins",
    )
    forward_score = st.sidebar.number_input(
        "Combined Score Threshold",
        value=900 if forward_entity_type != "Monomer (单体)" else 0,
        min_value=0, max_value=1000,
        help="Lower = more results, Higher = more confident links",
        key="forward_score",
    )
else:
    # ===== REVERSE MODE =====
    st.sidebar.subheader("Target Proteins")
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

    st.sidebar.subheader("Formula Filter (optional)")
    formula_text = st.sidebar.text_area(
        "TCM names to filter by (one per line, optional)",
        placeholder="熟地黄\n山药\n山茱萸",
        height=80,
        key="formula_text",
    )
    formula_id_type = st.sidebar.selectbox(
        "Formula ID Type",
        ["cn_name", "HVMID", "pinyin_name", "en_name", "latin_name"],
        key="formula_id_type",
    )

    reverse_score = st.sidebar.number_input(
        "Combined Score Threshold",
        value=0, min_value=0, max_value=1000,
        key="reverse_score",
    )

# ---- Common parameters ----
st.sidebar.markdown("---")
st.sidebar.subheader("Optimization Parameters")
opt_num = st.sidebar.number_input(
    "Optimization Iterations", value=100, min_value=10, max_value=5000,
)
opt_seed = st.sidebar.number_input(
    "Random Seed", value=138192, min_value=0,
)
tcm_opt_enabled = st.sidebar.checkbox(
    "TCM Combination Optimization", value=True,
)
fml_opt_enabled = st.sidebar.checkbox(
    "Formula Combination Optimization", value=True,
)

st.sidebar.markdown("---")
run_btn = st.sidebar.button(
    "🚀 Run Analysis", type="primary", use_container_width=True,
)


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
        st.error(f"Reverse analysis failed: {e}")
        return None

    return _build_network_results(result, formula_list, output_dir, is_forward=False)


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
    - "Name" → monomer (chemical name) query
    - "Formula Name" → search formula by name, then run forward analysis
    - everything else → standard from_tcm_or_formula
    """
    entities_list = list(entities_tuple)
    proteins_list = list(proteins_tuple) if proteins_tuple else None

    # --- Monomer: chemical name query ---
    if id_type == "Name":
        return _run_monomer_forward(
            chem_names=entities_list,
            proteins_list=proteins_list,
            score=score,
            output_dir=output_dir,
        )

    # --- Formula Name: resolve name → HVPID first ---
    if id_type == "Formula Name":
        from jnuherb._data import _load_table

        fml_all = _load_table("jnuherb_formula")
        all_hvpids = []
        for name in entities_list:
            mask = fml_all["name"].str.contains(name, na=False, case=False)
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
        st.error(f"Forward analysis failed: {e}")
        return None

    return _build_network_results(result, None, output_dir, is_forward=True)


def _run_monomer_forward(
    chem_names: list[str],
    proteins_list: list[str] | None,
    score: float,
    output_dir: str,
) -> dict | None:
    """Monomer-first forward analysis: chemical name → targets → TCMs → formulas."""
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

    # Step 1: Find chemicals by name
    chem = get_chemicals("Name", chem_names)
    if len(chem) == 0:
        st.error(f"No chemicals found for: {', '.join(chem_names)}")
        return None

    # Step 2: Get TCMs containing these chemicals
    tcm_chem_links = get_tcm_chem_links("HVCID", list(chem["HVCID"]))
    if len(tcm_chem_links) == 0:
        st.error("No TCMs found containing these chemicals.")
        return None

    tcm = get_tcm("HVMID", list(tcm_chem_links["HVMID"]))

    # Step 3: Get chemical-protein links
    from jnuherb.query import get_chem_protein_links, get_proteins
    chem_prot_links = get_chem_protein_links("HVCID", list(chem["HVCID"]), score)
    if len(chem_prot_links) == 0:
        st.error(f"No protein targets found at score >= {score}. Try lowering the threshold.")
        return None

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
        # Build from raw result
        tcl = api_result.get("tcm_chem_links")
        cpl = api_result.get("chem_protein_links")
        if tcl is not None and cpl is not None:
            network = pd.DataFrame({
                "SourceNode": list(tcl.iloc[:, 0]) + list(cpl.iloc[:, 0]),
                "TargetNode": list(tcl.iloc[:, 1]) + list(cpl.iloc[:, 1]),
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
# Run on button click
# ---------------------------------------------------------------------------
if run_btn:
    output_dir = Path(tempfile.gettempdir()) / "jnuherb_streamlit"
    output_dir.mkdir(parents=True, exist_ok=True)

    if is_forward:
        entities = parse_lines(forward_text)
        if not entities:
            st.error("Please enter at least one TCM/Formula/Monomer ID or name!")
        else:
            # --- Formula Name: resolve before cache ---
            entities_for_cache = entities
            id_type_for_cache = forward_id_type

            if forward_id_type == "Formula Name":
                from jnuherb._data import _load_table

                fml_df = _load_table("jnuherb_formula")
                all_hvpids = []
                found, not_found = [], []
                for name in entities:
                    # substring match on formula name
                    mask = fml_df["name"].str.contains(name, na=False, case=False)
                    matches = fml_df[mask]
                    if len(matches) > 0:
                        all_hvpids.extend(matches["HVPID"].tolist())
                        found.append(name)
                    else:
                        not_found.append(name)

                if not_found:
                    st.warning(
                        f"Formula name(s) not found: {', '.join(not_found)}. "
                        f"Try using HVPID directly or a different name."
                    )
                if not all_hvpids:
                    st.error("No matching formulas found. Check the names or try HVPID.")
                    st.stop()

                entities_for_cache = list(set(all_hvpids))
                id_type_for_cache = "HVMID"  # HVP prefix auto-detects as formula
                st.info(
                    f"Resolved to {len(entities_for_cache)} formula(s): "
                    f"{', '.join(entities_for_cache[:10])}"
                    f"{'...' if len(entities_for_cache) > 10 else ''}"
                )

            prot_filter = parse_lines(forward_proteins_text)
            st.session_state.result = _run_forward(
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
    else:
        proteins = parse_lines(reverse_text)
        if not proteins:
            st.error("Please enter at least one target protein ID!")
        else:
            fml_list = parse_lines(formula_text)
            st.session_state.result = _run_reverse(
                proteins_tuple=tuple(proteins),
                formula_tuple=tuple(fml_list) if fml_list else None,
                protein_id_type=reverse_id_type,
                formula_id_type=formula_id_type,
                score=float(reverse_score),
                num=int(opt_num),
                random_state=int(opt_seed),
                tcm_comp=tcm_opt_enabled,
                fml_comp=fml_opt_enabled,
                output_dir=str(output_dir),
            )

# ---------------------------------------------------------------------------
# Tabs
# ---------------------------------------------------------------------------
tabs = st.tabs([
    "📋 Overview",
    "💬 AI Consultation",
    "🌿 TCM Results",
    "⚗️ Chemical Results",
    "🧬 Protein Results",
    "🔗 Network Viz",
    "🔍 Formula Filter",
    "📊 Degree Stats",
    "⚙️ Optimization",
])

result = st.session_state.get("result")
is_fwd = result.get("is_forward", False) if result else False

# ---- Tab 1: Overview ----
with tabs[0]:
    if result is None:
        st.info("👈 Choose an analysis mode, enter query terms, and click **Run Analysis**.")
        st.markdown("""
        ### Quick Start Guide

        **Reverse mode** (Proteins → TCMs):
        1. Select *Reverse* mode in the sidebar.
        2. Enter target protein IDs (Ensembl IDs or gene names).
        3. Optionally filter by specific TCM formula.
        4. Click **Run Analysis**.

        **Forward mode** (TCMs / Formulas / Monomers → Proteins):
        1. Select *Forward* mode and choose an entity type.
        2. **TCM**: enter herb IDs or names (e.g. 柴胡, HVM0367).
        3. **Formula**: enter formula IDs or names (e.g. HVP1625, 右归丸).
        4. **Monomer**: enter chemical names (e.g. deoxyshikonin).
        5. Optionally restrict to specific proteins and set a score threshold.
        6. Click **Run Analysis**.

        ---
        **Example — Reverse (Ensembl IDs):**
        `ENSP00000381588` (TMPRSS2), `ENSP00000252519` (ALDH2)

        **Example — Reverse (gene names):**
        `TNF, IL6, TP53, PDCD1`

        **Example — Forward TCM:**
        `HVM0367, HVM1695` or `柴胡, 黄芩`

        **Example — Forward Formula:**
        `HVP1625` or `右归丸`

        **Example — Forward Monomer:**
        `deoxyshikonin, berberine`

        > ℹ️ The protein database contains ~14,200 proteins with known
        > TCM-related chemical interactions, not the full human proteome.
        > If your Ensembl IDs are not found, try using gene names instead.
        """)
    else:
        pr = result["protein_result"]
        ff = result["formula_filter"]

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

        st.markdown("---")
        st.subheader("Analysis Log")
        log_lines = [
            f"✅ Analysis completed!",
            f"Mode: {'Forward' if is_fwd else 'Reverse'}",
            f"TCMs: {len(pr.get('tcm', []))}",
            f"Chemicals: {len(pr.get('chem', []))}",
            f"Proteins: {len(pr.get('proteins', []))}",
            f"Network edges: {len(result['network'])}",
        ]
        if ff:
            log_lines.append(f"Formula TCMs found: {len(ff.get('tcms_found', []))}")
        st.code("\n".join(log_lines), language=None)

# ---- Tab 1: AI Consultation ----
with tabs[1]:
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

        **Examples:**
        - *头痛发热，恶寒无汗，脉浮紧*
        - *失眠多梦，心悸健忘，神疲乏力，舌淡苔薄*
        - *胃脘胀痛，嗳气频作，食欲不振，大便不畅*
        """)

        col1, col2 = st.columns([3, 1])
        with col1:
            symptoms_input = st.text_area(
                "Symptoms / Disease Description",
                placeholder="请描述症状...",
                height=100,
                key="ai_symptoms",
            )
        with col2:
            st.caption("💡 **Tips**")
            st.caption("- Include key signs")
            st.caption("- Mention tongue/pulse if known")
            st.caption("- More detail = better analysis")

        if st.button("🔍 Analyze Symptoms", type="primary", key="ai_analyze_btn"):
            if not symptoms_input.strip():
                st.warning("Please enter symptom description.")
            else:
                try:
                    from jnuherb.chat import consult, check_database

                    with st.spinner("Analyzing with DeepSeek AI..."):
                        ai_result = consult(symptoms_input.strip())

                    st.session_state.consultation = ai_result

                    # Cross-reference with database
                    with st.spinner("Cross-referencing with database..."):
                        db_check = check_database(
                            ai_result.get("herbs", []),
                            ai_result.get("formulas", []),
                        )
                    st.session_state.consultation["db_check"] = db_check
                    st.rerun()

                except ValueError as e:
                    st.error(
                        f"⚠️ {e}\n\n"
                        "Please set your DeepSeek API key below."
                    )
                except RuntimeError as e:
                    st.error(f"❌ {e}")

        # API key section (collapsible)
        with st.expander("⚙️ API Settings"):
            api_key_input = st.text_input(
                "DeepSeek API Key",
                type="password",
                placeholder="sk-...",
                help="Get your key at https://platform.deepseek.com",
                key="api_key_input",
            )
            if api_key_input:
                import os

                os.environ["DEEPSEEK_API_KEY"] = api_key_input.strip()
                st.success("✅ API key set for this session.")

            st.caption(
                "Your API key is only stored in memory for this session. "
                "To persist it, set the `DEEPSEEK_API_KEY` environment variable "
                "or add it to `.streamlit/secrets.toml`."
            )
    else:
        # Show consultation results
        ai = st.session_state.consultation

        c1, c2 = st.columns([2, 1])
        with c1:
            st.markdown(f"### 🏥 {ai.get('syndrome', 'TCM Analysis')}")
            st.markdown(f"**Analysis:** {ai.get('analysis', '')}")
            st.markdown(f"**Treatment Rationale:** {ai.get('explanation', '')}")
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
            st.subheader("🌿 Recommended Herbs")
            if herbs_found:
                herb_df_data = []
                for h in herbs_found:
                    herb_df_data.append({
                        "Herb": h["name"],
                        "HVMID": h["HVMID"],
                        "Pinyin": h.get("pinyin", ""),
                        "English": h.get("en_name", ""),
                    })
                st.dataframe(herb_df_data, use_container_width=True)

                # Button to analyze these herbs
                herb_names = [h["name"] for h in herbs_found]
                if st.button("🔬 Analyze These Herbs", key="ai_analyze_herbs"):
                    st.session_state.forward_preset = herb_names
                    st.session_state.forward_preset_type = "cn_name"
                    st.success(
                        f"Herb names copied! Switch to **Forward mode** "
                        f"→ **TCM** and paste these names to run analysis."
                    )
            if herbs_missing:
                st.caption(f"⚠️ Not in database: {', '.join(herbs_missing)}")

        with col_f:
            st.subheader("📜 Recommended Formulas")
            if formulas_found:
                fml_df_data = []
                for f in formulas_found:
                    fml_df_data.append({
                        "Formula": f["full_name"],
                        "HVPID": f["HVPID"],
                        "Source": (
                            f["source"][:60] + "..."
                            if len(f.get("source", "")) > 60
                            else f.get("source", "")
                        ),
                    })
                st.dataframe(fml_df_data, use_container_width=True)

                fml_hvpids = [f["HVPID"] for f in formulas_found]
                if st.button("🔬 Analyze These Formulas", key="ai_analyze_formulas"):
                    st.session_state.forward_preset = fml_hvpids
                    st.session_state.forward_preset_type = "HVMID"
                    st.success(
                        f"Formula IDs copied! Switch to **Forward mode** "
                        f"→ **Formula** and use HVPID to run analysis."
                    )
            if formulas_missing:
                st.caption(f"⚠️ Not in database: {', '.join(formulas_missing)}")

        st.markdown("---")
        with st.expander("📝 Raw AI Response"):
            st.text(ai.get("raw_response", ""))

# ---- Tab 2: TCM Results ----
with tabs[2]:
    if result is None:
        st.info("Run analysis first.")
    else:
        tcm_df = result["protein_result"].get("tcm")
        if tcm_df is not None and len(tcm_df) > 0:
            tcm_df = tcm_df.copy()
            disp_cols = [c for c in ["HVMID", "cn_name", "pinyin_name", "en_name", "Importance Score"]
                        if c in tcm_df.columns]
            tcm_disp = tcm_df[disp_cols].sort_values("Importance Score", ascending=False)

            st.dataframe(tcm_disp, use_container_width=True, height=400)

            c1, c2 = st.columns(2)
            with c1:
                fig = px.histogram(
                    tcm_disp, x="Importance Score", nbins=30,
                    title="Importance Score Distribution",
                    color_discrete_sequence=["#3c8dbc"],
                )
                st.plotly_chart(fig, use_container_width=True)
            with c2:
                top_n = tcm_disp.head(15)
                y_col = "cn_name" if "cn_name" in top_n.columns else top_n.columns[0]
                fig = px.bar(
                    top_n, x="Importance Score", y=y_col,
                    orientation="h", title="Top 15 TCMs",
                    color_discrete_sequence=["#3c8dbc"],
                )
                fig.update_layout(yaxis=dict(categoryorder="total ascending"))
                st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("No TCM results available.")

# ---- Tab 3: Chemical Results ----
with tabs[3]:
    if result is None:
        st.info("Run analysis first.")
    else:
        chem_df = result["protein_result"].get("chem")
        if chem_df is not None and len(chem_df) > 0:
            chem_df = chem_df.copy()
            disp_cols = [c for c in ["HVCID", "Name", "Molecular_Formula", "Importance Score"]
                        if c in chem_df.columns]
            chem_disp = chem_df[disp_cols].sort_values("Importance Score", ascending=False)

            st.dataframe(chem_disp, use_container_width=True, height=400)

            c1, c2 = st.columns(2)
            with c1:
                fig = px.histogram(
                    chem_disp, x="Importance Score", nbins=30,
                    title="Importance Score Distribution",
                    color_discrete_sequence=["#f39c12"],
                )
                st.plotly_chart(fig, use_container_width=True)
            with c2:
                top_n = chem_disp.head(15)
                y_col = "Name" if "Name" in top_n.columns else top_n.columns[0]
                fig = px.bar(
                    top_n, x="Importance Score", y=y_col,
                    orientation="h", title="Top 15 Chemicals",
                    color_discrete_sequence=["#f39c12"],
                )
                fig.update_layout(yaxis=dict(categoryorder="total ascending"))
                st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("No chemical results available.")

# ---- Tab 4: Protein Results ----
with tabs[4]:
    if result is None:
        st.info("Run analysis first.")
    else:
        prot_df = result["protein_result"].get("proteins")
        if prot_df is not None and len(prot_df) > 0:
            st.dataframe(prot_df, use_container_width=True, height=400)
        else:
            st.info("No protein results available.")

        st.subheader("Degree Summary")
        deg_sum = result.get("degree_summary")
        if deg_sum is not None:
            st.dataframe(deg_sum, use_container_width=True)

        # Also show formula info if in forward mode
        if is_fwd:
            fml_df = result["protein_result"].get("formula")
            if fml_df is not None and len(fml_df) > 0:
                st.subheader("Formulas")
                disp_cols = [c for c in ["HVPID", "name", "Source Document", "Importance Score"]
                            if c in fml_df.columns]
                st.dataframe(fml_df[disp_cols], use_container_width=True, height=300)

# ---- Tab 5: Network Visualization ----
with tabs[5]:
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
                # Build plotly network
                import plotly.graph_objects as go

                # Assign positions with a simple layout
                unique_nodes = list(
                    set(sub_edges["SourceNode"]) | set(sub_edges["TargetNode"])
                )
                node_idx = {n: i for i, n in enumerate(unique_nodes)}
                # Simple grid layout
                cols = int(np.ceil(np.sqrt(len(unique_nodes))))
                pos = {
                    n: (i % cols * 20, i // cols * 20)
                    for n, i in node_idx.items()
                }

                # Color by attribute
                attr_map = dict(zip(sub_nodes["Key"], sub_nodes["Attribute"]))
                color_map = {"TCM": "#1FA9E9", "Chemicals": "#f1c40f", "Proteins": "#e74c3c"}

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
                    color_map.get(attr_map.get(n, ""), "#95a5a6")
                    for n in unique_nodes if n in pos
                ]
                node_text = [
                    f"{n}<br>{attr_map.get(n, '?')}<br>Degree: "
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

                fig = go.Figure(
                    data=[edge_trace, node_trace],
                    layout=go.Layout(
                        title=f"Sub-network ({len(sub_edges)} edges, {len(set(unique_nodes))} nodes)",
                        showlegend=False,
                        hovermode="closest",
                        margin=dict(b=20, l=5, r=5, t=40),
                        xaxis=dict(showgrid=False, zeroline=False, visible=False),
                        yaxis=dict(showgrid=False, zeroline=False, visible=False),
                        height=600,
                    ),
                )
                st.plotly_chart(fig, use_container_width=True)

            st.info(
                "💡 For the **full** interactive network, export to "
                "**Cytoscape** using `Network.csv` and `Type.csv` in "
                f"the output directory. "
                f"Nodes: {len(node_df)}, Edges: {len(network_df)}."
            )

# ---- Tab 6: Formula Filter (reverse mode only) ----
with tabs[6]:
    if result is None:
        st.info("Run analysis first.")
    elif result.get("formula_filter") is None:
        if is_fwd:
            st.info("Formula filtering is available in **Reverse** mode. "
                     "Switch to Reverse mode and specify a formula to filter by.")
        else:
            st.info("No formula was specified. Enter TCM names in the sidebar "
                     "and re-run to see filtered results.")
    else:
        ff = result["formula_filter"]
        c1, c2, c3 = st.columns(3)
        c1.metric("Found", len(ff.get("tcms_found", [])))
        c2.metric("Missing", len(ff.get("tcms_missing", [])))
        c3.metric("Filtered Edges", len(ff["network"]))

        st.subheader("Filtered Network Edges")
        st.dataframe(ff["network"], use_container_width=True)
        st.subheader("Filtered Network Nodes")
        st.dataframe(ff["node"], use_container_width=True)

# ---- Tab 7: Degree Statistics ----
with tabs[7]:
    if result is None:
        st.info("Run analysis first.")
    else:
        node_deg = result.get("node_degree")
        deg_sum = result.get("degree_summary")

        if node_deg is not None and deg_sum is not None:
            st.subheader("Degree Summary by Attribute")
            st.dataframe(deg_sum, use_container_width=True)

            c1, c2 = st.columns(2)
            with c1:
                st.subheader("Degree Distribution Boxplot")
                fig, ax = plt.subplots(figsize=(6, 4))
                attrs = node_deg["Attribute"].unique()
                data_groups = [
                    node_deg[node_deg["Attribute"] == a]["Degree"].values
                    for a in attrs
                ]
                ax.boxplot(data_groups)
                ax.set_xticklabels(list(attrs))
                ax.set_title("Degree Distribution by Attribute")
                ax.set_xlabel("Node Attribute")
                ax.set_ylabel("Degree")
                st.pyplot(fig)
            with c2:
                st.subheader("Top Nodes by Degree")
                top20 = node_deg.nlargest(20, "Degree")
                fig = px.bar(
                    top20, x="Degree", y="Key",
                    orientation="h", color="Attribute",
                    title="Top 20 Nodes by Degree",
                )
                fig.update_layout(yaxis=dict(categoryorder="total ascending"))
                st.plotly_chart(fig, use_container_width=True)

            st.subheader("Full Node Degree Table")
            st.dataframe(
                node_deg.sort_values("Degree", ascending=False),
                use_container_width=True, height=400,
            )
        else:
            st.info("Degree statistics not available.")

# ---- Tab 8: Optimization ----
with tabs[8]:
    if result is None:
        st.info("Run analysis first.")
    elif is_fwd:
        st.info("Combination optimization is available in **Reverse** mode only. "
                 "Switch to Reverse mode to find optimal TCM/formula combinations.")
    else:
        pr = result["protein_result"]
        tcms_opt = pr.get("tcms")
        fmls_opt = pr.get("formulas")

        has_tcms = tcms_opt is not None and len(tcms_opt) > 0
        has_fmls = fmls_opt is not None and len(fmls_opt) > 0

        if not has_tcms and not has_fmls:
            st.info("No optimization results. Enable optimization and re-run.")
        else:
            c1, c2 = st.columns(2)
            with c1:
                st.subheader("TCM Combinations")
                if has_tcms:
                    st.dataframe(tcms_opt, use_container_width=True)
                    fig = px.bar(
                        tcms_opt.head(20), y="Boost",
                        title="TCM Combination Boost",
                        color_discrete_sequence=["#00a65a"],
                    )
                    fig.update_layout(xaxis_title="Combination Rank")
                    st.plotly_chart(fig, use_container_width=True)
                else:
                    st.info("No TCM combination results.")
            with c2:
                st.subheader("Formula Combinations")
                if has_fmls:
                    st.dataframe(fmls_opt, use_container_width=True)
                    fig = px.bar(
                        fmls_opt.head(20), y="Boost",
                        title="Formula Combination Boost",
                        color_discrete_sequence=["#3c8dbc"],
                    )
                    fig.update_layout(xaxis_title="Combination Rank")
                    st.plotly_chart(fig, use_container_width=True)
                else:
                    st.info("No formula combination results.")

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
            file_name="Network.csv",
            mime="text/csv",
            use_container_width=True,
        )
