#!/usr/bin/env python
"""Regression check: reproduce the three manuscript cases with the current code.

Cases (parameters from JNUHERB_SCI_manuscript / audit report):
  1. Gui Pi Decoction  HVP5811, score 900          (forward)
  2. TNF / IL6 / TP53, score 900                   (reverse, gene_name)
  3. Lung cancer MONDO:0008903, evidence >= 0.1, top-30 seeds, score 900
     (reverse; computed with BOTH web-default and manuscript optimization
     parameters because the first-ranked combination differs between them)

Writes a JSON report to stdout-friendly file for comparison with the
manuscript numbers.
"""

from __future__ import annotations

import json
import sys
import tempfile
import traceback
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from jnuherb.analysis import from_proteins, from_tcm_or_formula, resolve_disease_to_targets
from jnuherb.enrichment import run_enrichment

OUT = Path(__file__).resolve().parent.parent / "build" / "case_verification.json"


def _counts(res: dict) -> dict:
    def n(key):
        df = res.get(key)
        return 0 if df is None else len(df)

    return {
        "formulas": n("formula"),
        "tcm": n("tcm"),
        "chem": n("chem"),
        "proteins": n("proteins"),
    }


def _combo(res: dict, key: str) -> dict | None:
    df = res.get(key)
    if df is None or len(df) == 0:
        return None
    row = df.iloc[0].to_dict()
    return {k: (float(v) if isinstance(v, (int, float)) else str(v))
            for k, v in row.items()}


def _network_stats(path: Path) -> dict:
    net = path / "Network.csv"
    typ = path / "Type.csv"
    nodes = edges = 0
    if net.exists():
        edges = sum(1 for _ in open(net)) - 1
    if typ.exists():
        nodes = sum(1 for _ in open(typ)) - 1
    return {"nodes": nodes, "edges": edges}


def case_gui_pi() -> dict:
    out = Path(tempfile.mkdtemp(prefix="verify_gui_pi_"))
    res = from_tcm_or_formula(
        tcm_or_formula_id=["HVP5811"], score=900, id_type="HVMID",
        path=str(out), re=True, out_for_cytoscape=True, out_graph=True,
    )
    info = {"counts": _counts(res), "network": _network_stats(out)}

    genes = sorted(res["proteins"]["gene_name"].dropna().astype(str).unique().tolist())
    info["n_unique_genes_for_enrichment"] = len(genes)
    enr = run_enrichment(genes)
    for name, block in (("kegg", enr["kegg"]), ("disease", enr["do"])):
        top = block["top"]
        if len(top):
            row = top.iloc[0]
            info[f"{name}_top1"] = {
                "term": str(row.get("Term", "")),
                "overlap": str(row.get("Overlap", "")),
                "adjusted_p": float(row.get("Adjusted P-value", float("nan"))),
                "count": int(row.get("Count", 0)),
            }
        info[f"{name}_n_mapped"] = f"{block['n_mapped']} / {block['n_input']}"
    return info


def case_tnf() -> dict:
    out = Path(tempfile.mkdtemp(prefix="verify_tnf_"))
    res = from_proteins(
        proteins=["TNF", "IL6", "TP53"], score=900, id_type="gene_name",
        num=100, random_state=138192, path=str(out),
    )
    return {"counts": _counts(res), "network": _network_stats(out)}


def case_lung() -> dict:
    resolved = resolve_disease_to_targets("MONDO:0008903", min_score=0.1, max_targets=30)
    seeds = list(resolved["targets"])
    info = {
        "disease_id": resolved["disease_id"],
        "n_seeds_requested": 30,
        "n_seeds_returned": len(seeds),
        "n_total_associated": resolved["n_total"],
    }
    for label, num, seed in (("web_default", 100, 138192),
                             ("manuscript", 1000, 42)):
        out = Path(tempfile.mkdtemp(prefix=f"verify_lung_{label}_"))
        res = from_proteins(
            proteins=seeds, score=900, id_type="Ensembl_ID",
            num=num, random_state=seed, path=str(out),
        )
        info[label] = {
            "num": num, "random_state": seed,
            "counts": _counts(res),
            "tcm_combo_top1": _combo(res, "tcms"),
            "formula_combo_top1": _combo(res, "formulas"),
        }
    return info


CASES = {
    "gui_pi_decoction_forward": (
        case_gui_pi,
        {"counts": {"formulas": 1, "tcm": 8, "chem": 95, "proteins": 1347},
         "network": {"nodes": 1450, "edges": 4358},
         "kegg_top1_term": "Neuroactive ligand-receptor interaction",
         "kegg_top1_adj_p": 4.28e-53,
         "disease_top1_term": "Seizures",
         "disease_top1_adj_p": 4.05e-45},
    ),
    "tnf_il6_tp53_reverse": (
        case_tnf,
        {"counts": {"chem": 38, "tcm": 841, "formulas": 7194},
         "manuscript_formulas": 5780},  # known discrepancy, see report
    ),
    "lung_cancer_reverse": (
        case_lung,
        {"web_default_counts": {"proteins": 10, "chem": 30, "tcm": 745,
                                "formulas": 7125},
         "manuscript_tcm_combo_score": 0.5073067,
         "manuscript_tcm_combo_boost": 70.1916,
         "manuscript_formula_combo_score": 0.8315398,
         "manuscript_formula_combo_boost": 40.8606},
    ),
}


def main() -> None:
    report = {"generated_at": datetime.now(timezone.utc).isoformat(), "cases": {}}
    for name, (fn, expected) in CASES.items():
        print(f"[verify] running {name} ...", flush=True)
        try:
            got = fn()
            report["cases"][name] = {"status": "ok", "got": got,
                                     "expected": expected}
        except Exception:
            report["cases"][name] = {"status": "error",
                                     "traceback": traceback.format_exc()}
        print(f"[verify] {name} -> {report['cases'][name]['status']}", flush=True)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, ensure_ascii=False, indent=2))
    print(f"[verify] report written to {OUT}", flush=True)


if __name__ == "__main__":
    main()
