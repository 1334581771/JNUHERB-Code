"""Ingest HERB formula data into jnuherb as an INCREMENTAL, additive layer.

Principles (website already live at jnuherb.top):
  - jnuherb original rows are NEVER modified or removed.
  - Existing jnuherb formulas keep their HVPID; matching HERB formulas are only
    enriched in-place via the new *enrichment* columns (added at the END).
  - HERB formulas not already present (by normalized Chinese name) are APPENDED with
    new HVPIDs and marked provenance=herb.
  - Dedup: normalized Chinese name. jnuherb names win on collision.

Outputs (additive; new columns at end so old readers still work):
  - data_ingest/jnuherb_formula_enriched.csv   (copy-in replacement for data CSV)
  - data_ingest/jnuherb_formula.csv            (cleaned original-shaped view)
  - data_ingest/formula_syndrome_profiles.json (extended AI syndrome profiles)

New enrichment columns:
  Syndromes_zh, Syndromes_en, Indications_en, Formula_en_name, Formula_alias,
  Dosage_form, Formula_type, Formula_category, HERB_formula_id, Formula_provenance
"""
import json
import os
import re

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
JNU_FORMULA = os.path.join(ROOT, "src", "jnuherb", "data", "jnuherb_formula.csv")
HERB_FORMULA = os.path.join(ROOT, "data_raw", "herb", "HERB_formula_info_v2.txt")
OUT_DIR = os.path.join(ROOT, "data_ingest")
PROFILE_PATH = os.path.join(ROOT, "src", "jnuherb", "data", "formula_syndrome_profiles.json")

ENRICH_COLS = [
    "Syndromes_zh", "Syndromes_en", "Indications_en", "Formula_en_name",
    "Formula_alias", "Dosage_form", "Formula_type", "Formula_category",
    "HERB_formula_id", "Formula_provenance",
]


def norm_name(s):
    return re.sub(r"\s+", "", str(s)).lower()


def next_hvpids(start, count):
    return [f"HVP{i:05d}" for i in range(start, start + count)]


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    # load exactly as jnuherb does (dtype=str) so empty stays ""
    jf = pd.read_csv(JNU_FORMULA, dtype=str, keep_default_na=False)
    hf = pd.read_csv(HERB_FORMULA, sep=None, engine="python")

    for c in ENRICH_COLS:
        if c not in jf.columns:
            jf[c] = ""
    jf["Formula_provenance"] = "jnuherb"

    # index HERB by normalized cn name
    hf = hf.copy()
    hf["_norm"] = hf["Formula_cn_name"].map(norm_name)
    hf = hf[hf["_norm"].notna() & (hf["_norm"] != "") & (hf["_norm"] != "nan")]
    hf = hf.drop_duplicates("_norm", keep="first")
    herb_by_norm = hf.set_index("_norm")

    jf["_norm"] = jf["name"].map(norm_name)
    existing = set(jf["_norm"])

    # 1) in-place enrichment of matching jnuherb formulas
    matched = 0
    for idx, row in jf.iterrows():
        key = row["_norm"]
        if key in herb_by_norm.index:
            h = herb_by_norm.loc[key]
            jf.at[idx, "Syndromes_zh"] = "" if pd.isna(h["Syndromes_in_Chinese"]) else str(h["Syndromes_in_Chinese"])
            jf.at[idx, "Syndromes_en"] = "" if pd.isna(h["Syndromes_in_English"]) else str(h["Syndromes_in_English"])
            jf.at[idx, "Indications_en"] = "" if pd.isna(h["Indications_in_English"]) else str(h["Indications_in_English"])
            jf.at[idx, "Formula_en_name"] = "" if pd.isna(h["Formula_en_name"]) else str(h["Formula_en_name"])
            jf.at[idx, "Formula_alias"] = "" if pd.isna(h["Formula_alias_name"]) else str(h["Formula_alias_name"])
            jf.at[idx, "Dosage_form"] = "" if pd.isna(h["Dosage_form"]) else str(h["Dosage_form"])
            jf.at[idx, "Formula_type"] = "" if pd.isna(h["Type"]) else str(h["Type"])
            jf.at[idx, "Formula_category"] = "" if pd.isna(h["Category"]) else str(h["Category"])
            jf.at[idx, "HERB_formula_id"] = "" if pd.isna(h["Formula_id"]) else str(h["Formula_id"])
            jf.at[idx, "Formula_provenance"] = "jnuherb+herb"
            matched += 1

    # 2) append new HERB formulas not already present
    new_mask = ~herb_by_norm.index.isin(existing)
    new_rows = herb_by_norm[new_mask]
    start = jf["HVPID"].str.extract(r"HVP(\d+)")[0].astype(int).max() + 1
    new_hvpids = next_hvpids(start, len(new_rows))

    appended = 0
    for (key, h), hvp in zip(new_rows.iterrows(), new_hvpids):
        jf.loc[len(jf)] = {
            "HVPID": hvp,
            "name": "" if pd.isna(h["Formula_cn_name"]) else str(h["Formula_cn_name"]),
            "Prescription Composition": "",
            "Instructions": "",
            "Treatment Symptoms": "",
            "Source Document": "HERB" if str(h.get("Source", "")).upper().startswith("HERB") else str(h.get("Source", "") or ""),
            "Syndromes_zh": "" if pd.isna(h["Syndromes_in_Chinese"]) else str(h["Syndromes_in_Chinese"]),
            "Syndromes_en": "" if pd.isna(h["Syndromes_in_English"]) else str(h["Syndromes_in_English"]),
            "Indications_en": "" if pd.isna(h["Indications_in_English"]) else str(h["Indications_in_English"]),
            "Formula_en_name": "" if pd.isna(h["Formula_en_name"]) else str(h["Formula_en_name"]),
            "Formula_alias": "" if pd.isna(h["Formula_alias_name"]) else str(h["Formula_alias_name"]),
            "Dosage_form": "" if pd.isna(h["Dosage_form"]) else str(h["Dosage_form"]),
            "Formula_type": "" if pd.isna(h["Type"]) else str(h["Type"]),
            "Formula_category": "" if pd.isna(h["Category"]) else str(h["Category"]),
            "HERB_formula_id": "" if pd.isna(h["Formula_id"]) else str(h["Formula_id"]),
            "Formula_provenance": "herb",
            "_norm": key,
        }
        appended += 1

    # 3) fold syndrome/indication into Treatment Symptoms for AI keyword retrieval
    def fold(row):
        parts = [row["Treatment Symptoms"]]
        for c in ("Syndromes_zh", "Syndromes_en", "Indications_en"):
            v = row.get(c, "")
            if v:
                parts.append(v)
        return "; ".join([p for p in parts if p])
    jf["Treatment Symptoms"] = jf.apply(fold, axis=1)

    # 4) build extended syndrome profiles for AI
    profiles = []
    if os.path.exists(PROFILE_PATH):
        try:
            profiles = json.load(open(PROFILE_PATH, encoding="utf-8"))
        except Exception:
            profiles = []
    existing_prof = {p.get("name_zh", "") for p in profiles}

    def split_terms(s):
        if not s or pd.isna(s):
            return []
        parts = re.split(r"[;；]", str(s))
        return [p.strip() for p in parts if p.strip()]

    for _, row in jf.iterrows():
        name = row["name"]
        if not name or name in existing_prof:
            continue
        szh = split_terms(row["Syndromes_zh"])
        sen = split_terms(row["Syndromes_en"])
        if not szh and not sen:
            continue
        # take primary (first) syndrome as pattern
        profiles.append({
            "name_zh": name,
            "name_en": row.get("Formula_en_name", "") or "",
            "hvpid": row.get("HVPID", "") or "",
            "ba_gang": {},
            "si_zhen": {},
            "syndrome_pattern_zh": szh[0] if szh else "",
            "syndrome_pattern_en": sen[0] if sen else "",
            "key_symptoms": szh[1:] + sen[1:],
            "contraindications": [],
            "source_classic": row.get("Source Document", "") or "",
        })
        existing_prof.add(name)

    # outputs
    out_enriched = os.path.join(OUT_DIR, "jnuherb_formula_enriched.csv")
    jf.drop(columns=["_norm"]).to_csv(out_enriched, index=False)
    out_clean = os.path.join(OUT_DIR, "jnuherb_formula.csv")
    jf[["HVPID", "name", "Prescription Composition", "Instructions",
        "Treatment Symptoms", "Source Document"]].to_csv(out_clean, index=False)
    out_prof = os.path.join(OUT_DIR, "formula_syndrome_profiles.json")
    json.dump(profiles, open(out_prof, "w", encoding="utf-8"), ensure_ascii=False, indent=2)

    # stats
    total = len(jf)
    n_synd = (jf["Syndromes_zh"] != "").sum()
    n_ind = (jf["Indications_en"] != "").sum()
    print(f"jnuherb original rows: {len(existing)}")
    print(f"enriched in place (matched HERB): {matched}")
    print(f"appended new HERB formulas: {appended}")
    print(f"total formula rows now: {total}")
    print(f"  with Syndromes_zh: {n_synd} ({n_synd/total:.1%})")
    print(f"  with Indications_en: {n_ind} ({n_ind/total:.1%})")
    print(f"syndrome profiles total: {len(profiles)}")
    print(f"wrote: {out_enriched}")
    print(f"wrote: {out_clean}")
    print(f"wrote: {out_prof}")


if __name__ == "__main__":
    main()
