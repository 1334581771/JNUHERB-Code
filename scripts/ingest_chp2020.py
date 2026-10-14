"""Ingest Chinese Pharmacopoeia 2020 (CHP2020) 成方制剂 into jnuherb.

Source: data_ingest/chp2020_chengfang.json (1607 entries with 处方 /
功能与主治).  Incremental and additive, following ingest_herb_formula.py
conventions:

  - Matching rows (normalized Chinese name) get their empty
    ``Prescription Composition`` filled from the CHP 处方.  Per project
    decision the gram amounts (大黄56.6g -> 大黄) are stripped: jnuherb
    tracks composition for network analysis, not dosing.  用法与用量 is
    not imported either.
  - Non-matching entries are APPENDED with new HVPIDs and
    provenance=chp2020.
  - Composition herbs are matched against jnuherb_tcm (with 炮制-prefix
    stripping: 炒/炙/蜜炙/煅/炭/...) and new formula-TCM links are added
    so the enriched/newly-added formulas work in Forward analysis.
  - Existing rows are never deleted; non-empty fields are not overwritten.

Outputs (data_ingest/, copy into src/jnuherb/data/ then rebuild the DB):
  - jnuherb_formula_chp_enriched.csv
  - jnuherb_formula_tcm_links_chp.csv
  - chp2020_ingest_report.json
"""

from __future__ import annotations

import json
import os
import re

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "src", "jnuherb", "data")
CHP_JSON = os.path.join(ROOT, "data_ingest", "chp2020_chengfang.json")
OUT_DIR = os.path.join(ROOT, "data_ingest")

FORMULA_CSV = os.path.join(DATA, "jnuherb_formula.csv")
TCM_CSV = os.path.join(DATA, "jnuherb_tcm.csv")
LINKS_CSV = os.path.join(DATA, "jnuherb_formula_tcm_links.csv")

OUT_FORMULA = os.path.join(OUT_DIR, "jnuherb_formula_chp_enriched.csv")
OUT_LINKS = os.path.join(OUT_DIR, "jnuherb_formula_tcm_links_chp.csv")
OUT_REPORT = os.path.join(OUT_DIR, "chp2020_ingest_report.json")

# 炮制/修治前缀，长在前优先；单字前缀（醋/姜/酒/盐/生/焦…）只用于
# 精确匹配失败后的回溯（防止 姜黄->黄 这类误伤，精确命中优先）
PROCESS_PREFIXES = [
    "酒蒸", "酒炖", "酒炙", "酒炒", "蜜炙", "姜炙", "醋炙", "盐炙",
    "油炙", "麸炒", "土炒", "米炒", "煅淬", "炒", "炙", "煅", "蒸",
    "煨", "烫", "炖", "煮", "燀", "制",
]
SINGLE_PREFIXES = ["人工", "体外培育", "醋", "姜", "酒", "盐", "蜜", "油", "土", "麸", "生", "焦", "熟", "清"]
PROCESS_SUFFIXES = ["炭", "焦", "霜", "浓缩粉", "流浸膏", "浸膏", "细粉", "层粉", "提取物", "干膏", "总碱", "总苷", "总皂苷", "仁", "尾", "肉", "粉"]
PAREN_RE = re.compile(r"[（(][^（）()]{1,8}[)）]")
ERRATA_RE = re.compile(r"[■［］\[\]订正]")
MANUAL_HERB_ALIASES = {
    "藿香": "广藿香",
    "人工麝香": "麝香",
    "人工牛黄": "牛黄",
    "体外培育牛黄": "牛黄",
    "酒萸肉": "山茱萸",
    "山银花": "金银花",
    "瓜蒌皮": "瓜蒌",
    "诃子肉": "诃子",
    "光慈姑": "光慈菇",
    "熊胆粉": "熊胆",
    "虎掌南星": "天南星",
}

DOSE_RE = re.compile(r"\d+(?:\.\d+)?\s*(?:g|克|mg|ml|毫升|μg|ug)(?![a-zA-Z])")

DOSAGE_FORMS = [
    ("滴丸", "滴丸剂"), ("橡胶膏", "橡胶膏剂"), ("贴膏", "贴膏剂"),
    ("泡腾片", "泡腾片"), ("分散片", "分散片"), ("含片", "含片"),
    ("咀嚼片", "咀嚼片"), ("缓释片", "缓释片"), ("控释片", "控释片"),
    ("片", "片剂"), ("胶囊", "胶囊剂"), ("软胶囊", "软胶囊剂"),
    ("颗粒", "颗粒剂"), ("丸", "丸剂"), ("散", "散剂"), ("丹", "丹剂"),
    ("膏药", "膏药"), ("膏", "膏剂"), ("糖浆", "糖浆剂"),
    ("口服液", "合剂"), ("合剂", "合剂"), ("酒", "酒剂"), ("酊", "酊剂"),
    ("露", "露剂"), ("茶", "茶剂"), ("锭", "锭剂"), ("栓", "栓剂"),
    ("膜", "膜剂"), ("注射", "注射剂"), ("喷雾", "喷雾剂"),
    ("气雾", "气雾剂"), ("流浸膏", "流浸膏剂"), ("浸膏", "浸膏剂"),
    ("凝胶", "凝胶剂"), ("滴眼", "滴眼剂"), ("滴鼻", "滴鼻剂"),
    ("洗剂", "洗剂"), ("搽剂", "搽剂"), ("酊剂", "酊剂"),
]


def norm_name(s: str) -> str:
    return re.sub(r"\s+", "", str(s)).lower()


def strip_dose(text: str) -> tuple[str, list[str]]:
    """Remove gram/ml amounts from a CHP 处方 block.

    Returns (composition_str, leftover_lines) where leftover_lines still
    contain digits after cleaning (reported for manual review).
    """
    herbs, leftover = [], []
    for ln in str(text).split("\n"):
        ln = ln.strip()
        if not ln:
            continue
        cleaned = DOSE_RE.sub("", ln)
        cleaned = ERRATA_RE.sub("", cleaned)
        cleaned = re.sub(r"适量", "", cleaned).strip()
        cleaned = PAREN_RE.sub("", cleaned).strip()
        #  glued case like 茯苓149佩兰100g：去掉夹在药名间的数字
        cleaned = re.sub(r"(?<=[一-鿿])\d+(?:\.\d+)?(?=[一-鿿])", "", cleaned)
        if re.search(r"\d", cleaned):
            leftover.append(ln)
        if cleaned:
            herbs.append(cleaned)
    # 同味药去重保序
    seen, out = set(), []
    for h in herbs:
        if h not in seen:
            seen.add(h)
            out.append(h)
    return "、".join(out), leftover


def guess_dosage_form(title: str) -> str:
    for key, form in DOSAGE_FORMS:
        if key in title:
            return form
    return ""


def build_herb_index(tcm: pd.DataFrame) -> dict[str, str]:
    """cn_name -> HVMID, plus de-processed variants."""
    idx: dict[str, str] = {}
    for _, row in tcm.iterrows():
        name = str(row["cn_name"]).strip()
        if not name:
            continue
        idx.setdefault(name, row["HVMID"])
        # 去炮制前缀的裸名（炒酸枣仁 -> 酸枣仁）
        bare = name
        changed = True
        while changed and len(bare) >= 3:
            changed = False
            for pre in PROCESS_PREFIXES:
                if bare.startswith(pre) and len(bare) - len(pre) >= 2:
                    bare = bare[len(pre):]
                    changed = True
            for suf in PROCESS_SUFFIXES:
                if bare.endswith(suf) and len(bare) - len(suf) >= 2:
                    bare = bare[:-len(suf)]
                    changed = True
        if bare != name:
            idx.setdefault(bare, row["HVMID"])
    return idx


def match_herb(name: str, herb_idx: dict[str, str]) -> str | None:
    name = name.strip()
    if not name:
        return None
    if name in herb_idx:
        return herb_idx[name]
    if name in MANUAL_HERB_ALIASES and MANUAL_HERB_ALIASES[name] in herb_idx:
        return herb_idx[MANUAL_HERB_ALIASES[name]]
    # 处方里是炮制品名时逐级回溯：炮制前缀 -> 单字前缀 -> 后缀
    for pre in PROCESS_PREFIXES:
        if name.startswith(pre) and name[len(pre):] in herb_idx:
            return herb_idx[name[len(pre):]]
    for pre in SINGLE_PREFIXES:
        if name.startswith(pre) and len(name) - len(pre) >= 2:
            rest = name[len(pre):]
            if rest in herb_idx:
                return herb_idx[rest]
    for suf in PROCESS_SUFFIXES:
        if name.endswith(suf) and len(name) - len(suf) >= 2:
            rest = name[:-len(suf)]
            if rest in herb_idx:
                return herb_idx[rest]
            if rest in MANUAL_HERB_ALIASES and MANUAL_HERB_ALIASES[rest] in herb_idx:
                return herb_idx[MANUAL_HERB_ALIASES[rest]]
    return None


def segment_match(text: str, herb_idx: dict[str, str]) -> list[str]:
    """Greedy longest-match segmentation for glued herb names
    (e.g. 茯苓149佩兰 -> 茯苓+佩兰).  Returns matched HVMIDs."""
    res, i = [], 0
    while i < len(text):
        for L in range(min(7, len(text) - i), 1, -1):
            seg = text[i:i + L]
            if seg in herb_idx:
                res.append(herb_idx[seg])
                i += L
                break
        else:
            i += 1
    return res


def main() -> None:
    chp = json.load(open(CHP_JSON, encoding="utf-8"))
    jf = pd.read_csv(FORMULA_CSV, dtype=str, keep_default_na=False)
    tcm = pd.read_csv(TCM_CSV, dtype=str, keep_default_na=False)
    links = pd.read_csv(LINKS_CSV, dtype=str, keep_default_na=False)

    herb_idx = build_herb_index(tcm)
    jf["_norm"] = jf["name"].map(norm_name)
    name_to_idx = {}
    for i, n in enumerate(jf["_norm"]):
        name_to_idx.setdefault(n, i)

    report = {
        "chp_entries": len(chp),
        "composition_filled": 0,
        "composition_already_present": 0,
        "appended": 0,
        "links_added": 0,
        "herbs_unmatched": {},
        "leftover_dose_lines": {},
    }

    # ---- 1) enrich matching rows / collect append list -------------------
    to_append = []
    chp_compositions: dict[str, str] = {}  # norm name -> composition
    for d in chp:
        key = norm_name(d["title"])
        comp, leftover = strip_dose(d["处方"])
        chp_compositions[key] = comp
        if leftover:
            report["leftover_dose_lines"][d["title"]] = leftover
        if key in name_to_idx:
            idx = name_to_idx[key]
            if not str(jf.at[idx, "Prescription Composition"]).strip():
                jf.at[idx, "Prescription Composition"] = comp
                prov = str(jf.at[idx, "Formula_provenance"])
                jf.at[idx, "Formula_provenance"] = (
                    prov if "chp2020" in prov else (prov + "+chp2020").strip("+")
                )
                if not str(jf.at[idx, "Source Document"]).strip():
                    jf.at[idx, "Source Document"] = "Chinese Pharmacopoeia 2020 edition"
                report["composition_filled"] += 1
            else:
                report["composition_already_present"] += 1
        else:
            to_append.append((d, comp))

    # ---- 2) append new CHP formulas --------------------------------------
    if to_append:
        start = jf["HVPID"].str.extract(r"HVP(\d+)")[0].astype(int).max() + 1
        for offset, (d, comp) in enumerate(to_append):
            jf.loc[len(jf)] = {
                "HVPID": f"HVP{start + offset:05d}",
                "name": d["title"],
                "Prescription Composition": comp,
                "Instructions": "",
                "Treatment Symptoms": str(d.get("功能与主治", "")).strip(),
                "Source Document": "Chinese Pharmacopoeia 2020 edition",
                "Syndromes_zh": "",
                "Syndromes_en": "",
                "Indications_en": "",
                "Formula_en_name": "",
                "Formula_alias": "",
                "Dosage_form": guess_dosage_form(d["title"]),
                "Formula_type": "",
                "Formula_category": "Chinese patent medicines",
                "HERB_formula_id": "",
                "Formula_provenance": "chp2020",
                "_norm": norm_name(d["title"]),
            }
            report["appended"] += 1

    # ---- 3) formula-TCM links from composition ----------------------------
    existing_links = set(zip(links["HVPID"], links["HVMID"]))
    new_links = []
    for _, row in jf.iterrows():
        comp = str(row["Prescription Composition"]).strip()
        if not comp or "chp2020" not in str(row["Formula_provenance"]):
            continue
        for herb in re.split(r"[、,，;；\n]", comp):
            herb = herb.strip()
            if not herb:
                continue
            hvmid = match_herb(herb, herb_idx)
            if hvmid is not None:
                matched_ids = [hvmid]
            elif len(herb) >= 4 and re.search(r"[一-鿿]", herb):
                # 粘连药名（如 茯苓佩兰）做贪心最长匹配
                matched_ids = segment_match(herb, herb_idx)
                if matched_ids:
                    report.setdefault("segmented_herbs", {})[herb] = len(matched_ids)
            else:
                matched_ids = []
            if not matched_ids:
                report["herbs_unmatched"][herb] = report["herbs_unmatched"].get(herb, 0) + 1
                continue
            for hvmid in matched_ids:
                pair = (row["HVPID"], hvmid)
                if pair not in existing_links:
                    existing_links.add(pair)
                    new_links.append(pair)
                    report["links_added"] += 1

    if new_links:
        links = pd.concat(
            [links, pd.DataFrame(new_links, columns=["HVPID", "HVMID"])],
            ignore_index=True,
        )

    # ---- 4) write outputs --------------------------------------------------
    os.makedirs(OUT_DIR, exist_ok=True)
    jf = jf.drop(columns=["_norm"])
    jf.to_csv(OUT_FORMULA, index=False)
    links.to_csv(OUT_LINKS, index=False)
    with open(OUT_REPORT, "w", encoding="utf-8") as f:
        json.dump(report, f, ensure_ascii=False, indent=1)

    print(json.dumps({k: v for k, v in report.items()
                      if k not in ("herbs_unmatched", "leftover_dose_lines",)},
                     ensure_ascii=False, indent=1))
    print(f"unmatched herb names: {len(report['herbs_unmatched'])} "
          f"(top: {sorted(report['herbs_unmatched'].items(), key=lambda x: -x[1])[:8]})")
    print(f"leftover dose lines: {len(report['leftover_dose_lines'])}")
    print(f"\noutputs:\n  {OUT_FORMULA}\n  {OUT_LINKS}\n  {OUT_REPORT}")
    print("\nNext: copy the two CSVs into src/jnuherb/data/ and run "
          "`python -m jnuherb.build_db` to refresh jnuherb.db.")


if __name__ == "__main__":
    main()
