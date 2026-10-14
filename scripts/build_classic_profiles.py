"""Build / review syndrome profiles for famous classic formulas (经典方治理).

Two batches:

A) 补缺 — classic formulas missing from jnuherb entirely (e.g. 柴胡疏肝散).
   The LLM produces a full profile (GB/T pattern, ba_gang, key symptoms,
   contraindications, classical source) AND the composition; --apply then
   appends a real formula row (provenance=classic_llm) plus formula-TCM
   links so the classic works in Forward analysis too.

B) 复核 — formulas present in the DB whose profile pattern diverges from
   the syndrome they are typically queried by (e.g. 血府逐瘀汤 carries
   瘀阻胸胁证 but users ask with 气滞血瘀证).  The LLM re-picks the
   primary GB/T pattern from candidates, given the recorded indications.

Both batches are checkpointed and validated (patterns must canonicalize
to a GB/T standard term; compositions match the herb table with the
processing-prefix fallback from ingest_chp2020).

Usage
-----
    export DEEPSEEK_API_KEY=sk-...
    python scripts/build_classic_profiles.py --dry-run    # see the plan
    python scripts/build_classic_profiles.py              # LLM batches
    python scripts/build_classic_profiles.py --apply      # merge everything
    python -m jnuherb.build_db                            # then rebuild DB
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
sys.path.insert(0, os.path.join(ROOT, "src"))

from map_syndrome_to_gbt import candidates_for, extract_results  # noqa: E402
from ingest_chp2020 import build_herb_index, match_herb  # noqa: E402

DATA = os.path.join(ROOT, "src", "jnuherb", "data")
FORMULA_CSV = os.path.join(DATA, "jnuherb_formula.csv")
TCM_CSV = os.path.join(DATA, "jnuherb_tcm.csv")
LINKS_CSV = os.path.join(DATA, "jnuherb_formula_tcm_links.csv")
PROFILES = os.path.join(DATA, "formula_syndrome_profiles.json")
GBT_PATH = os.path.join(DATA, "gbt16751_syndromes.json")
CHECKPOINT = os.path.join(ROOT, "data_ingest", "classic_profile_checkpoint.jsonl")
BACKUP_DIR = os.path.join(ROOT, "data_ingest")

API_URL = "https://api.deepseek.com/chat/completions"

# 常用经典方清单（汤头歌诀级知名度）。库里已有档案者自动跳过。
CLASSIC_CANDIDATES = [
    "柴胡疏肝散", "补中益气汤", "玉屏风散", "半夏泻心汤", "平胃散",
    "藿香正气散", "二陈汤", "温胆汤", "清气化痰丸", "三子养亲汤",
    "止嗽散", "小青龙汤", "定喘汤", "苏子降气汤", "旋覆代赭汤",
    "补阳还五汤", "失笑散", "十灰散", "咳血方", "小蓟饮子",
    "槐花散", "黄土汤", "大承气汤", "大黄牡丹汤", "麻子仁丸",
    "济川煎", "大柴胡汤", "四逆散", "痛泻要方", "半夏厚朴汤",
    "天台乌药散", "暖肝煎", "厚朴温中汤", "枳实消痞丸", "越鞠丸",
    "金铃子散", "苏合香丸", "紫雪", "至宝丹", "知柏地黄丸",
    "杞菊地黄丸", "麦味地黄丸", "七宝美髯丹", "炙甘草汤", "当归四逆汤",
    "阳和汤", "大建中汤", "小建中汤", "吴茱萸汤", "四逆汤",
    "参附汤", "十全大补汤", "一贯煎", "百合固金汤", "桑杏汤",
    "清燥救肺汤", "杏苏散", "川芎茶调散", "独活寄生汤", "羚角钩藤汤",
    "镇肝熄风汤", "天麻钩藤饮", "大定风珠", "五苓散", "防己黄芪汤",
    "实脾散", "萆薢分清饮", "甘露消毒丹", "八正散", "三仁汤",
    "茵陈蒿汤", "五味消毒饮", "仙方活命饮", "导赤散", "龙胆泻肝汤",
]

# 复核批：已有档案但 pattern 与常见查询证型有争议者（可 --review 追加）
REVIEW_TARGETS = {
    "血府逐瘀汤": "临床常以 气滞血瘀证/胸中血瘀 检索",
    "金匮肾气丸": "临床常以 肾阳虚证 检索，档案现为 肾阳虚水泛证",
}


def norm(s: str) -> str:
    return re.sub(r"\s+", "", str(s)).lower()


PROMPT_CREATE = """你是中医方剂学专家。下面给出若干个中医经典方（出自经典医籍），其中部分在数据库中缺失，需要你补全完整档案。

为每个方剂输出：
1. syndrome_pattern_zh：主治病证的 GB/T 16751.2-2021 规范证型名（以"证"结尾），**必须从候选列表中挑选**，没有完全合适的输出空串；
2. syndrome_pattern_en：规范英文译名（Pattern 结尾）；
3. composition：组成药味（只写药名，炮制名用通用写法如 炙甘草，不要剂量，用、分隔）；
4. key_symptoms：3-6 个主证关键症状；
5. ba_gang：八纲归纳（yin_yang/interior_exterior/cold_heat/deficiency_excess 各一句）;
6. contraindications：禁忌慎用，无则空数组；
7. source_classic：出自医籍（如 医林改错/太平惠民和剂局方/伤寒论，不写《》）。

请严格按以下 JSON 格式回复（不要包含其他文字）：
{{"results": [{{"name": "<方名>", "syndrome_pattern_zh": "...", "syndrome_pattern_en": "...", "composition": "...", "key_symptoms": [...], "ba_gang": {{...}}, "contraindications": [...], "source_classic": "..."}}, ...]}}

待处理（共 {count} 个）：
{entries}"""

PROMPT_REVIEW = """你是中医标准化术语专家。下面每个方剂已有辨证档案，但其"典型证型"与临床常见检索证型存在偏差。请根据该方的功用主治重新核定其**最常用主证**的 GB/T 16751.2-2021 规范证型名：

1. syndrome_pattern_zh 必须从候选列表中挑选，没有合适的输出空串（保留原值）；
2. 只输出 syndrome_pattern_zh / syndrome_pattern_en 两个字段。

请严格按以下 JSON 格式回复（不要包含其他文字）：
{{"results": [{{"name": "<方名>", "syndrome_pattern_zh": "...", "syndrome_pattern_en": "..."}}, ...]}}

待处理（共 {count} 个）：
{entries}"""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--batch-size", type=int, default=10)
    ap.add_argument("--sleep", type=float, default=1.0)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--review", default="", help="追加复核方名，逗号分隔")
    ap.add_argument("--apply", action="store_true")
    args = ap.parse_args()

    import pandas as pd

    gbt = json.load(open(GBT_PATH, encoding="utf-8"))
    std_terms = gbt["standard_terms"]
    std_set = set(std_terms)
    profiles = json.load(open(PROFILES, encoding="utf-8"))
    pmap = {norm(p.get("name_zh", "")): p for p in profiles}
    jf = pd.read_csv(FORMULA_CSV, dtype=str, keep_default_na=False)
    db_names = {norm(n) for n in jf["name"]}

    review = dict(REVIEW_TARGETS)
    for item in (args.review or "").split(","):
        item = item.strip()
        if item:
            review.setdefault(item, "用户指定复核")

    # 补缺批：清单中既不在库、也无档案者
    create = [n for n in CLASSIC_CANDIDATES
              if norm(n) not in db_names and norm(n) not in pmap]
    # 复核批：在库且有档案者
    review_todo = {n: why for n, why in review.items() if norm(n) in pmap}

    print(f"补缺批: {len(create)} 个 | 复核批: {len(review_todo)} 个")
    for n in create[:40]:
        print(f"  [补] {n}")
    for n, why in review_todo.items():
        p = pmap[norm(n)]
        print(f"  [核] {n}: 现 pattern={p.get('syndrome_pattern_zh','')!r} ({why})")

    if args.dry_run:
        return 0

    # ---------------- --apply ----------------
    if args.apply:
        from jnuherb.chat import canonical_syndrome
        done = {}
        if os.path.exists(CHECKPOINT):
            for line in open(CHECKPOINT, encoding="utf-8"):
                r = json.loads(line)
                done[norm(r["name"])] = r

        # 备份
        for path in (PROFILES, FORMULA_CSV, LINKS_CSV):
            tag = os.path.basename(path).replace(".", "_")
            bak = os.path.join(BACKUP_DIR, f"{tag}.pre_classic.bak")
            if not os.path.exists(bak):
                import shutil
                shutil.copy2(path, bak)

        profiles_added = profiles_fixed = rows_added = links_added = rejected = 0
        links = pd.read_csv(LINKS_CSV, dtype=str, keep_default_na=False)
        existing_links = set(zip(links["HVPID"], links["HVMID"]))
        herb_idx = build_herb_index(pd.read_csv(TCM_CSV, dtype=str, keep_default_na=False))

        # 先处理补缺（需要新 HVPID）
        start = jf["HVPID"].str.extract(r"HVP(\d+)")[0].astype(int).max() + 1
        offset = 0
        for name in create:
            rec = done.get(norm(name))
            if not rec:
                continue
            pat = canonical_syndrome((rec.get("syndrome_pattern_zh") or "").strip())
            if pat not in std_set or not rec.get("composition"):
                rejected += 1
                continue
            hvpid = f"HVP{start + offset:05d}"
            offset += 1
            comp_raw = re.split(r"[、,，;；\s]+", rec["composition"].strip())
            comp_clean, unmatched = [], []
            for h in comp_raw:
                h = h.strip()
                if not h:
                    continue
                comp_clean.append(h)
            comp_str = "、".join(dict.fromkeys(comp_clean))
            jf.loc[len(jf)] = {
                "HVPID": hvpid, "name": name,
                "Prescription Composition": comp_str,
                "Instructions": "",
                "Treatment Symptoms": "、".join(rec.get("key_symptoms", [])),
                "Source Document": rec.get("source_classic", ""),
                "Syndromes_zh": pat, "Syndromes_en": rec.get("syndrome_pattern_en", ""),
                "Indications_en": "", "Formula_en_name": "", "Formula_alias": "",
                "Dosage_form": "汤剂", "Formula_type": "经典方",
                "Formula_category": "Classical decoctions",
                "HERB_formula_id": "", "Formula_provenance": "classic_llm",
            }
            for h in comp_clean:
                hvmid = match_herb(h, herb_idx)
                if hvmid is None:
                    unmatched.append(h)
                    continue
                pair = (hvpid, hvmid)
                if pair not in existing_links:
                    existing_links.add(pair)
                    links.loc[len(links)] = list(pair)
                    links_added += 1
            profiles.append({
                "name_zh": name, "name_en": "", "hvpid": hvpid,
                "ba_gang": rec.get("ba_gang", {}), "si_zhen": {},
                "syndrome_pattern_zh": pat,
                "syndrome_pattern_en": rec.get("syndrome_pattern_en", ""),
                "key_symptoms": rec.get("key_symptoms", []),
                "contraindications": rec.get("contraindications", []),
                "source_classic": rec.get("source_classic", ""),
                "pattern_source": "classic_llm",
            })
            profiles_added += 1
            rows_added += 1
            if unmatched:
                print(f"  [补] {name}: 未匹配药味 {unmatched}")

        # 再处理复核
        for name in review_todo:
            rec = done.get(norm(name))
            if not rec:
                continue
            pat = canonical_syndrome((rec.get("syndrome_pattern_zh") or "").strip())
            if not pat or pat not in std_set:
                rejected += 1
                continue
            p = pmap[norm(name)]
            if p.get("syndrome_pattern_zh") == pat:
                continue  # 无变化
            if "syndrome_pattern_zh_raw" not in p:
                p["syndrome_pattern_zh_raw"] = p.get("syndrome_pattern_zh", "")
            p["syndrome_pattern_zh"] = pat
            p["syndrome_pattern_en"] = rec.get("syndrome_pattern_en", p.get("syndrome_pattern_en", ""))
            p["pattern_source"] = "classic_review"
            profiles_fixed += 1
            print(f"  [核] {name}: {p['syndrome_pattern_zh_raw']} -> {pat}")

        with open(PROFILES, "w", encoding="utf-8") as f:
            json.dump(profiles, f, ensure_ascii=False, indent=1)
        jf.to_csv(FORMULA_CSV, index=False)
        links.to_csv(LINKS_CSV, index=False)
        print(f"\nprofiles +{profiles_added} / 修正 {profiles_fixed} | "
              f"formula rows +{rows_added} | links +{links_added} | rejected {rejected}")
        print("Next: python -m jnuherb.build_db")
        return 0

    # ---------------- LLM batches ----------------
    done_names = set()
    if os.path.exists(CHECKPOINT):
        done_names = {norm(json.loads(l)["name"]) for l in open(CHECKPOINT, encoding="utf-8") if l.strip()}
    todo_create = [n for n in create if norm(n) not in done_names]
    todo_review = [n for n in review_todo if norm(n) not in done_names]
    print(f"to create: {len(todo_create)} | to review: {len(todo_review)}")

    api_key = os.environ.get("DEEPSEEK_API_KEY", "")
    if not api_key:
        print("ERROR: DEEPSEEK_API_KEY not set")
        return 1
    import requests

    def run_batches(todo, prompt_tpl, block_fn, tag):
        n = (len(todo) + args.batch_size - 1) // args.batch_size
        with open(CHECKPOINT, "a", encoding="utf-8") as ck:
            for bi in range(n):
                batch = todo[bi * args.batch_size:(bi + 1) * args.batch_size]
                prompt = prompt_tpl.format(count=len(batch),
                                           entries="\n\n".join(block_fn(x) for x in batch))
                raw = None
                for attempt in range(3):
                    try:
                        resp = requests.post(API_URL,
                                             headers={"Authorization": f"Bearer {api_key}"},
                                             json={"model": "deepseek-chat",
                                                   "messages": [{"role": "user", "content": prompt}],
                                                   "temperature": 0.2, "max_tokens": 4096},
                                             timeout=120)
                        resp.raise_for_status()
                        raw = resp.json()["choices"][0]["message"]["content"]
                        break
                    except Exception as e:  # noqa: BLE001
                        print(f"  [{tag}] batch {bi+1}/{n} attempt {attempt+1} failed: {e}")
                        time.sleep(5 * (attempt + 1))
                if raw is None:
                    continue
                by_name = {(r.get("name") or "").strip(): r
                           for r in (extract_results(raw) or []) if isinstance(r, dict)}
                got = 0
                for x in batch:
                    r = by_name.get(x)
                    if not r:
                        continue
                    rec = {"name": x}
                    rec.update({k: r.get(k, "") for k in (
                        "syndrome_pattern_zh", "syndrome_pattern_en", "composition",
                        "key_symptoms", "ba_gang", "contraindications", "source_classic")})
                    ck.write(json.dumps(rec, ensure_ascii=False) + "\n")
                    got += 1
                print(f"  [{tag}] batch {bi+1}/{n}: {got}/{len(batch)}")
                time.sleep(args.sleep)

    def create_block(name):
        return f"- 方剂：{name}"

    def review_block(name):
        p = pmap[norm(name)]
        row = jf[jf["name"].map(norm) == norm(name)]
        ind = str(row.iloc[0].get("Treatment Symptoms", ""))[:150] if len(row) else ""
        cur = p.get("syndrome_pattern_zh", "")
        cands = candidates_for(cur or "气血", std_terms, top_n=6)
        return (f"- 方剂：{name}（现 pattern：{cur}）\n"
                f"  主治/症状：{ind}\n"
                f"  候选：{' / '.join(cands)}")

    if todo_create:
        run_batches(todo_create, PROMPT_CREATE, create_block, "补")
    if todo_review:
        run_batches(todo_review, PROMPT_REVIEW, review_block, "核")
    print("Done. Review the checkpoint, then run with --apply.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
