"""Generate / fix formula-syndrome profiles for CHP2020 entries with the LLM.

Targets (see chp2020_ingest_report):
  - CHP entries with NO profile yet (newly appended formulas, etc.)
  - CHP entries whose profile pattern is not a GB/T 16751.2 standard term
    (e.g. 抗胆碱药…, 暑湿内蕴寒邪束表证, 太阳少阳合病证)

The prompt feeds the CHP 功能与主治 (rich syndrome/indication text) plus the
composition, and requires the syndrome_pattern_zh to be chosen from GB/T
standard candidates.  Checkpointed + validated before --apply.

Usage
-----
    export DEEPSEEK_API_KEY=sk-...
    python scripts/build_chp_profiles.py               # run (resumable)
    python scripts/build_chp_profiles.py --apply       # merge into profiles
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
from ingest_chp2020 import strip_dose  # noqa: E402

DATA = os.path.join(ROOT, "src", "jnuherb", "data")
CHP_JSON = os.path.join(ROOT, "data_ingest", "chp2020_chengfang.json")
PROFILES = os.path.join(DATA, "formula_syndrome_profiles.json")
GBT_PATH = os.path.join(DATA, "gbt16751_syndromes.json")
CHECKPOINT = os.path.join(ROOT, "data_ingest", "chp_profile_checkpoint.jsonl")
BACKUP = os.path.join(ROOT, "data_ingest", "formula_syndrome_profiles.backup.json")

API_URL = "https://api.deepseek.com/chat/completions"

PROMPT = """你是中医标准化术语专家。下面给出若干个《中国药典》2020年版成方制剂的【功能与主治】和【处方】。请为每个制剂构建辨证档案：

1. syndrome_pattern_zh：该方主治病证的**GB/T 16751.2-2021 规范证型名**（以"证"结尾），必须从给出的候选列表中挑选，没有完全合适的可以输出空串"";
2. syndrome_pattern_en：该证型的规范英文译名（Pattern 结尾）;
3. ba_gang：八纲归纳（yin_yang/interior_exterior/cold_heat/deficiency_excess 四个字段，各一句话）;
4. key_symptoms：3-6 个支持该证型的关键症状（从功能主治的"症见"部分提取）;
5. contraindications：禁忌/慎用情形，没有则空数组。

请严格按以下 JSON 格式回复（不要包含其他文字）：
{{"results": [{{"name": "<制剂名>", "syndrome_pattern_zh": "...", "syndrome_pattern_en": "...", "ba_gang": {{"yin_yang": "...", "interior_exterior": "...", "cold_heat": "...", "deficiency_excess": "..."}}, "key_symptoms": [...], "contraindications": [...]}}, ...]}}

待处理制剂（共 {count} 个）：
{entries}"""


def norm(s: str) -> str:
    return re.sub(r"\s+", "", str(s)).lower()


def entry_block(d: dict, current_pat: str, std_terms: list[str]) -> str:
    lines = [f"- 制剂：{d['title']}"]
    if current_pat:
        lines.append(f"  现有证型（待规范）：{current_pat}")
    zx = str(d.get("功能与主治", "")).strip()
    if zx:
        lines.append(f"  功能与主治：{zx[:220]}")
    cf, _ = strip_dose(str(d.get("处方", "")))
    if cf:
        lines.append(f"  处方：{cf[:120]}")
    cands = candidates_for(current_pat or zx[:30], std_terms, top_n=6)
    lines.append(f"  GB/T候选证型：{' / '.join(cands)}")
    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--batch-size", type=int, default=8)
    ap.add_argument("--sleep", type=float, default=1.0)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--apply", action="store_true")
    args = ap.parse_args()

    gbt = json.load(open(GBT_PATH, encoding="utf-8"))
    std_terms = gbt["standard_terms"]
    std_set = set(std_terms)

    profiles = json.load(open(PROFILES, encoding="utf-8"))
    pmap = {norm(p.get("name_zh", "")): p for p in profiles}
    chp = json.load(open(CHP_JSON, encoding="utf-8"))

    # 目标：无 profile 或 pattern 未国标化的 CHP 条目
    todo = []
    for d in chp:
        key = norm(d["title"])
        p = pmap.get(key)
        if p is None:
            todo.append((d, ""))
        else:
            pat = (p.get("syndrome_pattern_zh") or "").strip()
            if pat and pat not in std_set:
                todo.append((d, pat))
    print(f"targets: {len(todo)} (missing profile: "
          f"{sum(1 for _, c in todo if not c)}, non-standard pattern: "
          f"{sum(1 for _, c in todo if c)})")

    if args.dry_run:
        for d, cur in todo[:5]:
            print(entry_block(d, cur, std_terms))
            print()
        return 0

    # ---- --apply ----------------------------------------------------------
    if args.apply:
        from jnuherb.chat import canonical_syndrome
        done = {}
        if os.path.exists(CHECKPOINT):
            for line in open(CHECKPOINT, encoding="utf-8"):
                r = json.loads(line)
                done[norm(r["name"])] = r
        if not os.path.exists(BACKUP):
            json.dump(profiles, open(BACKUP, "w", encoding="utf-8"),
                      ensure_ascii=False, indent=1)
        added = updated = rejected = 0
        for d, cur in todo:
            rec = done.get(norm(d["title"]))
            if not rec:
                continue
            pat = (rec.get("syndrome_pattern_zh") or "").strip()
            # 幻觉防护：只接受标准术语（或其同义词可解析到标准术语）
            if pat and canonical_syndrome(pat) not in std_set:
                rejected += 1
                continue
            key = norm(d["title"])
            if key in pmap:
                p = pmap[key]
                if "syndrome_pattern_zh_raw" not in p:
                    p["syndrome_pattern_zh_raw"] = p.get("syndrome_pattern_zh", "")
                if pat:
                    p["syndrome_pattern_zh"] = pat
                    p["syndrome_pattern_en"] = rec.get("syndrome_pattern_en", "")
                    p["pattern_source"] = "chp2020_llm"
                    if rec.get("ba_gang"):
                        p["ba_gang"] = rec["ba_gang"]
                    if rec.get("key_symptoms"):
                        p["key_symptoms"] = rec["key_symptoms"]
                    if rec.get("contraindications"):
                        p["contraindications"] = rec["contraindications"]
                    updated += 1
            else:
                if not pat:
                    rejected += 1
                    continue
                profiles.append({
                    "name_zh": d["title"],
                    "name_en": "",
                    "hvpid": "",
                    "ba_gang": rec.get("ba_gang", {}),
                    "si_zhen": {},
                    "syndrome_pattern_zh": pat,
                    "syndrome_pattern_en": rec.get("syndrome_pattern_en", ""),
                    "key_symptoms": rec.get("key_symptoms", []),
                    "contraindications": rec.get("contraindications", []),
                    "source_classic": f"《中国药典》2020年版 一部 第{d.get('pageNum','')}页",
                    "pattern_source": "chp2020_llm",
                })
                added += 1
        with open(PROFILES, "w", encoding="utf-8") as f:
            json.dump(profiles, f, ensure_ascii=False, indent=1)
        print(f"profiles added: {added}, updated: {updated}, rejected: {rejected}")
        return 0

    # ---- run LLM batches ---------------------------------------------------
    done_names = set()
    if os.path.exists(CHECKPOINT):
        done_names = {norm(json.loads(l)["name"]) for l in open(CHECKPOINT, encoding="utf-8") if l.strip()}
    todo = [(d, c) for d, c in todo if norm(d["title"]) not in done_names]
    print(f"to process now: {len(todo)}")

    api_key = os.environ.get("DEEPSEEK_API_KEY", "")
    if not api_key:
        print("ERROR: DEEPSEEK_API_KEY not set")
        return 1

    import requests

    n_batches = (len(todo) + args.batch_size - 1) // args.batch_size
    with open(CHECKPOINT, "a", encoding="utf-8") as ck:
        for bi in range(n_batches):
            batch = todo[bi * args.batch_size : (bi + 1) * args.batch_size]
            prompt = PROMPT.format(
                count=len(batch),
                entries="\n\n".join(entry_block(d, c, std_terms) for d, c in batch),
            )
            raw = None
            for attempt in range(3):
                try:
                    resp = requests.post(
                        API_URL,
                        headers={"Authorization": f"Bearer {api_key}"},
                        json={"model": "deepseek-chat",
                              "messages": [{"role": "user", "content": prompt}],
                              "temperature": 0.2, "max_tokens": 4096},
                        timeout=120,
                    )
                    resp.raise_for_status()
                    raw = resp.json()["choices"][0]["message"]["content"]
                    break
                except Exception as e:  # noqa: BLE001
                    print(f"  batch {bi+1}/{n_batches} attempt {attempt+1} failed: {e}")
                    time.sleep(5 * (attempt + 1))
            if raw is None:
                continue
            results = extract_results(raw) or []
            by_name = {(r.get("name") or "").strip(): r for r in results if isinstance(r, dict)}
            got = 0
            for d, _ in batch:
                r = by_name.get(d["title"])
                if not r:
                    continue
                ck.write(json.dumps({
                    "name": d["title"],
                    "syndrome_pattern_zh": (r.get("syndrome_pattern_zh") or "").strip(),
                    "syndrome_pattern_en": (r.get("syndrome_pattern_en") or "").strip(),
                    "ba_gang": r.get("ba_gang", {}),
                    "key_symptoms": r.get("key_symptoms", []),
                    "contraindications": r.get("contraindications", []),
                }, ensure_ascii=False) + "\n")
                got += 1
            print(f"  batch {bi+1}/{n_batches}: {got}/{len(batch)} parsed")
            time.sleep(args.sleep)
    print("Done. Run with --apply to merge into profiles.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
