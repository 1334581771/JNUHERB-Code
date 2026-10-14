"""Map clinical (non-GB/T) syndrome names to GB/T 16751.2-2021 standard terms.

The profiles carry ~530 unique syndrome names that are neither standard terms
nor official aliases (e.g. 食积证, 痰热蕴肺证, 湿热黄疸证).  This script asks
the DeepSeek LLM to pick the closest GB/T standard term for each, choosing
ONLY from a short candidate list pre-computed by token overlap, and abstaining
(empty string) when nothing fits.  Results are checkpointed and only applied
after validation against the official term list.

Usage
-----
    export DEEPSEEK_API_KEY=sk-...
    python scripts/map_syndrome_to_gbt.py                 # run (resumable)
    python scripts/map_syndrome_to_gbt.py --apply         # merge into profiles
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GBT_PATH = os.path.join(ROOT, "src", "jnuherb", "data", "gbt16751_syndromes.json")
PROFILES = os.path.join(ROOT, "src", "jnuherb", "data", "formula_syndrome_profiles.json")
UNMATCHED = os.path.join(ROOT, "data_ingest", "unmatched_patterns.json")
CHECKPOINT = os.path.join(ROOT, "data_ingest", "gbt_map_checkpoint.jsonl")
BACKUP = os.path.join(ROOT, "data_ingest", "formula_syndrome_profiles.backup.json")
OUT_JSON = os.path.join(ROOT, "src", "jnuherb", "data", "gbt16751_syndromes.json")

API_URL = "https://api.deepseek.com/chat/completions"

# Token lists borrowed from jnuherb.chat for candidate pre-scoring
sys.path.insert(0, os.path.join(ROOT, "src"))
from jnuherb.chat import _ZANGFU_TOKENS, _BINGJI_TOKENS  # noqa: E402

PROMPT = """你是中医标准化术语专家，精通 GB/T 16751.2-2021《中医临床诊疗术语 第2部分：证候》。

下面给出若干【临床惯用证型名】，每个都附有从标准中预选的候选标准术语。请为每个惯用名选择**临床含义最接近**的标准术语：

规则：
1. 只能从该惯用名对应的候选列表中选择，不要自行创造术语；
2. 选择临床主治病机最贴近者，不必字面相同（如 食积证 → 饮食积滞证，痰热蕴肺证 → 痰热壅肺证）；
3. 候选列表中没有合适对应时，standard 输出空字符串 ""，不要勉强匹配；
4. 注意区分：证型名中的脏腑定位（心/肝/脾/肺/肾/胃…）和病机（气虚/血瘀/湿热/阴虚…）都要尽量对上。

请严格按以下 JSON 格式回复（不要包含其他文字）：
{{"results": [{{"name": "<临床惯用名>", "standard": "<选中的标准术语或空串>"}}, ...]}}

待处理（共 {count} 个）：
{entries}"""


def candidates_for(name: str, std_terms: list[str], top_n: int = 6) -> list[str]:
    """Pre-select likely standard terms by syndrome-token + char overlap."""
    toks = [t for t in _ZANGFU_TOKENS + _BINGJI_TOKENS if t in name]
    name_chars = set(name)
    scored = []
    for t in std_terms:
        s = 0.0
        for tok in toks:
            if tok in t:
                s += 3.0
        s += 0.15 * len(name_chars & set(t))
        if s > 0:
            scored.append((s, t))
    scored.sort(reverse=True)
    picks = [t for _, t in scored[:top_n]]
    # 兜底：保证列表非空，给模型一点选择余地
    if not picks:
        picks = std_terms[:top_n]
    return picks


def extract_results(text: str) -> list[dict] | None:
    m = re.search(r"```(?:json)?\s*([\s\S]*?)```", text)
    cands = [m.group(1)] if m else []
    m2 = re.search(r"\{[\s\S]*\}", text)
    if m2:
        cands.append(m2.group(0))
    cands.append(text)
    for c in cands:
        try:
            obj = json.loads(c)
        except json.JSONDecodeError:
            continue
        if isinstance(obj, dict) and isinstance(obj.get("results"), list):
            return obj["results"]
        if isinstance(obj, list):
            return obj
    return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--batch-size", type=int, default=25)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--sleep", type=float, default=1.0)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--apply", action="store_true")
    args = ap.parse_args()

    gbt = json.load(open(GBT_PATH, encoding="utf-8"))
    std_terms = gbt["standard_terms"]
    std_set = set(std_terms)

    if args.apply:
        done: dict[str, str] = {}
        if os.path.exists(CHECKPOINT):
            for line in open(CHECKPOINT, encoding="utf-8"):
                r = json.loads(line)
                if r.get("standard") and r["standard"] in std_set:
                    done[r["name"]] = r["standard"]
        gbt.setdefault("clinical_alias_to_standard", {})
        gbt["clinical_alias_to_standard"].update(done)
        with open(OUT_JSON, "w", encoding="utf-8") as f:
            json.dump(gbt, f, ensure_ascii=False, indent=1)

        profiles = json.load(open(PROFILES, encoding="utf-8"))
        if not os.path.exists(BACKUP):
            json.dump(profiles, open(BACKUP, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
        fixed = 0
        for p in profiles:
            pat = (p.get("syndrome_pattern_zh") or "").strip()
            if pat in done:
                if "syndrome_pattern_zh_raw" not in p:
                    p["syndrome_pattern_zh_raw"] = pat
                p["syndrome_pattern_zh"] = done[pat]
                p["pattern_source"] = "gbt_clinical_alias"
                fixed += 1
        with open(PROFILES, "w", encoding="utf-8") as f:
            json.dump(profiles, f, ensure_ascii=False, indent=1)
        print(f"clinical aliases added: {len(done)} | profiles corrected: {fixed}")
        return 0

    unmatched = json.load(open(UNMATCHED, encoding="utf-8"))
    done = set()
    if os.path.exists(CHECKPOINT):
        done = {json.loads(l)["name"] for l in open(CHECKPOINT, encoding="utf-8") if l.strip()}
    todo = [n for n in unmatched if n not in done]
    if args.limit:
        todo = todo[: args.limit]
    print(f"unique unmatched: {len(unmatched)} | done: {len(done)} | to process: {len(todo)}")

    if args.dry_run:
        for n in todo[:10]:
            print(f"  {n} -> {candidates_for(n, std_terms)}")
        return 0

    api_key = os.environ.get("DEEPSEEK_API_KEY", "")
    if not api_key:
        print("ERROR: DEEPSEEK_API_KEY not set")
        return 1

    import requests

    n_batches = (len(todo) + args.batch_size - 1) // args.batch_size
    with open(CHECKPOINT, "a", encoding="utf-8") as ck:
        for bi in range(n_batches):
            batch = todo[bi * args.batch_size : (bi + 1) * args.batch_size]
            entries = "\n".join(
                f"- {n}\n  候选：{' / '.join(candidates_for(n, std_terms))}"
                for n in batch
            )
            prompt = PROMPT.format(count=len(batch), entries=entries)
            raw = None
            for attempt in range(3):
                try:
                    resp = requests.post(
                        API_URL,
                        headers={"Authorization": f"Bearer {api_key}"},
                        json={"model": "deepseek-chat",
                              "messages": [{"role": "user", "content": prompt}],
                              "temperature": 0.1, "max_tokens": 4096},
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
            by_name = {(r.get("name") or "").strip(): (r.get("standard") or "").strip()
                       for r in results if isinstance(r, dict)}
            got = 0
            for n in batch:
                std = by_name.get(n, "")
                if std and std not in std_set:
                    std = ""  # 幻觉防护：只接受官方列表内的术语
                ck.write(json.dumps({"name": n, "standard": std}, ensure_ascii=False) + "\n")
                got += 1
            print(f"  batch {bi+1}/{n_batches}: {got}/{len(batch)} written")
            time.sleep(args.sleep)
    print("Done. Run with --apply to merge.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
