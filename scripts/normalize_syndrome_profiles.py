"""Batch-normalize formula syndrome profiles with the DeepSeek LLM.

Background
----------
formula_syndrome_profiles.json was seeded largely from "功效/治法" fields, so
many entries carry an effect description (补血 / 利水渗湿 / 寒下 ...) in the
``syndrome_pattern_zh`` slot instead of a standard TCM syndrome name (证型名).
The AI validation layer in jnuherb.chat cannot match those, so this script
asks the LLM to rewrite each affected entry to a standard syndrome name
(following GB/T 16751.2-style naming, always ending with "证").

What it does
------------
1. Selects entries whose ``syndrome_pattern_zh`` is empty or does NOT end
   with "证" (the 证-suffixed entries are left untouched).
2. Sends them to DeepSeek in batches (default 20/call) with whatever
   context the profile already has (name, original pattern, source_classic,
   key_symptoms, ba_gang, si_zhen).
3. Streams results to a checkpoint JSONL so runs are resumable and cheap
   to re-run after a failure.
4. ``--apply`` merges the checkpoint back into the profiles JSON:
   - ``syndrome_pattern_zh``/``_en`` are replaced with the normalized values;
   - the old value is preserved in ``syndrome_pattern_zh_raw``;
   - ``pattern_source`` is set to ``llm_normalized`` (``manual`` otherwise);
   - the original file is backed up next to the output first.

Usage
-----
    export DEEPSEEK_API_KEY=sk-...
    python scripts/normalize_syndrome_profiles.py                 # run (resumable)
    python scripts/normalize_syndrome_profiles.py --dry-run       # selection only
    python scripts/normalize_syndrome_profiles.py --apply         # merge checkpoint
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_INPUT = os.path.join(ROOT, "src", "jnuherb", "data", "formula_syndrome_profiles.json")
DEFAULT_CHECKPOINT = os.path.join(ROOT, "data_ingest", "profile_norm_checkpoint.jsonl")
DEFAULT_OUTPUT = os.path.join(ROOT, "data_ingest", "formula_syndrome_profiles_normalized.json")
DEFAULT_BACKUP = os.path.join(ROOT, "data_ingest", "formula_syndrome_profiles.backup.json")

API_URL = "https://api.deepseek.com/chat/completions"

PROMPT_TEMPLATE = """你是中医标准化术语专家。下面每个方剂档案中的"原典型证型"字段填写的其实是【功效/治法描述】（如 补血、利水渗湿、寒下），而不是规范证型名。

请为每个方剂推断其**主治的标准证型名**（依据中医诊断学/方剂学规范，参照 GB/T 16751.2《中医临床诊疗术语 证候部分》的命名习惯），要求：
1. 证型名必须以"证"字结尾，使用规范结构：脏腑/部位 + 病机 + 证（如：心脾两虚证、脾胃气虚证、湿热黄疸证、风寒束表证、瘀阻脑络证）；
2. 一个方剂只输出一个最主要的【主证】；若该方主治两个并列的经典证候（如归脾汤主治心脾气血两虚与脾不统血），选临床应用最广者；
3. 结合方剂名、原功效描述、出处和关键症状综合判断，不要只同义改写原功效词；
4. 确实无法推断规范证型时，syndrome_pattern_zh 输出空字符串 ""，不要编造；
5. syndrome_pattern_en 给出该证型的规范英文译名（如 Heart-Spleen Qi and Blood Deficiency Pattern）。

原功效描述不要丢弃：把其中蕴含的病机信息（如"补血"→血虚、"利水渗湿"→水湿内停）翻译进证型名。

请严格按以下 JSON 格式回复（不要包含其他文字）：
{{"results": [{{"name": "<方剂中文名>", "syndrome_pattern_zh": "<规范证型名>", "syndrome_pattern_en": "<English pattern name>"}}, ...]}}

待处理方剂（共 {count} 个）：
{entries}"""


def needs_norm(pat: str) -> bool:
    """Entries with an empty pattern or a non-证 suffix need normalization."""
    pat = (pat or "").strip()
    return not pat or not pat.endswith("证")


def entry_context(p: dict) -> str:
    """Compact one-entry context for the prompt."""
    lines = [f"- 方剂：{p.get('name_zh', '')}（{p.get('name_en', '')}）"]
    pat = (p.get("syndrome_pattern_zh") or "").strip()
    if pat:
        lines.append(f"  原典型证型（实为功效描述）：{pat}")
    src = (p.get("source_classic") or "").strip()
    if src:
        lines.append(f"  出处/来源：{src[:120]}")
    ks = [s for s in (p.get("key_symptoms") or []) if s][:6]
    if ks:
        lines.append(f"  关键症状：{'、'.join(ks)}")
    bg = p.get("ba_gang") or {}
    if bg:
        lines.append(
            "  八纲："
            + f"阴阳={bg.get('yin_yang', '')}, 表里={bg.get('interior_exterior', '')}, "
            + f"寒热={bg.get('cold_heat', '')}, 虚实={bg.get('deficiency_excess', '')}"
        )
    return "\n".join(lines)


def extract_json_array(text: str) -> list[dict] | None:
    """Pull the results array out of the model response."""
    text = text.strip()
    # fenced ```json ... ```
    m = re.search(r"```(?:json)?\s*([\s\S]*?)```", text)
    candidates = [m.group(1)] if m else []
    candidates.append(text)
    m2 = re.search(r"\{[\s\S]*\}", text)
    if m2:
        candidates.append(m2.group(0))
    for cand in candidates:
        try:
            obj = json.loads(cand)
        except json.JSONDecodeError:
            continue
        if isinstance(obj, dict) and isinstance(obj.get("results"), list):
            return obj["results"]
        if isinstance(obj, list):
            return obj
    return None


def call_deepseek(api_key: str, model: str, prompt: str, timeout: int = 120) -> str:
    import requests

    resp = requests.post(
        API_URL,
        headers={
            "Authorization": f"Bearer {api_key}",
            "Content-Type": "application/json",
        },
        json={
            "model": model,
            "messages": [{"role": "user", "content": prompt}],
            "temperature": 0.2,
            "max_tokens": 4096,
        },
        timeout=timeout,
    )
    resp.raise_for_status()
    return resp.json()["choices"][0]["message"]["content"]


def load_checkpoint(path: str) -> dict[str, dict]:
    """name_zh -> normalized record, from prior runs."""
    done: dict[str, dict] = {}
    if not os.path.exists(path):
        return done
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue
            done[rec["name"]] = rec
    return done


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--input", default=DEFAULT_INPUT)
    ap.add_argument("--checkpoint", default=DEFAULT_CHECKPOINT)
    ap.add_argument("--output", default=DEFAULT_OUTPUT)
    ap.add_argument("--backup", default=DEFAULT_BACKUP)
    ap.add_argument("--api-key", default=os.environ.get("DEEPSEEK_API_KEY", ""))
    ap.add_argument("--model", default="deepseek-chat")
    ap.add_argument("--batch-size", type=int, default=20)
    ap.add_argument("--limit", type=int, default=0, help="process at most N entries (debug)")
    ap.add_argument("--sleep", type=float, default=1.0, help="seconds between API calls")
    ap.add_argument("--dry-run", action="store_true", help="only report selection")
    ap.add_argument("--apply", action="store_true", help="merge checkpoint into output JSON")
    args = ap.parse_args()

    with open(args.input, encoding="utf-8") as f:
        profiles = json.load(f)

    # ---- --apply: merge checkpoint back into the profiles JSON ------------
    if args.apply:
        done = load_checkpoint(args.checkpoint)
        if not done:
            print("Checkpoint is empty; nothing to apply.")
            return 1
        with open(args.input, encoding="utf-8") as f:
            original = json.load(f)  # pristine snapshot for the backup
        changed = skipped = 0
        for p in profiles:
            rec = done.get(p.get("name_zh", ""))
            if not rec:
                continue
            new_pat = (rec.get("syndrome_pattern_zh") or "").strip()
            if not new_pat or not new_pat.endswith("证"):
                skipped += 1
                continue
            if "syndrome_pattern_zh_raw" not in p:
                p["syndrome_pattern_zh_raw"] = p.get("syndrome_pattern_zh", "")
            p["syndrome_pattern_zh"] = new_pat
            p["syndrome_pattern_en"] = rec.get("syndrome_pattern_en", "")
            p["pattern_source"] = "llm_normalized"
            changed += 1
        for p in profiles:
            p.setdefault("pattern_source", "manual")
        if not os.path.exists(args.backup):
            with open(args.backup, "w", encoding="utf-8") as f:
                json.dump(original, f, ensure_ascii=False, indent=1)
        with open(args.output, "w", encoding="utf-8") as f:
            json.dump(profiles, f, ensure_ascii=False, indent=1)
        print(f"Applied {changed} normalizations ({skipped} rejected), "
              f"backup -> {args.backup}\noutput -> {args.output}")
        return 0

    # ---- selection ---------------------------------------------------------
    todo = [p for p in profiles if needs_norm(p.get("syndrome_pattern_zh", ""))]
    if args.limit:
        todo = todo[: args.limit]
    done = load_checkpoint(args.checkpoint)
    todo = [p for p in todo if p.get("name_zh", "") not in done]

    print(f"total profiles: {len(profiles)} | need normalization: "
          f"{sum(1 for p in profiles if needs_norm(p.get('syndrome_pattern_zh', '')))} "
          f"| already done: {len(done)} | to process now: {len(todo)}")

    if args.dry_run:
        for p in todo[: args.batch_size]:
            print("  -", p.get("name_zh"), "|", (p.get("syndrome_pattern_zh") or "<empty>").strip())
        return 0

    if not args.api_key:
        print("ERROR: no DeepSeek API key. Set DEEPSEEK_API_KEY or pass --api-key.")
        return 1

    os.makedirs(os.path.dirname(args.checkpoint), exist_ok=True)
    n_batches = (len(todo) + args.batch_size - 1) // args.batch_size
    ok = bad = 0
    with open(args.checkpoint, "a", encoding="utf-8") as ck:
        for bi in range(n_batches):
            batch = todo[bi * args.batch_size : (bi + 1) * args.batch_size]
            prompt = PROMPT_TEMPLATE.format(
                count=len(batch),
                entries="\n\n".join(entry_context(p) for p in batch),
            )
            for attempt in range(3):
                try:
                    raw = call_deepseek(args.api_key, args.model, prompt)
                    break
                except Exception as e:  # noqa: BLE001 - resumable, keep going
                    print(f"  batch {bi + 1}/{n_batches} attempt {attempt + 1} failed: {e}")
                    time.sleep(5 * (attempt + 1))
            else:
                bad += len(batch)
                continue
            results = extract_json_array(raw)
            by_name = {
                (r.get("name") or "").strip(): r
                for r in (results or [])
                if isinstance(r, dict)
            }
            got = 0
            for p in batch:
                name = p.get("name_zh", "")
                r = by_name.get(name)
                if not r:
                    continue
                ck.write(json.dumps({
                    "name": name,
                    "syndrome_pattern_zh": (r.get("syndrome_pattern_zh") or "").strip(),
                    "syndrome_pattern_en": (r.get("syndrome_pattern_en") or "").strip(),
                }, ensure_ascii=False) + "\n")
                got += 1
            ok += got
            bad += len(batch) - got
            print(f"  batch {bi + 1}/{n_batches}: {got}/{len(batch)} parsed "
                  f"(cumulative ok={ok}, bad={bad})")
            time.sleep(args.sleep)

    print(f"Done. ok={ok} bad={bad}. Checkpoint: {args.checkpoint}")
    print("Review the checkpoint, then run with --apply to write the merged JSON.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
