"""Fetch TCMSP herb->ingredient / ingredient->target / target->disease relations.

Source: old.tcmsp-e.com single-herb detail pages (tcmspsearch.php?qr=<en_name>).
Resumable: caches each herb's three grid payloads as JSON under data_raw/tcmsp/relation_cache/.
Re-running only fetches herbs not yet cached. Duplicate English names are skipped and
logged to data_raw/tcmsp/ambiguous_herbs.csv to avoid cross-contamination.

Outputs (long-format, one row per relation):
  data_raw/tcmsp/tcmsp_herb_ingredient.csv   (herb_en_name, MOL_ID, molecule_name, ob, dl, ...)
  data_raw/tcmsp/tcmsp_ingredient_target.csv (MOL_ID, molecule_name, target_ID, target_name, drugbank_ID, validated, SVM_score, RF_score)
  data_raw/tcmsp/tcmsp_target_disease.csv    (target_ID, target_name, disease_ID, disease_name)
"""
import json
import os
import re
import time
import urllib.parse
import urllib.request

import pandas as pd

BASE = "https://old.tcmsp-e.com"
RAW = "data_raw/tcmsp"
CACHE = os.path.join(RAW, "relation_cache")
os.makedirs(CACHE, exist_ok=True)

HEADERS = {
    "User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
    "(KHTML, like Gecko) Chrome/120.0 Safari/537.36"
}


def get(url, retries=3, timeout=90):
    last = None
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers=HEADERS)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(2 * (i + 1))
    raise RuntimeError(f"GET failed {url}: {last}")


def get_token():
    html = get(f"{BASE}/tcmspsearch.php")
    m = re.search(r'name="token"\s+value="([0-9a-f]+)"', html)
    if not m:
        raise RuntimeError("token not found")
    return m.group(1)


def extract_grid(html, grid_id):
    idx = html.find(f'$("#{grid_id}").kendoGrid')
    if idx < 0:
        return None
    seg = html[idx:]
    m = re.search(r"data:\s*(\[)", seg)
    if not m:
        return None
    start = m.start(1)
    depth = 0
    end = -1
    in_str = False
    esc = False
    for i in range(start, len(seg)):
        c = seg[i]
        if in_str:
            if esc:
                esc = False
            elif c == "\\":
                esc = True
            elif c == '"':
                in_str = False
            continue
        if c == '"':
            in_str = True
        elif c == "[":
            depth += 1
        elif c == "]":
            depth -= 1
            if depth == 0:
                end = i
                break
    if end < 0:
        return None
    return seg[start : end + 1]


def fetch_herb(token, en_name):
    q = urllib.parse.quote(en_name)
    url = f"{BASE}/tcmspsearch.php?qr={q}&qsr=herb_en_name&token={token}"
    html = get(url)
    out = {}
    for gid in ("grid", "grid2", "grid3"):
        js = extract_grid(html, gid)
        out[gid] = json.loads(js) if js else []
    return out


def main():
    herbs = pd.read_csv(os.path.join(RAW, "tcmsp_herbs.csv"))
    # duplicate English names cannot be disambiguated by this endpoint -> skip
    dup_mask = herbs.duplicated("herb_en_name", keep=False)
    ambiguous = herbs[dup_mask][["herb_cn_name", "herb_pinyin", "herb_en_name"]]
    ambiguous.to_csv(os.path.join(RAW, "ambiguous_herbs.csv"), index=False)
    herbs = herbs[~dup_mask].reset_index(drop=True)
    print(f"herbs to fetch: {len(herbs)} (ambiguous skipped: {len(ambiguous)})")

    token = get_token()
    herb_ing, ing_targ, targ_dis = [], [], []
    t0 = time.time()
    for i, row in herbs.iterrows():
        en = row["herb_en_name"]
        cache_path = os.path.join(CACHE, f"{row['herb_pinyin']}.json")
        if os.path.exists(cache_path):
            with open(cache_path, encoding="utf-8") as f:
                data = json.load(f)
        else:
            data = fetch_herb(token, en)
            with open(cache_path, "w", encoding="utf-8") as f:
                json.dump(data, f, ensure_ascii=False)
            time.sleep(0.6)

        for rec in data["grid"]:
            herb_ing.append({"herb_en_name": en, **rec})
        for rec in data["grid2"]:
            ing_targ.append({"herb_en_name": en, **rec})
        for rec in data["grid3"]:
            targ_dis.append({"herb_en_name": en, **rec})

        if (i + 1) % 50 == 0:
            print(f"  {i+1}/{len(herbs)}  elapsed {time.time()-t0:.0f}s")

    pd.DataFrame(herb_ing).to_csv(os.path.join(RAW, "tcmsp_herb_ingredient.csv"), index=False)
    pd.DataFrame(ing_targ).to_csv(os.path.join(RAW, "tcmsp_ingredient_target.csv"), index=False)
    pd.DataFrame(targ_dis).to_csv(os.path.join(RAW, "tcmsp_target_disease.csv"), index=False)
    print(f"herb_ingredient={len(herb_ing)}  ingredient_target={len(ing_targ)}  target_disease={len(targ_dis)}")


if __name__ == "__main__":
    main()
