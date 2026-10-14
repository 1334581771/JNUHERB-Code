"""Ingest CTD / Open Targets gene–disease data into jnuherb (staged).

Stages (run in order; each is idempotent and writes to data_ingest/):

  1   OT direct associations + CTD curated associations
      → jnuherb_protein_disease_links.csv   (filtered to jnuherb proteins)
  2   CTD aggregate associations (3.1 GB gz, streamed line-by-line)
      → merged into the same links CSV
  3   OT indirect (EFO-propagated) associations — optional
  filter  Drop EFO measurement/trait terms (therapeuticAreas == EFO:0001444
      only) — half the OT index is metabolite/phenotype traits, not diseases
  gmt Regenerate disease_ctdot.gmt from the links CSV
      (top-N genes per disease; replaces the 90-term do.gmt for enrichment)
  apply  Copy the links CSV into src/jnuherb/data/ (gzip-compressed, with
      .bak backup)

External reference needed for stage 1:
  data_ingest/raw/gene2ensembl.gz   (NCBI; maps ENSG ↔ ENSP ↔ GeneID)

Raw inputs (data_ingest/raw/):
  CTD_genes_diseases.csv.gz, CTD_curated_genes_diseases.csv.gz
  ot_direct/*.parquet, ot_indirect/*.parquet   (indirect = optional stage 3)

Output schema (one row per protein–disease–source):
  Ensembl_ID, gene_name, disease_id, disease_name, source,
  evidence_type, score, pmids

Design notes:
  * OT parquet is read column-projected (timeseries struct skipped) and
    batch-filtered before ever materialising a full partition in memory.
  * Everything is filtered to the 14,184 jnuherb proteins (ENSP) first,
    which collapses OT-direct from ~7.2 M rows to a few hundred thousand.
  * Disease vocabularies are NOT merged: CTD keeps MESH:/OMIM: ids, OT
    keeps EFO/Orphanet/MONDO ids as-is; disease_id is always prefixed.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

RAW = ROOT / "data_ingest" / "raw"
OUT = ROOT / "data_ingest"
DATA = ROOT / "src" / "jnuherb" / "data"

LINKS_CSV = OUT / "jnuherb_protein_disease_links.csv"
REPORT = OUT / "ctd_ot_ingest_report.json"
GMT_OUT = DATA / "enrichment" / "disease_ctdot.gmt"

OT_SCORE_MIN = 0.05          # OT direct/indirect association score cutoff
CTD_INFER_MIN = 40.0         # CTD aggregate InferenceScore cutoff (0–100+)
# CTD inference rows (weakest evidence: chemical-proxied) are capped per
# disease — the two consumers of this table are (a) Reverse disease seeding,
# which takes top-N by score, and (b) the enrichment GMT (top-N per disease),
# so deep per-disease tails are dead weight: 3.7 M rows / ~2 GB pandas memory
# for zero functional gain.
CTD_INFER_TOP_PER_DISEASE = 100
GENE2ENSEMBL = RAW / "gene2ensembl.gz"

# OT diseaseId arrives as EFO_0000180 (underscore); normalise to the CURIE
# form EFO:0000180 so the whole table uses one "PREFIX:id" convention.
# HP_/OBA_/VT_ entries are phenotypes/traits, not diseases — dropped.
_OT_DISEASE_PREFIXES = ("EFO:", "MONDO:", "Orphanet:")


def _norm_disease_id(raw: str) -> str | None:
    d = raw.strip().replace("_", ":", 1)
    return d if d.startswith(_OT_DISEASE_PREFIXES) else None


def load_protein_index() -> tuple[set[str], dict[str, str], dict[str, str]]:
    """(ENSP set, {ENSP: gene_name}, {GENE_SYMBOL: ENSP}).

    The symbol index upper-cases gene_name — CTD/OT both speak HGNC
    symbols, jnuherb stores them as-is, so the join is case-normalised.
    """
    import pandas as pd
    df = pd.read_csv(DATA / "jnuherb_proteins.csv", dtype=str, keep_default_na=False)
    df = df[df["Ensembl_ID"].str.startswith("ENSP")]
    sym2ensp: dict[str, str] = {}
    for ensp, sym in zip(df["Ensembl_ID"], df["gene_name"]):
        sym2ensp.setdefault(sym.upper(), ensp)
    return set(df["Ensembl_ID"]), dict(zip(df["Ensembl_ID"], df["gene_name"])), sym2ensp


def load_ot_target_symbols() -> dict[str, str]:
    """ENSG → approvedSymbol from the OT target-index parquet.

    data_ingest/raw/ot_targets/*.parquet (download "Target index" from the
    same Open Targets page as the association data).  Aborts stage 1 when
    absent — OT rows cannot be mapped without it.
    """
    import pyarrow.parquet as pq
    files = sorted(RAW.glob("ot_targets*.parquet"))
    tgt_dir = RAW / "ot_targets"
    if tgt_dir.is_dir():
        files += sorted(tgt_dir.glob("*.parquet"))
    if not files:
        sys.exit(f"missing OT target index — download 'Target index' from the "
                 f"Open Targets downloads page into {RAW}/ot_targets/")
    out: dict[str, str] = {}
    for f in files:
        t = pq.read_table(str(f), columns=["id", "approvedSymbol"])
        for i, s in zip(t["id"].to_pylist(), t["approvedSymbol"].to_pylist()):
            out[i] = s
    return out


def load_gene2ensembl(ensp_set: set[str]) -> tuple[dict[str, str], dict[str, str]]:
    """Build ENSG→ENSP and GeneID→ENSP maps restricted to jnuherb proteins.

    gene2ensembl.gz column layout (tab-separated):
      tax_id, GeneID, Ensembl_gene_identifier, RNA_acc, Ensembl_rna,
      protein_acc.version, Ensembl_protein_identifier
    """
    ensg2ensp: dict[str, str] = {}
    geneid2ensp: dict[str, str] = {}
    if not GENE2ENSEMBL.exists():
        print("note: gene2ensembl.gz absent — ENSG join needs the OT target "
              "index fallback, GeneID join disabled")
        return ensg2ensp, geneid2ensp
    try:
        with gzip.open(GENE2ENSEMBL, "rt") as fh:
            fh.readline()  # header
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 7 or parts[0] != "9606":   # human only
                    continue
                ensg = parts[2]
                ensp = parts[6].split(".")[0]
                if ensp in ensp_set:
                    ensg2ensp.setdefault(ensg, ensp)
                    geneid2ensp.setdefault(parts[1], ensp)
    except (EOFError, gzip.BadGzipFile):
        print("note: gene2ensembl.gz incomplete/corrupt — using whatever "
              "was parsed before the truncation point")
    return ensg2ensp, geneid2ensp


def build_ensg2ensp(sym2ensp: dict[str, str],
                    ensg2ensp_g2e: dict[str, str]) -> dict[str, str]:
    """ENSG → ENSP: prefer the direct gene2ensembl map; fill any gaps via
    the OT target index (ENSG→approvedSymbol → jnuherb symbol→ENSP)."""
    ensg2ensp = dict(ensg2ensp_g2e)
    missing_sample = None
    try:
        ensg2sym = load_ot_target_symbols()
    except SystemExit:
        ensg2sym = {}
    bridged = 0
    for ensg, sym in ensg2sym.items():
        if ensg in ensg2ensp:
            continue
        ensp = sym2ensp.get(sym.upper())
        if ensp:
            ensg2ensp[ensg] = ensp
            bridged += 1
        elif missing_sample is None:
            missing_sample = (ensg, sym)
    print(f"ENSG→ENSP: {len(ensg2ensp_g2e):,} via gene2ensembl"
          + (f", +{bridged:,} via symbol bridge" if bridged else ""))
    if missing_sample:
        print(f"  e.g. OT target {missing_sample[0]} ({missing_sample[1]}) "
              "not in jnuherb — expected for paralogs/non-human imports")
    return ensg2ensp


def load_ot_disease_names() -> dict[str, str]:
    """ disease_id → name from an OT disease-index file if present.

    Accepts either the OT Platform ``diseases`` parquet export
    (data_ingest/raw/ot_diseases.parquet, columns id/name) or a plain
    two-column TSV (disease_id<TAB>name).  Returns {} when absent —
    links still work, only disease_name stays empty.
    """
    pq_path = RAW / "ot_diseases.parquet"
    tsv_path = RAW / "ot_diseases.tsv"
    if pq_path.exists():
        import pyarrow.parquet as pq
        t = pq.read_table(str(pq_path), columns=["id", "name"])
        return {_norm_disease_id(i) or i: n
                for i, n in zip(t["id"].to_pylist(), t["name"].to_pylist())}
    if tsv_path.exists():
        out = {}
        with open(tsv_path, encoding="utf-8") as fh:
            for line in fh:
                did, _, name = line.rstrip("\n").partition("\t")
                if did:
                    out[_norm_disease_id(did) or did] = name
        return out
    print("note: no ot_diseases.parquet/tsv found — OT rows get empty "
          "disease_name (disease lookup by name will only cover CTD terms)")
    return {}


def stage1() -> None:
    """OT direct + CTD curated → links CSV."""
    import pyarrow.parquet as pq
    import pyarrow.compute as pc

    ensp_set, ensp2gene, sym2ensp = load_protein_index()
    ensg2ensp_g2e, geneid2ensp = load_gene2ensembl(ensp_set)
    ensg2ensp = build_ensg2ensp(sym2ensp, ensg2ensp_g2e)
    ot_names = load_ot_disease_names()
    print(f"jnuherb proteins: {len(ensp_set)} | GeneID mapped: "
          f"{len(geneid2ensp):,} | OT disease names: {len(ot_names)}")

    report: dict[str, int] = {}
    rows: dict[tuple[str, str, str], dict] = {}

    # ---- OT direct --------------------------------------------------------
    cols = ["diseaseId", "targetId", "associationScore", "evidenceCount"]
    ot_seen = ot_kept = ot_dropped_prefix = 0
    for f in sorted(RAW.glob("ot_direct/*.parquet")):
        try:
            t = pq.read_table(f, columns=cols)
        except Exception as e:  # noqa: BLE001
            print(f"  skip unreadable partition {f.name}: {e}")
            continue
        ot_seen += t.num_rows
        t = t.filter(pc.greater_equal(t["associationScore"], OT_SCORE_MIN))
        d_ids = t["diseaseId"].to_pylist()
        t_ids = t["targetId"].to_pylist()
        scores = t["associationScore"].to_pylist()
        evs = t["evidenceCount"].to_pylist()
        for d_raw, g, s, ev in zip(d_ids, t_ids, scores, evs):
            d = _norm_disease_id(d_raw)
            if d is None:
                ot_dropped_prefix += 1
                continue
            ensp = ensg2ensp.get(g)
            if ensp is None:
                continue
            key = (ensp, d, "OT")
            if key in rows and rows[key]["score"] >= s:
                continue
            rows[key] = {
                "Ensembl_ID": ensp, "gene_name": ensp2gene[ensp],
                "disease_id": d, "disease_name": ot_names.get(d, ""),
                "source": "OT", "evidence_type": f"direct;ev={ev}",
                "score": round(s, 4), "pmids": "",
            }
            ot_kept += 1
    report["ot_direct_rows_seen"] = ot_seen
    report["ot_direct_rows_kept"] = ot_kept
    report["ot_dropped_phenotype_terms"] = ot_dropped_prefix
    report["ot_diseases_with_names"] = sum(1 for r in rows.values()
                                           if r["source"] == "OT" and r["disease_name"])
    print(f"OT direct: {ot_seen:,} seen → {ot_kept:,} kept "
          f"(score≥{OT_SCORE_MIN}, jnuherb targets only, "
          f"{ot_dropped_prefix:,} phenotype terms dropped)")

    # ---- CTD curated ------------------------------------------------------
    ctd_kept = 0
    with gzip.open(RAW / "CTD_curated_genes_diseases.csv.gz", "rt") as fh:
        rdr = csv.reader(filter(lambda ln: not ln.startswith("#"), fh))
        next(rdr, None)  # column-name line
        for parts in rdr:
            if len(parts) < 7:
                continue
            symbol, geneid, dname, did, dev, _omim, pmids = parts[:7]
            # symbol join is primary; GeneID join is the fallback for rows
            # whose symbol is an alias jnuherb does not carry
            ensp = sym2ensp.get(symbol.upper()) or geneid2ensp.get(geneid)
            if ensp is None:
                continue
            key = (ensp, did, "CTD")
            if key in rows:
                continue
            rows[key] = {
                "Ensembl_ID": ensp, "gene_name": ensp2gene[ensp],
                "disease_id": did, "disease_name": dname,
                "source": "CTD", "evidence_type": dev or "curated",
                "score": "", "pmids": pmids,
            }
            ctd_kept += 1
    report["ctd_curated_rows_kept"] = ctd_kept
    print(f"CTD curated: {ctd_kept:,} kept")

    _write_rows(rows, report, stage="1")


def stage2() -> None:
    """CTD aggregate (3.1 GB) — stream-filtered merge into links CSV.

    CTD aggregate carries MULTIPLE rows per (gene, disease) pair (one per
    supporting inference chemical), so per-disease selection must dedup by
    GENE, not by row: a per-disease {ensp: best_score} dict accumulates
    during the single streaming pass (O(unique pairs) memory ≈ 1.6 M
    entries), then each disease is truncated to its top-K genes at flush.
    """
    ensp_set, ensp2gene, sym2ensp = load_protein_index()
    _ensg2ensp_g2e, geneid2ensp = load_gene2ensembl(ensp_set)

    rows = _read_existing()
    before = len(rows)
    seen = kept = 0
    # disease_id -> {ensp: best_inference_score}
    infer_best: dict[str, dict[str, float]] = {}
    # disease_id -> {ensp: (dname, ev_label, score_str, pmids)} for the
    # row that produced the best score (kept for the flush)
    infer_meta: dict[str, dict[str, tuple]] = {}

    def _infer_note(did: str, ensp: str, inf: float,
                    dname: str, pmids: str) -> None:
        best = infer_best.setdefault(did, {})
        if inf <= best.get(ensp, -1.0):
            return
        best[ensp] = inf
        infer_meta.setdefault(did, {})[ensp] = (
            dname, f"inference:{inf:.0f}", round(inf, 2), pmids)

    path = RAW / "CTD_genes_diseases.csv.gz"
    with gzip.open(path, "rt") as fh:
        rdr = csv.reader(filter(lambda ln: not ln.startswith("#"), fh))
        next(rdr, None)
        for parts in rdr:
            seen += 1
            if seen % 10_000_000 == 0:
                print(f"  scanned {seen:,} rows, curated kept {kept:,}, "
                      f"diseases {len(infer_best):,} ...", flush=True)
            if len(parts) < 9:
                continue
            symbol, geneid, dname, did, dev, _chem, iscore, _omim, pmids = parts[:9]
            dev = dev.strip()
            try:
                inf = float(iscore) if iscore else 0.0
            except ValueError:
                inf = 0.0
            ensp = sym2ensp.get(symbol.upper()) or geneid2ensp.get(geneid)
            if ensp is None:
                continue
            if dev:
                # curated/direct evidence rows always pass
                key = (ensp, did, "CTD")
                if key in rows:
                    continue
                rows[key] = {
                    "Ensembl_ID": ensp, "gene_name": ensp2gene[ensp],
                    "disease_id": did, "disease_name": dname,
                    "source": "CTD", "evidence_type": dev,
                    "score": "", "pmids": pmids,
                }
                kept += 1
            elif inf >= CTD_INFER_MIN:
                _infer_note(did, ensp, inf, dname, pmids)

    # flush per-disease top-K unique genes into the link table
    import heapq as _hq
    infer_kept = 0
    for did, best in infer_best.items():
        if len(best) > CTD_INFER_TOP_PER_DISEASE:
            top = _hq.nlargest(CTD_INFER_TOP_PER_DISEASE, best.items(),
                               key=lambda kv: kv[1])
        else:
            top = list(best.items())
        for ensp, _inf in top:
            key = (ensp, did, "CTD")
            if key in rows:
                continue
            dname, ev, sc, pmids = infer_meta[did][ensp]
            rows[key] = {
                "Ensembl_ID": ensp, "gene_name": ensp2gene[ensp],
                "disease_id": did, "disease_name": dname,
                "source": "CTD", "evidence_type": ev,
                "score": sc, "pmids": pmids,
            }
            infer_kept += 1
    kept += infer_kept

    report = {"ctd_aggregate_rows_seen": seen,
              "ctd_aggregate_rows_kept": kept,
              "ctd_inference_kept": infer_kept,
              "ctd_diseases": len(infer_best),
              "rows_before": before}
    print(f"CTD aggregate: {seen:,} scanned → {kept:,} kept "
          f"({infer_kept:,} inference, top {CTD_INFER_TOP_PER_DISEASE} "
          f"genes/disease, score≥{CTD_INFER_MIN:.0f}, "
          f"{len(infer_best):,} diseases)")
    _write_rows(rows, report, stage="2")


def stage3_ot_indirect() -> None:
    """Optional: OT indirect (EFO-propagated) — needs corrupt partitions fixed."""
    import pyarrow.parquet as pq
    import pyarrow.compute as pc
    ensp_set, ensp2gene, sym2ensp = load_protein_index()
    ensg2ensp_g2e, _gid = load_gene2ensembl(ensp_set)
    ensg2ensp = build_ensg2ensp(sym2ensp, ensg2ensp_g2e)
    ot_names = load_ot_disease_names()
    rows = _read_existing()
    before = len(rows)
    kept = bad = 0
    cols = ["diseaseId", "targetId", "associationScore", "evidenceCount"]
    for f in sorted(RAW.glob("ot_indirect/*.parquet")):
        try:
            t = pq.read_table(f, columns=cols)
        except Exception:  # noqa: BLE001
            bad += 1
            continue
        t = t.filter(pc.greater_equal(t["associationScore"], 0.3))
        for d_raw, g, s, ev in zip(t["diseaseId"].to_pylist(),
                                   t["targetId"].to_pylist(),
                                   t["associationScore"].to_pylist(),
                                   t["evidenceCount"].to_pylist()):
            d = _norm_disease_id(d_raw)
            if d is None:
                continue
            ensp = ensg2ensp.get(g)
            if ensp is None or (ensp, d, "OT") in rows:
                continue
            rows[(ensp, d, "OT-indirect")] = {
                "Ensembl_ID": ensp, "gene_name": ensp2gene[ensp],
                "disease_id": d, "disease_name": ot_names.get(d, ""),
                "source": "OT-indirect", "evidence_type": f"indirect;ev={ev}",
                "score": round(s, 4), "pmids": "",
            }
            kept += 1
    report = {"ot_indirect_rows_kept": kept, "unreadable_partitions": bad,
              "rows_before": before}
    print(f"OT indirect: {kept:,} kept, {bad} unreadable partitions")
    _write_rows(rows, report, stage="3")


def filter_measurement_terms() -> None:
    """Drop EFO measurement/trait terms from the links table.

    Over half of the OT disease index (25,334 / 45,896 terms) consists of
    metabolite/phenotype MEASUREMENT traits ("1,2-dilinoleoyl-GPC
    measurement", "blood pressure measurement"), not diseases — OT
    propagates associations to them via EFO.  Their only therapeutic area
    is EFO:0001444 (measurement); terms that ALSO belong to a real
    disease therapeutic area are kept.  Without this filter the "disease"
    count is inflated by ~5,300 phantom terms and enrichment/lookup get
    polluted with trait sets.
    """
    import numpy as np
    import pyarrow.parquet as pq
    f = RAW / "ot_diseases.parquet"
    if not f.exists():
        sys.exit(f"missing {f} — required for the measurement-term filter")
    idx = pq.read_table(f, columns=["id", "therapeuticAreas"]).to_pandas()
    black: set[str] = set()
    for _i, _r in idx.iterrows():
        ta = _r["therapeuticAreas"]
        ta_set = {str(x) for x in np.atleast_1d(ta)} if ta is not None else set()
        if ta_set == {"EFO_0001444"}:
            black.add(_norm_disease_id(_r["id"]) or _r["id"])

    rows = _read_existing()
    before = len(rows)
    dropped = {k for k in rows if k[1] in black}
    for k in dropped:
        del rows[k]
    report = {"measurement_terms_in_index": len(black),
              "rows_dropped": len(dropped), "rows_before": before}
    print(f"measurement-trait filter: dropped {len(dropped):,} rows "
          f"({len({k[1] for k in dropped})} trait terms), "
          f"{len(rows):,} rows remain")
    _write_rows(rows, report, stage="filter")


def gen_gmt(top_n: int = 500) -> None:
    """ disease_ctdot.gmt: top-N scored genes per disease, gene symbols.

    Genes are TAB-separated on each line: the GMT convention (MSigDB
    format) that gseapy actually parses — a single space-separated field
    is read as ONE gene name and enrichment silently finds no overlap.
    """
    import pandas as pd
    src = LINKS_CSV
    if not src.exists():
        # working copy was cleaned up after apply; regenerate from the
        # packaged (gzip-compressed) table instead
        src = DATA / "jnuherb_protein_disease_links.csv.gz"
    if not src.exists():
        sys.exit(f"missing {LINKS_CSV} (and packaged fallback) — run stage 1 first")
    df = pd.read_csv(src, dtype=str, keep_default_na=False)
    df["score_num"] = pd.to_numeric(df["score"], errors="coerce").fillna(0.5)
    # OT rows have no disease_name yet → fall back to disease_id for the term
    df["term"] = df["disease_name"].where(df["disease_name"] != "",
                                          df["disease_id"])
    n_terms = 0
    with open(GMT_OUT, "w", encoding="utf-8") as out:
        for (term, did), grp in df.groupby(["term", "disease_id"]):
            grp = grp.sort_values("score_num", ascending=False).head(top_n)
            genes = grp["gene_name"].dropna().astype(str)
            genes = [g for g in dict.fromkeys(genes) if g]
            if len(genes) < 5:          # tiny sets enrich noisily
                continue
            out.write(f"{term} ({did})\t{did}\t{'\t'.join(genes)}\n")
            n_terms += 1
    print(f"wrote {GMT_OUT}: {n_terms} disease terms (top {top_n} genes each)")


def apply() -> None:
    import gzip
    import shutil
    if not LINKS_CSV.exists():
        sys.exit(f"missing {LINKS_CSV} — run stage 1 first")
    dest = DATA / "jnuherb_protein_disease_links.csv.gz"
    plain = DATA / "jnuherb_protein_disease_links.csv"
    bak = OUT / "jnuherb_protein_disease_links.pre_apply.bak"
    if plain.exists() and not bak.exists():
        shutil.copy2(plain, bak)
    # ship gzip-compressed: 146 MB CSV → ~25 MB, under git's 100 MB limit.
    # pandas/read_csv and build_db transparently decompress by extension.
    with open(LINKS_CSV, "rb") as fin, \
            gzip.open(dest, "wb", compresslevel=6) as fout:
        shutil.copyfileobj(fin, fout)
    if plain.exists():
        plain.unlink()
    print(f"installed {dest}")
    if GMT_OUT.exists():
        print(f"disease_ctdot.gmt already generated at {GMT_OUT}")
    print("Next: python -m jnuherb.build_db")


def _read_existing() -> dict[tuple[str, str, str], dict]:
    rows: dict[tuple[str, str, str], dict] = {}
    if LINKS_CSV.exists():
        with open(LINKS_CSV, encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                rows[(r["Ensembl_ID"], r["disease_id"], r["source"])] = r
    return rows


def _write_rows(rows: dict, report: dict, stage: str) -> None:
    fields = ["Ensembl_ID", "gene_name", "disease_id", "disease_name",
              "source", "evidence_type", "score", "pmids"]
    with open(LINKS_CSV, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in rows.values():
            w.writerow(r)
    report["total_rows"] = len(rows)
    prev = {}
    if REPORT.exists():
        prev = json.load(open(REPORT, encoding="utf-8"))
    prev[f"stage{stage}"] = report
    json.dump(prev, open(REPORT, "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print(f"total links: {len(rows):,} → {LINKS_CSV}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("stage",
                    choices=["1", "2", "3", "filter", "gmt", "apply"])
    ap.add_argument("--top-n", type=int, default=500, help="genes per disease in GMT")
    args = ap.parse_args()
    if args.stage == "1":
        stage1()
    elif args.stage == "2":
        stage2()
    elif args.stage == "3":
        stage3_ot_indirect()
    elif args.stage == "filter":
        filter_measurement_terms()
    elif args.stage == "gmt":
        gen_gmt(args.top_n)
    else:
        apply()
    return 0


if __name__ == "__main__":
    sys.exit(main())
