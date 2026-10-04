"""Build Ensembl protein ID -> UniProt ID mapping for the JNUHERB database.

The website's protein table only stores Ensembl IDs (ENSP...), but predicted
protein structures are fetched from AlphaFold DB, which is keyed by UniProt
accessions. This script uses the UniProt ID Mapping API (async job,
Ensembl_Protein -> UniProtKB, 2000 IDs per job) and writes a two-column CSV
that ships with the package so the web app never needs external calls for
the mapping at runtime.

Usage:
    python scripts/build_uniprot_mapping.py

Output:
    src/jnuherb/data/ensembl_uniprot.csv   (Ensembl_ID,UniProt_ID)
"""
from pathlib import Path
import csv
import time
import urllib.parse
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
PROTEINS_CSV = ROOT / "src" / "jnuherb" / "data" / "jnuherb_proteins.csv"
OUT_CSV = ROOT / "src" / "jnuherb" / "data" / "ensembl_uniprot.csv"

IDMAP_RUN = "https://rest.uniprot.org/idmapping/run"
IDMAP_STATUS = "https://rest.uniprot.org/idmapping/status/{job}"
IDMAP_RESULT = "https://rest.uniprot.org/idmapping/uniprotkb/results/stream/{job}?format=tsv"

BATCH = 1000
POLL_INTERVAL = 3
POLL_TIMEOUT = 600


def read_ensembl_ids() -> list[str]:
    with open(PROTEINS_CSV, newline="") as f:
        return [row["Ensembl_ID"] for row in csv.DictReader(f)]


def http_get(url: str, timeout: int = 120) -> str:
    with urllib.request.urlopen(url, timeout=timeout) as resp:
        return resp.read().decode()


def submit_job(ids: list[str]) -> str:
    body = urllib.parse.urlencode({
        "from": "Ensembl_Protein",
        "to": "UniProtKB",
        "ids": ",".join(ids),
    }).encode()
    req = urllib.request.Request(
        IDMAP_RUN, data=body,
        headers={"Content-Type": "application/x-www-form-urlencoded"},
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=120) as resp:
        import json
        return json.loads(resp.read().decode())["jobId"]


def wait_job(job_id: str) -> None:
    deadline = time.time() + POLL_TIMEOUT
    while time.time() < deadline:
        import json
        status = json.loads(http_get(IDMAP_STATUS.format(job=job_id), timeout=30))
        # Small jobs often finish before the first poll: the status endpoint
        # then 302-redirects (urllib follows it) straight to the results.
        if "results" in status:
            return
        if status.get("jobStatus") == "FINISHED":
            return
        if "jobStatus" in status and status["jobStatus"] not in ("RUNNING", "NEW"):
            raise SystemExit(f"[mapping] job {job_id}: {status}")
        time.sleep(POLL_INTERVAL)
    raise SystemExit(f"[mapping] job {job_id} timed out")


def fetch_results(job_id: str) -> dict[str, str]:
    """Return {Ensembl_ID: UniProt accession}, reviewed entries preferred."""
    tsv = http_get(IDMAP_RESULT.format(job=job_id), timeout=300)
    best: dict[str, tuple[int, str]] = {}
    rows = csv.DictReader(tsv.splitlines(), delimiter="\t")
    for row in rows:
        ens, entry = row.get("From", ""), row.get("Entry", "")
        if not ens or not entry:
            continue
        reviewed = 0 if row.get("Reviewed") == "reviewed" else 1
        if ens not in best or reviewed < best[ens][0]:
            best[ens] = (reviewed, entry)
    return {k: v[1] for k, v in best.items()}


def main() -> None:
    ids = read_ensembl_ids()
    print(f"[mapping] {len(ids)} Ensembl IDs from {PROTEINS_CSV.name}")

    mapping: dict[str, str] = {}
    for start in range(0, len(ids), BATCH):
        chunk = ids[start:start + BATCH]
        for attempt in range(3):
            try:
                job = submit_job(chunk)
                wait_job(job)
                got = fetch_results(job)
                mapping.update(got)
                break
            except Exception as e:  # noqa: BLE001 - retry transient failures
                print(f"[mapping] batch {start} attempt {attempt + 1} failed: {e}")
                time.sleep(5 * (attempt + 1))
        else:
            raise SystemExit(f"[mapping] batch starting at {start} failed 3 times")
        print(f"[mapping] {start + len(chunk)}/{len(ids)} queried, "
              f"{len(mapping)} mapped so far")

    with open(OUT_CSV, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Ensembl_ID", "UniProt_ID"])
        for ens_id in ids:  # keep database order, stable diffs
            if ens_id in mapping:
                w.writerow([ens_id, mapping[ens_id]])
    print(f"[mapping] wrote {len(mapping)} mappings -> {OUT_CSV}")


if __name__ == "__main__":
    main()
