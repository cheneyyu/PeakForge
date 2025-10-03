#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${ROOT}/data"
MANIFEST="${DATA_DIR}/encode_manifest.tsv"
mkdir -p "${DATA_DIR}"

: "${DRY_RUN:=1}"

cat <<'BANNER'
This helper fetches a matched 2 vs 2 ENCODE ChIP-seq dataset (hg38):
  - GM12878 H3K27ac (ENCSR000AKP)
  - K562 H3K27ac (ENCSR000AKO)

By default the script runs in DRY-RUN mode and only prints the curl
commands.  Set DRY_RUN=0 to download the files.
BANNER

python3 - <<'PY' "${DATA_DIR}" "${MANIFEST}" "${DRY_RUN}"
import json
import sys
import urllib.parse
import urllib.request
from pathlib import Path

BASE = "https://www.encodeproject.org"
HEADERS = {
    "User-Agent": "PeakForge-fetcher/1.0 (+https://github.com/)",
    "Accept": "application/json",
}

def fetch_json(url: str) -> dict:
    req = urllib.request.Request(url, headers=HEADERS)
    # Respect proxies from the environment
    opener = urllib.request.build_opener()
    with opener.open(req) as resp:
        if resp.status != 200:
            raise RuntimeError(f"Failed to fetch {url} (HTTP {resp.status})")
        return json.load(resp)

def iter_files(experiment: str) -> list[dict]:
    params = {
        "type": "File",
        "dataset": experiment,
        "output_type": "alignments",
        "assembly": "GRCh38",
        "file_format": "bam",
        "status": "released",
        "limit": "all",
        "format": "json",
        "field": [
            "accession",
            "href",
            "md5sum",
            "file_format",
            "file_type",
            "replicate.title",
            "replicate.biological_replicate_number",
            "replicate.technical_replicate_number",
            "preferred_default",
        ],
    }
    query = urllib.parse.urlencode(params, doseq=True)
    url = f"{BASE}/search/?{query}"
    data = fetch_json(url)
    return data["@graph"]

def select_replicates(records: list[dict], count: int = 2) -> list[dict]:
    sorted_records = sorted(
        records,
        key=lambda rec: (
            rec.get("replicate", {}).get("biological_replicate_number", 99),
            rec.get("replicate", {}).get("technical_replicate_number", 99),
            not rec.get("preferred_default", False),
            rec.get("accession"),
        ),
    )
    selected = []
    seen = set()
    for rec in sorted_records:
        rep = rec.get("replicate") or {}
        bio = rep.get("biological_replicate_number")
        if bio is None or bio in seen:
            continue
        selected.append(rec)
        seen.add(bio)
        if len(selected) == count:
            break
    if len(selected) < count:
        raise RuntimeError(f"Expected >= {count} replicates but found {len(selected)}")
    return selected

def main(data_dir: Path, manifest_path: Path, dry_run: bool) -> None:
    experiments = [
        ("GM12878", "ENCSR000AKP"),
        ("K562", "ENCSR000AKO"),
    ]
    rows = ["sample\tcondition\taccession\tdestination\tmd5sum\turl"]
    for condition, experiment in experiments:
        files = select_replicates(iter_files(experiment))
        for rec in files:
            rep = rec.get("replicate") or {}
            bio = rep.get("biological_replicate_number")
            accession = rec["accession"]
            ext = rec.get("file_format", "bam")
            dest = data_dir / f"{condition}_rep{bio}.{ext}"
            download_url = f"{BASE}/files/{accession}/@@download/{accession}.{ext}"
            rows.append(
                "\t".join(
                    [
                        f"{condition}_rep{bio}",
                        condition,
                        accession,
                        str(dest),
                        rec.get("md5sum", ""),
                        download_url,
                    ]
                )
            )
            cmd = [
                "curl",
                "-L",
                "-H",
                "Accept: application/json",
                "-o",
                str(dest),
                download_url,
            ]
            print("[command]", " ".join(cmd))
            if not dry_run:
                data_dir.mkdir(parents=True, exist_ok=True)
                import subprocess

                subprocess.run(cmd, check=True)
    manifest_path.write_text("\n".join(rows) + "\n")
    print(f"Manifest written to {manifest_path}")

if __name__ == "__main__":
    data_dir = Path(sys.argv[1])
    manifest_path = Path(sys.argv[2])
    dry_run = bool(int(sys.argv[3]))
    main(data_dir, manifest_path, dry_run)
PY

if [[ "${DRY_RUN}" == "1" ]]; then
  echo "(dry-run) No files were downloaded. Re-run with DRY_RUN=0 to fetch data."
fi

python3 - <<'PY' "${DATA_DIR}"
from pathlib import Path
metadata = Path("example/data/metadata.tsv")
if metadata.exists():
    print(f"Metadata sheet already present: {metadata}")
else:
    print("WARNING: metadata.tsv is missing; create it before running the pipeline.")
PY
