#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${ROOT}/data"
MANIFEST="${DATA_DIR}/encode_manifest.tsv"
mkdir -p "${DATA_DIR}"

: "${DRY_RUN:=1}"

cat <<'BANNER'
This helper fetches a matched 2 vs 2 ENCODE MYC ChIP-seq dataset (hg38):
  - K562 MYC (ENCFF975ETI, ENCFF380OWL)
  - HepG2 MYC (ENCFF315AUW, ENCFF987GJQ)

By default the script runs in DRY-RUN mode and only prints the curl
commands.  Set DRY_RUN=0 to download the files.
BANNER

python3 - <<'PY' "${DATA_DIR}" "${MANIFEST}" "${DRY_RUN}"
import subprocess
import sys
from pathlib import Path

FILES = [
    {
        "sample": "K562_rep1",
        "condition": "K562",
        "accession": "ENCFF975ETI",
        "url": "https://www.encodeproject.org/files/ENCFF975ETI/@@download/ENCFF975ETI.bam",
    },
    {
        "sample": "K562_rep2",
        "condition": "K562",
        "accession": "ENCFF380OWL",
        "url": "https://www.encodeproject.org/files/ENCFF380OWL/@@download/ENCFF380OWL.bam",
    },
    {
        "sample": "HepG2_rep1",
        "condition": "HepG2",
        "accession": "ENCFF315AUW",
        "url": "https://www.encodeproject.org/files/ENCFF315AUW/@@download/ENCFF315AUW.bam",
    },
    {
        "sample": "HepG2_rep2",
        "condition": "HepG2",
        "accession": "ENCFF987GJQ",
        "url": "https://www.encodeproject.org/files/ENCFF987GJQ/@@download/ENCFF987GJQ.bam",
    },
]


def main(data_dir: Path, manifest_path: Path, dry_run: bool) -> None:
    data_dir.mkdir(parents=True, exist_ok=True)
    rows = ["sample\tcondition\taccession\tdestination\tmd5sum\turl"]
    for entry in FILES:
        dest = data_dir / f"{entry['sample']}.bam"
        rows.append(
            "\t".join(
                [
                    entry["sample"],
                    entry["condition"],
                    entry["accession"],
                    str(dest),
                    "",
                    entry["url"],
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
            entry["url"],
        ]
        print("[command]", " ".join(cmd))
        if not dry_run:
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
