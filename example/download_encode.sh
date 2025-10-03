#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${ROOT}/.." && pwd)"
DATA_DIR="${ROOT}/data"
MANIFEST="${DATA_DIR}/encode_manifest.tsv"
METADATA="${DATA_DIR}/metadata.tsv"
mkdir -p "${DATA_DIR}"

: "${DRY_RUN:=1}"

cat <<'BANNER'
This helper fetches a matched 2 vs 2 ENCODE MYC ChIP-seq dataset (hg38):
  - K562 MYC (ENCFF975ETI, ENCFF380OWL)
  - HepG2 MYC (ENCFF315AUW, ENCFF987GJQ)

By default the script runs in DRY-RUN mode and only prints the curl
commands.  Set DRY_RUN=0 to download the files.
BANNER

python3 - <<'PY' "${DATA_DIR}" "${MANIFEST}" "${METADATA}" "${PROJECT_ROOT}" "${DRY_RUN}"
import subprocess
import sys
import shutil
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


def write_metadata(metadata_path: Path, entries: list[dict[str, str]], *, base_dir: Path, project_root: Path) -> None:
    rows = ["sample\tcondition\tbam"]
    try:
        rel_base = base_dir.relative_to(project_root)
    except ValueError:
        rel_base = Path("example") / "data"
    for entry in entries:
        bam_path = base_dir / f"{entry['sample']}.bam"
        try:
            bam_rel = bam_path.relative_to(project_root)
        except ValueError:
            bam_rel = rel_base / f"{entry['sample']}.bam"
        rows.append("\t".join([entry["sample"], entry["condition"], str(bam_rel)]))
    metadata_path.write_text("\n".join(rows) + "\n")
    print(f"Metadata sheet written to {metadata_path}")


def main(data_dir: Path, manifest_path: Path, metadata_path: Path, project_root: Path, dry_run: bool) -> None:
    data_dir.mkdir(parents=True, exist_ok=True)
    rows = ["sample\tcondition\taccession\tdestination\tmd5sum\turl"]
    downloaded: list[Path] = []
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
        if dest.exists() and not dry_run:
            print(f"[skip] {dest} already exists; reusing existing file.")
            downloaded.append(dest)
            continue
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
            downloaded.append(dest)

    manifest_path.write_text("\n".join(rows) + "\n")
    print(f"Manifest written to {manifest_path}")

    if dry_run:
        print("(dry-run) Skipping metadata generation and BAM indexing.")
        return

    write_metadata(metadata_path, FILES, base_dir=data_dir, project_root=project_root)

    samtools = shutil.which("samtools")
    if samtools is None:
        print("WARNING: samtools not found on PATH; BAM indices (.bai) were not created.")
        return

    for path in downloaded:
        print(f"[command] {samtools} index {path}")
        subprocess.run([samtools, "index", str(path)], check=True)

    print("BAM indexing complete.")


if __name__ == "__main__":
    data_dir = Path(sys.argv[1])
    manifest_path = Path(sys.argv[2])
    metadata_path = Path(sys.argv[3])
    project_root = Path(sys.argv[4])
    dry_run = bool(int(sys.argv[5]))
    main(data_dir, manifest_path, metadata_path, project_root, dry_run)
PY

if [[ "${DRY_RUN}" == "1" ]]; then
  echo "(dry-run) No files were downloaded. Re-run with DRY_RUN=0 to fetch data."
elif [[ ! -f "${METADATA}" ]]; then
  echo "WARNING: Failed to create metadata sheet at ${METADATA}" >&2
fi
