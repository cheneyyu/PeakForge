#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${ROOT}/.." && pwd)"
DATA_DIR="${ROOT}/data"
DATA2_DIR="${ROOT}/data2"
MANIFEST="${DATA_DIR}/encode_manifest.tsv"
MANIFEST_3V3="${DATA_DIR}/encode_manifest_3v3.tsv"
METADATA="${DATA_DIR}/metadata.tsv"
METADATA_1V1="${DATA_DIR}/metadata_1v1.tsv"
METADATA_3V3="${DATA_DIR}/metadata_3v3.tsv"
mkdir -p "${DATA_DIR}" "${DATA2_DIR}"

: "${DRY_RUN:=1}"
: "${INCLUDE_THIRD_REPLICATES:=0}"

cat <<'BANNER'
This helper fetches a matched ENCODE MYC ChIP-seq benchmark (hg38)
together with the corresponding common-input controls:
  - K562 MYC (ENCFF975ETI, ENCFF380OWL)
  - HepG2 MYC (ENCFF315AUW, ENCFF987GJQ)
  - K562 input (ENCFF500EAO)
  - HepG2 input (ENCFF062SCH)

Set DRY_RUN=0 to download the files.

Set INCLUDE_THIRD_REPLICATES=1 to also fetch:
  - K562 MYC rep3 (ENCFF953SVM)
  - HepG2 MYC rep3 (ENCFF439SHL)

Those two extra BAMs are only needed for the optional 3-fold held-out validation.
BANNER

python3 - <<'PY' "${DATA_DIR}" "${DATA2_DIR}" "${MANIFEST}" "${MANIFEST_3V3}" "${METADATA}" "${METADATA_1V1}" "${METADATA_3V3}" "${PROJECT_ROOT}" "${DRY_RUN}" "${INCLUDE_THIRD_REPLICATES}"
import shutil
import subprocess
import sys
from pathlib import Path

CORE_FILES = [
    {
        "sample": "K562_rep1",
        "condition": "K562",
        "role": "chip",
        "accession": "ENCFF975ETI",
        "filename": "K562_rep1.bam",
        "subdir": "data",
        "url": "https://www.encodeproject.org/files/ENCFF975ETI/@@download/ENCFF975ETI.bam",
    },
    {
        "sample": "K562_rep2",
        "condition": "K562",
        "role": "chip",
        "accession": "ENCFF380OWL",
        "filename": "K562_rep2.bam",
        "subdir": "data",
        "url": "https://www.encodeproject.org/files/ENCFF380OWL/@@download/ENCFF380OWL.bam",
    },
    {
        "sample": "HepG2_rep1",
        "condition": "HepG2",
        "role": "chip",
        "accession": "ENCFF315AUW",
        "filename": "HepG2_rep1.bam",
        "subdir": "data",
        "url": "https://www.encodeproject.org/files/ENCFF315AUW/@@download/ENCFF315AUW.bam",
    },
    {
        "sample": "HepG2_rep2",
        "condition": "HepG2",
        "role": "chip",
        "accession": "ENCFF987GJQ",
        "filename": "HepG2_rep2.bam",
        "subdir": "data",
        "url": "https://www.encodeproject.org/files/ENCFF987GJQ/@@download/ENCFF987GJQ.bam",
    },
    {
        "sample": "K562_input",
        "condition": "K562",
        "role": "input",
        "accession": "ENCFF500EAO",
        "filename": "K562_input.bam",
        "subdir": "data",
        "url": "https://www.encodeproject.org/files/ENCFF500EAO/@@download/ENCFF500EAO.bam",
    },
    {
        "sample": "HepG2_input",
        "condition": "HepG2",
        "role": "input",
        "accession": "ENCFF062SCH",
        "filename": "HepG2_input.bam",
        "subdir": "data",
        "url": "https://www.encodeproject.org/files/ENCFF062SCH/@@download/ENCFF062SCH.bam",
    },
]

EXTRA_FILES = [
    {
        "sample": "K562_rep3",
        "condition": "K562",
        "role": "chip",
        "accession": "ENCFF953SVM",
        "filename": "ENCFF953SVM.bam",
        "subdir": "data2",
        "url": "https://www.encodeproject.org/files/ENCFF953SVM/@@download/ENCFF953SVM.bam",
    },
    {
        "sample": "HepG2_rep3",
        "condition": "HepG2",
        "role": "chip",
        "accession": "ENCFF439SHL",
        "filename": "ENCFF439SHL.bam",
        "subdir": "data2",
        "url": "https://www.encodeproject.org/files/ENCFF439SHL/@@download/ENCFF439SHL.bam",
    },
]


def local_path(entry: dict[str, str], *, data_dir: Path, data2_dir: Path) -> Path:
    base = data_dir if entry["subdir"] == "data" else data2_dir
    return base / entry["filename"]


def rel_path(entry: dict[str, str]) -> str:
    return str(Path("example") / entry["subdir"] / entry["filename"])


def write_manifest(path: Path, entries: list[dict[str, str]], *, data_dir: Path, data2_dir: Path) -> None:
    rows = ["sample\tcondition\trole\taccession\tdestination\tmd5sum\turl"]
    for entry in entries:
        rows.append(
            "\t".join(
                [
                    entry["sample"],
                    entry["condition"],
                    entry["role"],
                    entry["accession"],
                    str(local_path(entry, data_dir=data_dir, data2_dir=data2_dir)),
                    "",
                    entry["url"],
                ]
            )
        )
    path.write_text("\n".join(rows) + "\n")
    print(f"Manifest written to {path}")


def write_metadata_tables(
    metadata_path: Path,
    metadata_1v1_path: Path,
    metadata_3v3_path: Path,
    *,
    include_third_replicates: bool,
) -> None:
    controls = {
        entry["condition"]: rel_path(entry)
        for entry in CORE_FILES
        if entry["role"] == "input"
    }

    core_chip = [entry for entry in CORE_FILES if entry["role"] == "chip"]

    rows_main = ["sample\tcondition\tbam\tcontrol_bam"]
    for entry in core_chip:
        rows_main.append(
            "\t".join(
                [
                    entry["sample"],
                    entry["condition"],
                    rel_path(entry),
                    controls[entry["condition"]],
                ]
            )
        )
    metadata_path.write_text("\n".join(rows_main) + "\n")
    print(f"Metadata sheet written to {metadata_path}")

    rows_1v1 = [
        "sample\tcondition\tbam\tcontrol_bam",
        "\t".join(["K562_rep1", "K562", "example/data/K562_rep1.bam", controls["K562"]]),
        "\t".join(["HepG2_rep1", "HepG2", "example/data/HepG2_rep1.bam", controls["HepG2"]]),
    ]
    metadata_1v1_path.write_text("\n".join(rows_1v1) + "\n")
    print(f"Metadata sheet written to {metadata_1v1_path}")

    if include_third_replicates:
        rows_3v3 = rows_main[:]
        for entry in EXTRA_FILES:
            rows_3v3.append(
                "\t".join(
                    [
                        entry["sample"],
                        entry["condition"],
                        rel_path(entry),
                        controls[entry["condition"]],
                    ]
                )
            )
        metadata_3v3_path.write_text("\n".join(rows_3v3) + "\n")
        print(f"Metadata sheet written to {metadata_3v3_path}")


def download_entries(entries: list[dict[str, str]], *, data_dir: Path, data2_dir: Path, dry_run: bool) -> list[Path]:
    downloaded: list[Path] = []
    for entry in entries:
        dest = local_path(entry, data_dir=data_dir, data2_dir=data2_dir)
        dest.parent.mkdir(parents=True, exist_ok=True)
        if dest.exists() and not dry_run:
            print(f"[skip] {dest} already exists; reusing existing file.")
            downloaded.append(dest)
            continue
        cmd = ["curl", "-L", "-H", "Accept: application/json", "-o", str(dest), entry["url"]]
        print("[command]", " ".join(cmd))
        if not dry_run:
            subprocess.run(cmd, check=True)
            downloaded.append(dest)
    return downloaded


def index_bams(paths: list[Path]) -> None:
    if not paths:
        return

    samtools = shutil.which("samtools")
    if samtools is None:
        print("WARNING: samtools not found on PATH; BAM indices (.bai) were not created.")
        return

    for path in paths:
        print(f"[command] {samtools} index {path}")
        subprocess.run([samtools, "index", str(path)], check=True)

    print("BAM indexing complete.")


def main(
    data_dir: Path,
    data2_dir: Path,
    manifest_path: Path,
    manifest_3v3_path: Path,
    metadata_path: Path,
    metadata_1v1_path: Path,
    metadata_3v3_path: Path,
    dry_run: bool,
    include_third_replicates: bool,
) -> None:
    all_entries = CORE_FILES + EXTRA_FILES
    write_manifest(manifest_path, CORE_FILES, data_dir=data_dir, data2_dir=data2_dir)
    write_manifest(manifest_3v3_path, all_entries, data_dir=data_dir, data2_dir=data2_dir)

    downloaded = download_entries(CORE_FILES, data_dir=data_dir, data2_dir=data2_dir, dry_run=dry_run)
    extra_downloaded: list[Path] = []
    if include_third_replicates:
        extra_downloaded = download_entries(EXTRA_FILES, data_dir=data_dir, data2_dir=data2_dir, dry_run=dry_run)

    if dry_run:
        print("(dry-run) Skipping metadata generation and BAM indexing.")
        return

    write_metadata_tables(
        metadata_path,
        metadata_1v1_path,
        metadata_3v3_path,
        include_third_replicates=include_third_replicates,
    )
    index_bams(downloaded + extra_downloaded)


if __name__ == "__main__":
    data_dir = Path(sys.argv[1])
    data2_dir = Path(sys.argv[2])
    manifest_path = Path(sys.argv[3])
    manifest_3v3_path = Path(sys.argv[4])
    metadata_path = Path(sys.argv[5])
    metadata_1v1_path = Path(sys.argv[6])
    metadata_3v3_path = Path(sys.argv[7])
    _project_root = Path(sys.argv[8])
    dry_run = bool(int(sys.argv[9]))
    include_third_replicates = bool(int(sys.argv[10]))
    main(
        data_dir,
        data2_dir,
        manifest_path,
        manifest_3v3_path,
        metadata_path,
        metadata_1v1_path,
        metadata_3v3_path,
        dry_run,
        include_third_replicates,
    )
PY

if [[ "${DRY_RUN}" == "1" ]]; then
  echo "(dry-run) No files were downloaded. Re-run with DRY_RUN=0 to fetch data."
elif [[ ! -f "${METADATA}" || ! -f "${METADATA_1V1}" ]]; then
  echo "WARNING: Failed to create metadata sheets under ${DATA_DIR}" >&2
elif [[ "${INCLUDE_THIRD_REPLICATES}" == "1" && ! -f "${METADATA_3V3}" ]]; then
  echo "WARNING: Failed to create 3v3 metadata sheet at ${METADATA_3V3}" >&2
fi
