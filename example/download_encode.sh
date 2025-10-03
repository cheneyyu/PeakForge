#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${ROOT}/data"
MANIFEST="${DATA_DIR}/encode_manifest.tsv"
mkdir -p "${DATA_DIR}"

: "${DRY_RUN:=1}"

cat <<'BANNER'
This helper fetches a matched 2 vs 2 ENCODE MYC ChIP-seq dataset (HepG2, GRCh38):
  - Experiment A (replicates 1 & 2): ENCFF975ETI, ENCFF380OWL
  - Experiment B (replicates 1 & 2): ENCFF315AUW, ENCFF987GJQ

By default the script runs in DRY-RUN mode and only prints the curl
commands.  Set DRY_RUN=0 to download the files (~11 GB total).
BANNER

FILES=$'HepG2_MYCa_rep1\tHepG2_MYCa\tENCFF975ETI\thttps://www.encodeproject.org/files/ENCFF975ETI/@@download/ENCFF975ETI.bam\n'
FILES+=$'HepG2_MYCa_rep2\tHepG2_MYCa\tENCFF380OWL\thttps://www.encodeproject.org/files/ENCFF380OWL/@@download/ENCFF380OWL.bam\n'
FILES+=$'HepG2_MYCb_rep1\tHepG2_MYCb\tENCFF315AUW\thttps://www.encodeproject.org/files/ENCFF315AUW/@@download/ENCFF315AUW.bam\n'
FILES+=$'HepG2_MYCb_rep2\tHepG2_MYCb\tENCFF987GJQ\thttps://www.encodeproject.org/files/ENCFF987GJQ/@@download/ENCFF987GJQ.bam\n'

rows=("sample\tcondition\taccession\tdestination\tmd5sum\turl")
while IFS=$'\t' read -r sample condition accession url; do
  [[ -z "${sample}" ]] && continue
  dest="${DATA_DIR}/${sample}.bam"
  rows+=("${sample}\t${condition}\t${accession}\t${dest}\t\t${url}")
  cmd=(
    curl -L
    -H "Accept: application/json"
    -o "${dest}"
    "${url}"
  )
  printf '[command] %s\n' "${cmd[*]}"
  if [[ "${DRY_RUN}" != "1" ]]; then
    mkdir -p "${DATA_DIR}"
    "${cmd[@]}"
  fi
done <<< "${FILES}"

printf '%s\n' "${rows[@]}" > "${MANIFEST}"
printf 'Manifest written to %s\n' "${MANIFEST}"

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
