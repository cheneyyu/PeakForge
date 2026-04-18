#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${ROOT}/.." && pwd)"
PYTHON_BIN="${PROJECT_ROOT}/.venv/bin/python"
METADATA="${ROOT}/data/metadata_1v1.tsv"
RESULTS_DIR="${RESULTS_DIR:-${ROOT}/results_1v1}"
THREADS="${THREADS:-16}"
USE_CONTROLS="${USE_CONTROLS:-1}"

if [[ ! -x "${PYTHON_BIN}" ]]; then
  PYTHON_BIN="$(command -v python3)"
fi

if [[ ! -f "${METADATA}" ]]; then
  echo "Metadata sheet not found: ${METADATA}" >&2
  exit 1
fi

mkdir -p "${RESULTS_DIR}" "${RESULTS_DIR}/peaks"

RUN_METADATA="${METADATA}"
if [[ "${USE_CONTROLS}" != "1" ]]; then
  RUN_METADATA="${RESULTS_DIR}/metadata.no_controls.tsv"
  awk 'BEGIN{FS=OFS="\t"} NR==1{for(i=1;i<=NF;i++) idx[$i]=i; print; next} {if(("control_bam" in idx)) $idx["control_bam"]=""; if(("input_bam" in idx)) $idx["input_bam"]=""; print}' "${METADATA}" > "${RUN_METADATA}"
fi

echo "[peakforge] Running 1v1 example -> ${RESULTS_DIR}" \
  && (
    cd "${PROJECT_ROOT}"
    "${PYTHON_BIN}" "${PROJECT_ROOT}/chipdiff.py" \
      tsvmode "${RUN_METADATA}" \
      --output-dir "${RESULTS_DIR}" \
      --peak-dir "${RESULTS_DIR}/peaks" \
      --min-overlap 1 \
      --macs2-genome hs \
      --threads "${THREADS}"
  )

echo "Results written to ${RESULTS_DIR}" \
  && echo "Inspect ${RESULTS_DIR}/differential_results.tsv for MARS outputs."
