#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${ROOT}/.." && pwd)"
METADATA="${ROOT}/data/metadata_1v1.tsv"
RESULTS_DIR="${RESULTS_DIR:-${ROOT}/results_1v1}"
THREADS="${THREADS:-4}"

if [[ ! -f "${METADATA}" ]]; then
  echo "Metadata sheet not found: ${METADATA}" >&2
  exit 1
fi

mkdir -p "${RESULTS_DIR}" "${RESULTS_DIR}/peaks"

echo "[chipdiff] Running 1v1 example -> ${RESULTS_DIR}" \
  && python "${PROJECT_ROOT}/chipdiff.py" \
    --metadata "${METADATA}" \
    --output-dir "${RESULTS_DIR}" \
    --peak-dir "${RESULTS_DIR}/peaks" \
    --min-overlap 1 \
    --macs2-genome hs \
    --threads "${THREADS}"

echo "Results written to ${RESULTS_DIR}" \
  && echo "Inspect ${RESULTS_DIR}/differential_results.tsv for MARS outputs."
