#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${ROOT}/.." && pwd)"
DATA_DIR="${ROOT}/data"
RESULTS_DIR="${ROOT}/results"
THREADS="${THREADS:-1}"

BW_A="${DATA_DIR}/demo_sample_A.bw"
BW_B="${DATA_DIR}/demo_sample_B.bw"
BED_PATH="${DATA_DIR}/demo_regions.bed"

if [[ ! -f "${BW_A}" || ! -f "${BW_B}" || ! -f "${BED_PATH}" ]]; then
  echo "Generating synthetic demo tracks in ${DATA_DIR}" >&2
  python "${ROOT}/make_demo_tracks.py" --output "${DATA_DIR}"
fi

mkdir -p "${RESULTS_DIR}"

python "${PROJECT_ROOT}/peak_shape.py" \
  --bigwig-a "${BW_A}" \
  --bigwig-b "${BW_B}" \
  --bed "${BED_PATH}" \
  --core 400 \
  --flank 800 2000 \
  --out "${RESULTS_DIR}" \
  --threads "${THREADS}"

echo "Results written to ${RESULTS_DIR}"
