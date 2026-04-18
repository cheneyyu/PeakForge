#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${ROOT}/.." && pwd)"
PYTHON_BIN="${PROJECT_ROOT}/.venv/bin/python"
DATA_DIR="${ROOT}/data"
METADATA="${DATA_DIR}/metadata.tsv"
RESULTS_ROOT="${RESULTS_ROOT:-${ROOT}/results}"
THREADS="${THREADS:-16}"
USE_CONTROLS="${USE_CONTROLS:-1}"

if [[ ! -x "${PYTHON_BIN}" ]]; then
  PYTHON_BIN="$(command -v python3)"
fi

if [[ ! -f "${METADATA}" ]]; then
  echo "Metadata sheet not found: ${METADATA}" >&2
  exit 1
fi

mkdir -p "${RESULTS_ROOT}"

RUN_METADATA="${METADATA}"
if [[ "${USE_CONTROLS}" != "1" ]]; then
  RUN_METADATA="${RESULTS_ROOT}/metadata.no_controls.tsv"
  awk 'BEGIN{FS=OFS="\t"} NR==1{for(i=1;i<=NF;i++) idx[$i]=i; print; next} {if(("control_bam" in idx)) $idx["control_bam"]=""; if(("input_bam" in idx)) $idx["input_bam"]=""; print}' "${METADATA}" > "${RUN_METADATA}"
fi

# Parse metadata.tsv to build sample -> bam and condition mappings.
declare -A SAMPLE_BAMS
declare -A CONDITION_SAMPLES
declare -A CONDITION_SEEN
declare -a CONDITION_ORDER

while IFS=$'\t' read -r sample condition bam _rest; do
  if [[ -z "${sample}" || "${sample}" == sample ]]; then
    continue
  fi
  SAMPLE_BAMS["${sample}"]="${bam}"
  CONDITION_SAMPLES["${condition}"]+=" ${sample}"
  if [[ -z "${CONDITION_SEEN[${condition}]:-}" ]]; then
    CONDITION_ORDER+=("${condition}")
    CONDITION_SEEN["${condition}"]=1
  fi
done < "${RUN_METADATA}"

if [[ "${#CONDITION_ORDER[@]}" -ne 2 ]]; then
  echo "Expected exactly 2 conditions but found ${#CONDITION_ORDER[@]}" >&2
  exit 1
fi

COND_A="${CONDITION_ORDER[0]}"
COND_B="${CONDITION_ORDER[1]}"
read -ra REPS_A <<<"${CONDITION_SAMPLES[${COND_A}]}"
read -ra REPS_B <<<"${CONDITION_SAMPLES[${COND_B}]}"

if [[ "${#REPS_A[@]}" -lt 2 || "${#REPS_B[@]}" -lt 2 ]]; then
  echo "Need at least two replicates per condition for 2v2 example" >&2
  exit 1
fi

run_chipdiff() {
  local metadata_path="$1"
  local output_dir="$2"
  shift 2
  mkdir -p "${output_dir}"
  echo "[peakforge] ${metadata_path} -> ${output_dir}"
  (
    cd "${PROJECT_ROOT}"
    "${PYTHON_BIN}" "${PROJECT_ROOT}/chipdiff.py" \
      tsvmode "${metadata_path}" \
      --output-dir "${output_dir}" \
      --peak-dir "${output_dir}/peaks" \
      --min-overlap 2 \
      --macs2-genome hs \
      --threads "${THREADS}" \
      "$@"
  )
}

MAIN_DIR="${RESULTS_ROOT}/2v2"
run_chipdiff "${RUN_METADATA}" "${MAIN_DIR}"

echo "All analyses completed. Review ${MAIN_DIR} for results."
echo "Consensus peaks saved to ${MAIN_DIR}/consensus_peaks.bed (pass via --consensus-peaks for follow-up runs)."
