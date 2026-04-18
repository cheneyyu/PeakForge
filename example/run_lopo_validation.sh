#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${ROOT}/.." && pwd)"
ENV_BIN="${PROJECT_ROOT}/.venv/bin"
PYTHON_BIN="${ENV_BIN}/python"
RESULTS_ROOT="${RESULTS_ROOT:-${ROOT}/results_3v3}"
THREADS="${THREADS:-16}"
USE_CONTROLS="${USE_CONTROLS:-1}"
METADATA_DIR="${RESULTS_ROOT}/metadata"
SUMMARY_DIR="${RESULTS_ROOT}/summary"
MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/peakforge-mpl}"

if [[ ! -x "${PYTHON_BIN}" ]]; then
  PYTHON_BIN="$(command -v python3)"
fi

if [[ -z "${PYTHON_BIN}" ]]; then
  echo "Python environment not found" >&2
  exit 1
fi

export PATH="${ENV_BIN}:${PATH}"
export MPLCONFIGDIR
mkdir -p "${MPLCONFIGDIR}"

declare -A CONDITIONS=(
  [K562_rep1]="K562"
  [K562_rep2]="K562"
  [K562_rep3]="K562"
  [HepG2_rep1]="HepG2"
  [HepG2_rep2]="HepG2"
  [HepG2_rep3]="HepG2"
)

declare -A BAMS=(
  [K562_rep1]="${ROOT}/data/K562_rep1.bam"
  [K562_rep2]="${ROOT}/data/K562_rep2.bam"
  [K562_rep3]="${ROOT}/data2/ENCFF953SVM.bam"
  [HepG2_rep1]="${ROOT}/data/HepG2_rep1.bam"
  [HepG2_rep2]="${ROOT}/data/HepG2_rep2.bam"
  [HepG2_rep3]="${ROOT}/data2/ENCFF439SHL.bam"
)

declare -A CONTROLS=(
  [K562_rep1]="${ROOT}/data/K562_input.bam"
  [K562_rep2]="${ROOT}/data/K562_input.bam"
  [K562_rep3]="${ROOT}/data/K562_input.bam"
  [HepG2_rep1]="${ROOT}/data/HepG2_input.bam"
  [HepG2_rep2]="${ROOT}/data/HepG2_input.bam"
  [HepG2_rep3]="${ROOT}/data/HepG2_input.bam"
)

declare -A PEAKS=(
  [K562_rep1]="${ROOT}/results/2v2/peaks/K562_rep1_peaks.narrowPeak"
  [K562_rep2]="${ROOT}/results/2v2/peaks/K562_rep2_peaks.narrowPeak"
  [K562_rep3]=""
  [HepG2_rep1]="${ROOT}/results/2v2/peaks/HepG2_rep1_peaks.narrowPeak"
  [HepG2_rep2]="${ROOT}/results/2v2/peaks/HepG2_rep2_peaks.narrowPeak"
  [HepG2_rep3]=""
)

if [[ "${USE_CONTROLS}" != "1" ]]; then
  PEAKS[K562_rep1]=""
  PEAKS[K562_rep2]=""
  PEAKS[K562_rep3]=""
  PEAKS[HepG2_rep1]=""
  PEAKS[HepG2_rep2]=""
  PEAKS[HepG2_rep3]=""
fi

require_file() {
  local path="$1"
  if [[ ! -f "${path}" ]]; then
    echo "Missing required file: ${path}" >&2
    exit 1
  fi
}

ensure_bam_index() {
  local bam="$1"
  if [[ -f "${bam}.bai" || -f "${bam%.bam}.bai" ]]; then
    return
  fi
  echo "[lopo] Indexing ${bam}"
  samtools index -@ "${THREADS}" "${bam}"
}

write_metadata() {
  local output_path="$1"
  shift
  mkdir -p "$(dirname "${output_path}")"
  printf "sample\tcondition\tbam\tcontrol_bam\tpeaks\tpeak_type\n" > "${output_path}"
  local sample peak_path control_bam
  for sample in "$@"; do
    peak_path="${PEAKS[${sample}]}"
    if [[ -n "${peak_path}" && ! -f "${peak_path}" ]]; then
      peak_path=""
    fi
    control_bam="${CONTROLS[${sample}]}"
    if [[ "${USE_CONTROLS}" != "1" ]]; then
      control_bam=""
    fi
    printf "%s\t%s\t%s\t%s\t%s\tnarrow\n" \
      "${sample}" \
      "${CONDITIONS[${sample}]}" \
      "${BAMS[${sample}]}" \
      "${control_bam}" \
      "${peak_path}" >> "${output_path}"
  done
}

run_chipdiff() {
  local metadata_path="$1"
  local output_dir="$2"
  local min_overlap="$3"
  local consensus_path="${4:-}"

  if [[ -f "${output_dir}/differential_results.tsv" ]]; then
    echo "[lopo] Reusing existing results in ${output_dir}"
    return
  fi

  mkdir -p "${output_dir}" "${output_dir}/peaks"
  echo "[lopo] ${metadata_path} -> ${output_dir}"

  local cmd=(
    "${PYTHON_BIN}" "${PROJECT_ROOT}/chipdiff.py"
    tsvmode "${metadata_path}"
    --output-dir "${output_dir}"
    --peak-dir "${output_dir}/peaks"
    --peak-type narrow
    --peak-extension 250
    --min-overlap "${min_overlap}"
    --macs2-genome hs
    --threads "${THREADS}"
  )

  if [[ -n "${consensus_path}" ]]; then
    cmd+=(--consensus-peaks "${consensus_path}")
  fi

  "${cmd[@]}"
}

refresh_peak_cache() {
  local output_dir="$1"
  shift
  local sample peak_path
  for sample in "$@"; do
    peak_path="${output_dir}/peaks/${sample}_peaks.narrowPeak"
    if [[ -f "${peak_path}" ]]; then
      PEAKS["${sample}"]="${peak_path}"
    fi
  done
}

prepare_fold() {
  local fold="$1"
  local -n gold_samples_ref="$2"
  local -n heldout_samples_ref="$3"

  local gold_metadata="${METADATA_DIR}/fold${fold}_gold.tsv"
  local heldout_metadata="${METADATA_DIR}/fold${fold}_heldout.tsv"
  local gold_dir="${RESULTS_ROOT}/fold${fold}_gold_2v2"
  local heldout_dir="${RESULTS_ROOT}/fold${fold}_heldout_1v1"

  write_metadata "${gold_metadata}" "${gold_samples_ref[@]}"
  run_chipdiff "${gold_metadata}" "${gold_dir}" 2
  refresh_peak_cache "${gold_dir}" "${gold_samples_ref[@]}"

  write_metadata "${heldout_metadata}" "${heldout_samples_ref[@]}"
  run_chipdiff "${heldout_metadata}" "${heldout_dir}" 1 "${gold_dir}/consensus_peaks.bed"
}

main() {
  mkdir -p "${RESULTS_ROOT}" "${METADATA_DIR}" "${SUMMARY_DIR}"

  local sample
  for sample in "${!BAMS[@]}"; do
    require_file "${BAMS[${sample}]}"
  done
  if [[ "${USE_CONTROLS}" == "1" ]]; then
    for sample in "${!CONTROLS[@]}"; do
      require_file "${CONTROLS[${sample}]}"
    done
    for sample in K562_rep1 K562_rep2 HepG2_rep1 HepG2_rep2; do
      require_file "${PEAKS[${sample}]}"
    done
  fi

  ensure_bam_index "${BAMS[K562_rep3]}"
  ensure_bam_index "${BAMS[HepG2_rep3]}"

  local fold1_gold=(K562_rep2 K562_rep3 HepG2_rep2 HepG2_rep3)
  local fold1_heldout=(K562_rep1 HepG2_rep1)
  local fold2_gold=(K562_rep1 K562_rep3 HepG2_rep1 HepG2_rep3)
  local fold2_heldout=(K562_rep2 HepG2_rep2)
  local fold3_gold=(K562_rep1 K562_rep2 HepG2_rep1 HepG2_rep2)
  local fold3_heldout=(K562_rep3 HepG2_rep3)

  prepare_fold 1 fold1_gold fold1_heldout
  prepare_fold 2 fold2_gold fold2_heldout
  prepare_fold 3 fold3_gold fold3_heldout

  "${PYTHON_BIN}" "${ROOT}/summarize_lopo_validation.py" \
    --results-root "${RESULTS_ROOT}" \
    --output-dir "${SUMMARY_DIR}"

  echo "[lopo] Summary written to ${SUMMARY_DIR}"
}

main "$@"
