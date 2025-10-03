#!/usr/bin/env bash
set -euo pipefail

show_help() {
  cat <<'USAGE'
Usage: run_example2.sh [options]

Run the PeakForge pipeline using already-downloaded BAM files and optional
pre-computed peak calls. This helper creates a temporary metadata sheet for
you so you can launch the full pipeline without re-downloading or
re-indexing inputs.

Required arguments:
  --condition-a NAME         Label for condition A (e.g. K562)
  --a-bams FILE [FILE ...]   One or more BAMs belonging to condition A
  --condition-b NAME         Label for condition B (e.g. HepG2)
  --b-bams FILE [FILE ...]   One or more BAMs belonging to condition B

Optional arguments:
  --a-peaks FILE [FILE ...]  Peak files aligned with --a-bams (summits.bed,
                             narrowPeak, or broadPeak). When provided, MACS2
                             is skipped for those samples.
  --b-peaks FILE [FILE ...]  Peak files aligned with --b-bams.
  --output-dir DIR           Output directory (default: example/results/example2)
  --threads N                Threads for multiBamSummary (default: 16)
  --min-overlap N            Minimum samples required for consensus peaks (default: 2)
  --peak-type TYPE           Default peak type when calling MACS2 (default: narrow)
  --summit-extension BP      Summit extension/half-window size (default: 250)
  --macs2-genome G           Genome size string passed to MACS2 (default: hs)
  -h, --help                 Show this help message and exit

Examples:
  # 1 vs 1 with provided summits
  bash run_example2.sh \
    --condition-a K562 \
    --a-bams data/K562_rep1.bam \
    --a-peaks results/2v2/peaks/K562_rep1_summits.bed \
    --condition-b HepG2 \
    --b-bams data/HepG2_rep1.bam \
    --b-peaks results/2v2/peaks/HepG2_rep1_summits.bed

  # 2 vs 2 using replicate peaks
  bash run_example2.sh \
    --condition-a K562 \
    --a-bams data/K562_rep1.bam data/K562_rep2.bam \
    --a-peaks results/2v2/peaks/K562_rep1_summits.bed results/2v2/peaks/K562_rep2_summits.bed \
    --condition-b HepG2 \
    --b-bams data/HepG2_rep1.bam data/HepG2_rep2.bam \
    --b-peaks results/2v2/peaks/HepG2_rep1_summits.bed results/2v2/peaks/HepG2_rep2_summits.bed
USAGE
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

COND_A=""
COND_B=""
OUTPUT_DIR="${SCRIPT_DIR}/results/example2"
THREADS=16
MIN_OVERLAP=2
PEAK_TYPE="narrow"
SUMMIT_EXTENSION=250
MACS2_GENOME="hs"
A_BAMS=()
B_BAMS=()
A_PEAKS=()
B_PEAKS=()

if [[ $# -eq 0 ]]; then
  show_help
  exit 1
fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    --condition-a)
      COND_A="$2"
      shift 2
      ;;
    --condition-b)
      COND_B="$2"
      shift 2
      ;;
    --a-bams)
      shift
      while [[ $# -gt 0 && $1 != --* ]]; do
        A_BAMS+=("$1")
        shift
      done
      ;;
    --b-bams)
      shift
      while [[ $# -gt 0 && $1 != --* ]]; do
        B_BAMS+=("$1")
        shift
      done
      ;;
    --a-peaks)
      shift
      while [[ $# -gt 0 && $1 != --* ]]; do
        A_PEAKS+=("$1")
        shift
      done
      ;;
    --b-peaks)
      shift
      while [[ $# -gt 0 && $1 != --* ]]; do
        B_PEAKS+=("$1")
        shift
      done
      ;;
    --output-dir)
      OUTPUT_DIR="$2"
      shift 2
      ;;
    --threads)
      THREADS="$2"
      shift 2
      ;;
    --min-overlap)
      MIN_OVERLAP="$2"
      shift 2
      ;;
    --peak-type)
      PEAK_TYPE="$2"
      shift 2
      ;;
    --summit-extension)
      SUMMIT_EXTENSION="$2"
      shift 2
      ;;
    --macs2-genome)
      MACS2_GENOME="$2"
      shift 2
      ;;
    -h|--help)
      show_help
      exit 0
      ;;
    --*)
      echo "Unknown option: $1" >&2
      show_help >&2
      exit 1
      ;;
    *)
      echo "Unexpected argument: $1" >&2
      show_help >&2
      exit 1
      ;;
  esac
done

if [[ -z "${COND_A}" || -z "${COND_B}" ]]; then
  echo "Both --condition-a and --condition-b are required" >&2
  exit 1
fi

if [[ ${#A_BAMS[@]} -eq 0 || ${#B_BAMS[@]} -eq 0 ]]; then
  echo "You must supply at least one BAM for each condition" >&2
  exit 1
fi

if [[ ${#A_PEAKS[@]} -gt 0 && ${#A_PEAKS[@]} -ne ${#A_BAMS[@]} ]]; then
  echo "Number of --a-peaks entries must match --a-bams" >&2
  exit 1
fi

if [[ ${#B_PEAKS[@]} -gt 0 && ${#B_PEAKS[@]} -ne ${#B_BAMS[@]} ]]; then
  echo "Number of --b-peaks entries must match --b-bams" >&2
  exit 1
fi

mkdir -p "${OUTPUT_DIR}" "${OUTPUT_DIR}/peaks"

CMD=(
  python "${PROJECT_ROOT}/chipdiff.py" runmode
  --condition-a "${COND_A}"
  --a-bams "${A_BAMS[@]}"
  --condition-b "${COND_B}"
  --b-bams "${B_BAMS[@]}"
  --output-dir "${OUTPUT_DIR}"
  --peak-dir "${OUTPUT_DIR}/peaks"
  --min-overlap "${MIN_OVERLAP}"
  --peak-type "${PEAK_TYPE}"
  --summit-extension "${SUMMIT_EXTENSION}"
  --macs2-genome "${MACS2_GENOME}"
  --threads "${THREADS}"
)

if [[ ${#A_PEAKS[@]} -gt 0 ]]; then
  CMD+=(--a-peaks "${A_PEAKS[@]}")
fi

if [[ ${#B_PEAKS[@]} -gt 0 ]]; then
  CMD+=(--b-peaks "${B_PEAKS[@]}")
fi

echo "[example2] Launching PeakForge..." >&2
"${CMD[@]}"

echo "[example2] Results available in ${OUTPUT_DIR}" >&2
