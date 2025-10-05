#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
OUT_DIR="${1:-${PROJECT_ROOT}/results/prior_demo}"

mkdir -p "${OUT_DIR}"
PRIOR_DIR="${OUT_DIR}/prior_inputs"
mkdir -p "${PRIOR_DIR}"

PRIOR_BED="${PRIOR_DIR}/demo_prior.bed"
cat > "${PRIOR_BED}" <<'EOF'
chr1 100000 100400
chr1 200000 200450
chr1 300000 300520
EOF

PRIOR_SHAPE="${PRIOR_DIR}/demo_shape.json"
cat <<'EOF' > "${PRIOR_SHAPE}"
{
  "FWHM": {"mean": 320.0, "std": 85.0},
  "core_flank_ratio": {"mean": 4.2, "std": 1.3},
  "centroid": {"mean": 0.0, "std": 25.0},
  "skewness": {"mean": 0.0, "std": 0.6}
}
EOF

MANIFEST="${PRIOR_DIR}/prior_manifest.json"
cat <<EOF > "${MANIFEST}"
{
  "prior_bed": "${PRIOR_BED}",
  "prior_stats": "${PRIOR_SHAPE}",
  "prior_weight": 0.4
}
EOF

echo "[prior-demo] Running PeakForge differential pipeline with priors"
"${PROJECT_ROOT}/peakforge" tsvmode "${PROJECT_ROOT}/example/data/metadata_1v1.tsv" \
  --output-dir "${OUT_DIR}" \
  --prior-manifest "${MANIFEST}" \
  --prior-weight 0.4

PEAK_SHAPE_DIR="${PROJECT_ROOT}/example/peak_shape"
PEAK_SHAPE_DATA="${PEAK_SHAPE_DIR}/data"
if [[ ! -f "${PEAK_SHAPE_DATA}/demo_sample_A.bw" || ! -f "${PEAK_SHAPE_DATA}/demo_sample_B.bw" ]]; then
  echo "[prior-demo] Generating synthetic peak shape tracks" >&2
  python "${PEAK_SHAPE_DIR}/make_demo_tracks.py" --output "${PEAK_SHAPE_DATA}"
fi

SHAPE_OUT="${OUT_DIR}/shape"
mkdir -p "${SHAPE_OUT}"

echo "[prior-demo] Running peak shape comparison with priors"
"${PROJECT_ROOT}/peakforge" peakshape \
  --bigwig-a "${PEAK_SHAPE_DATA}/demo_sample_A.bw" \
  --bigwig-b "${PEAK_SHAPE_DATA}/demo_sample_B.bw" \
  --bed "${PEAK_SHAPE_DATA}/demo_regions.bed" \
  --core 400 \
  --flank 800 2000 \
  --out "${SHAPE_OUT}" \
  --threads 1 \
  --prior-shape "${PRIOR_SHAPE}" \
  --prior-weight 0.4

echo "[prior-demo] Prior-aware outputs available in ${OUT_DIR}"

