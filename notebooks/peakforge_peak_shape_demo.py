# %%
"""PeakForge peak shape demo analysis.

This notebook-style script generates the synthetic bigWig tracks that ship with
PeakForge, runs the `peak_shape` module to compute shape descriptors, and
visualises the resulting metrics.  Each section is separated by `# %%` markers so
that the file can be imported into Jupyter or executed with tools that understand
cell-style comments.
"""

from __future__ import annotations

from pathlib import Path
import importlib.util
import json
import sys
from typing import Dict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyBigWig

# %%
# Paths and configuration
REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
EXAMPLE_DIR = REPO_ROOT / "example" / "peak_shape"
DATA_DIR = EXAMPLE_DIR / "data"
OUTPUT_DIR = EXAMPLE_DIR / "notebook_results"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
PLOTS_DIR = OUTPUT_DIR / "notebook_plots"
PLOTS_DIR.mkdir(parents=True, exist_ok=True)

import peak_shape

plt.style.use("seaborn-v0_8")

# %%
# Helper to import the demo track generator without turning the directory into a package
def load_make_demo_tracks():
    spec = importlib.util.spec_from_file_location(
        "make_demo_tracks", EXAMPLE_DIR / "make_demo_tracks.py"
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


make_demo_tracks = load_make_demo_tracks()

# %%
# Generate the synthetic bigWig tracks and BED intervals
bw_a, bw_b = make_demo_tracks.create_demo_tracks(DATA_DIR)
bed_path = DATA_DIR / "demo_regions.bed"

print("Demo tracks ready:")
print(f"  Sample A: {bw_a}")
print(f"  Sample B: {bw_b}")
print(f"  Regions : {bed_path}")

# %%
# Execute the peak shape workflow using the generated tracks
args = peak_shape.parse_args(
    [
        "--bigwig-a",
        str(bw_a),
        "--bigwig-b",
        str(bw_b),
        "--bed",
        str(bed_path),
        "--core",
        "500",
        "--flank",
        "1000",
        "3000",
        "--threads",
        "1",
        "--out",
        str(OUTPUT_DIR),
    ]
)
peak_shape.run_peak_shape(args)

results_path = OUTPUT_DIR / "peak_shape.tsv"
print(f"Metrics table written to {results_path}")

# %%
# Load the metrics table
metrics = pd.read_csv(results_path, sep="\t")
metrics

# %%
# Summarise the delta metrics between samples A and B
delta_columns = [
    "delta_FWHM",
    "delta_core_flank",
    "delta_centroid",
    "delta_skewness",
]
delta_summary = metrics[delta_columns].describe().T
print("Delta metric summary (B - A):")
delta_summary

# %%
# Visualise relationships between the delta metrics
fig, ax = plt.subplots(figsize=(8, 6))
scatter = ax.scatter(
    metrics["delta_core_flank"],
    metrics["delta_centroid"],
    c=metrics["delta_FWHM"],
    cmap="coolwarm",
    s=80,
    edgecolor="black",
)
ax.set_title("PeakForge peak shape deltas")
ax.set_xlabel("Δ core:flank ratio (B - A)")
ax.set_ylabel("Δ centroid (B - A, bp)")
cbar = fig.colorbar(scatter, ax=ax, label="Δ FWHM (bp)")
fig.tight_layout()
fig.savefig(PLOTS_DIR / "delta_relationships.png", dpi=300)
plt.show()

# %%
# Compute an outlier score similar to the CLI heatmap logic
z_scores = []
for column in delta_columns:
    values = metrics[column].to_numpy(dtype=float)
    mask = np.isfinite(values)
    z = np.zeros_like(values)
    if np.any(mask):
        subset = values[mask]
        std = subset.std(ddof=0)
        if std > 0:
            mean = subset.mean()
            z[mask] = (subset - mean) / std
    z_scores.append(np.abs(z))
outlier_score = np.sum(z_scores, axis=0)
metrics["outlier_score"] = outlier_score

metrics_sorted = metrics.sort_values("outlier_score", ascending=False)
print("Top peaks by combined delta z-score:")
metrics_sorted[["peak_id", "outlier_score", *delta_columns]].head(5)

# %%
# Plot the normalised signal profiles for the top three outliers
def fetch_signal(path: Path, chrom: str, start: int, end: int) -> np.ndarray:
    with pyBigWig.open(str(path)) as bw:
        values = bw.values(chrom, start, end, numpy=True)
    signal = np.asarray(values, dtype=float)
    signal = np.nan_to_num(signal, nan=0.0, posinf=0.0, neginf=0.0)
    signal[signal < 0] = 0.0
    return signal


def normalise_signal(signal: np.ndarray) -> np.ndarray:
    total = signal.sum()
    if not np.isfinite(total) or total <= 0:
        return np.zeros_like(signal)
    return signal / total


def plot_peak_signals(row: pd.Series) -> None:
    chrom = row["chrom"]
    start = int(row["start"])
    end = int(row["end"])
    span = end - start
    signal_a = fetch_signal(Path(bw_a), chrom, start, end)
    signal_b = fetch_signal(Path(bw_b), chrom, start, end)
    norm_a = normalise_signal(signal_a)
    norm_b = normalise_signal(signal_b)

    positions = np.linspace(-span / 2, span / 2, num=norm_a.size)
    plt.figure(figsize=(8, 4))
    plt.plot(positions, norm_a, label="Sample A", linewidth=2)
    plt.plot(positions, norm_b, label="Sample B", linewidth=2)
    plt.title(f"Normalised signal profiles for {row['peak_id']}")
    plt.xlabel("Position relative to region centre (bp)")
    plt.ylabel("Normalised signal")
    plt.legend()
    plt.tight_layout()
    plot_path = PLOTS_DIR / f"signal_{row['peak_id'].replace(':', '_').replace('-', '_')}.png"
    plt.savefig(plot_path, dpi=300)
    plt.show()


for _, peak_row in metrics_sorted.head(3).iterrows():
    plot_peak_signals(peak_row)

# %%
# Persist a JSON snapshot of the top peaks for downstream reporting
top_records: Dict[str, Dict[str, float]] = {}
for _, row in metrics_sorted.head(10).iterrows():
    top_records[row["peak_id"]] = {
        "delta_FWHM": float(row["delta_FWHM"]),
        "delta_core_flank": float(row["delta_core_flank"]),
        "delta_centroid": float(row["delta_centroid"]),
        "delta_skewness": float(row["delta_skewness"]),
        "outlier_score": float(row["outlier_score"]),
    }

json_path = OUTPUT_DIR / "top_peak_deltas.json"
json_path.write_text(json.dumps(top_records, indent=2))
print(f"Saved top peak delta summary to {json_path}")

# %%
print("Notebook-style analysis complete.")
