#!/usr/bin/env python3
"""Peak shape analysis for paired bigWig samples.

This module computes several shape descriptors (FWHM, core:flank ratio,
centroid shift, and skewness) for each region provided in a BED file using two
bigWig tracks. The script also generates summary plots and an outlier heatmap
to visualise differences between the samples.
"""

from __future__ import annotations

import argparse
import atexit
import logging
import math
import os
import threading
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import pyBigWig
from matplotlib import pyplot as plt
from scipy import stats


LOGGER = logging.getLogger(__name__)


@dataclass
class Interval:
    """Container for BED interval information."""

    chrom: str
    start: int
    end: int
    name: str

    @property
    def peak_id(self) -> str:
        return f"{self.chrom}:{self.start}-{self.end}"

    @property
    def length(self) -> int:
        return self.end - self.start


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Peak shape analysis")
    parser.add_argument("--bigwig-a", required=True, help="Path to sample A bigWig")
    parser.add_argument("--bigwig-b", required=True, help="Path to sample B bigWig")
    parser.add_argument("--bed", required=True, help="BED file with regions")
    parser.add_argument(
        "--core",
        type=int,
        default=500,
        help="Half-width (bp) of the core window around the region centre",
    )
    parser.add_argument(
        "--flank",
        type=int,
        nargs=2,
        metavar=("MIN", "MAX"),
        default=(1000, 3000),
        help="Absolute distance range (bp) for the flank window",
    )
    parser.add_argument(
        "--fwhm-threshold",
        type=float,
        default=0.5,
        help="Fraction of the peak maximum used for FWHM determination",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of worker threads for peak processing",
    )
    parser.add_argument(
        "--out",
        default="results/shape",
        help="Output directory for tables and plots",
    )
    return parser.parse_args()


def load_bed(path: str) -> List[Interval]:
    intervals: List[Interval] = []
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) < 3:
                LOGGER.warning("Skipping malformed BED line: %s", line.strip())
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            if end <= start:
                LOGGER.warning("Skipping zero/negative length interval: %s", line.strip())
                continue
            name = parts[3] if len(parts) > 3 else f"{chrom}:{start}-{end}"
            intervals.append(Interval(chrom=chrom, start=start, end=end, name=name))
    if not intervals:
        raise ValueError(f"No valid intervals found in BED file: {path}")
    return intervals


BIGWIG_PATHS: Dict[str, str] = {}
THREAD_LOCAL = threading.local()
OPEN_HANDLES: List[pyBigWig.pyBigWig] = []
OPEN_HANDLES_LOCK = threading.Lock()


def _register_handle(handle: pyBigWig.pyBigWig) -> None:
    with OPEN_HANDLES_LOCK:
        if not any(existing is handle for existing in OPEN_HANDLES):
            OPEN_HANDLES.append(handle)


def _close_handles() -> None:
    with OPEN_HANDLES_LOCK:
        for handle in OPEN_HANDLES:
            try:
                handle.close()
            except Exception:  # pragma: no cover - cleanup best effort
                continue
        OPEN_HANDLES.clear()


atexit.register(_close_handles)


def get_bigwig(sample: str) -> pyBigWig.pyBigWig:
    handles = getattr(THREAD_LOCAL, "handles", None)
    if handles is None:
        handles = {}
        THREAD_LOCAL.handles = handles
    bw = handles.get(sample)
    if bw is None:
        bw = pyBigWig.open(BIGWIG_PATHS[sample])
        handles[sample] = bw
        _register_handle(bw)
    return bw


def fetch_signal(sample: str, interval: Interval) -> Optional[np.ndarray]:
    bw = get_bigwig(sample)
    try:
        values = bw.values(interval.chrom, interval.start, interval.end, numpy=True)
    except RuntimeError:
        LOGGER.warning(
            "%s missing coverage for %s (interval outside of chrom)", sample, interval.peak_id
        )
        return None
    if values is None:
        return None
    signal = np.asarray(values, dtype=float)
    if signal.size == 0:
        return None
    signal = np.nan_to_num(signal, nan=0.0, posinf=0.0, neginf=0.0)
    signal[signal < 0] = 0.0
    if np.all(signal == 0):
        return None
    return signal


def compute_metrics(
    signal: Optional[np.ndarray],
    *,
    core_half_width: int,
    flank_range: Tuple[int, int],
    fwhm_threshold: float,
) -> Dict[str, Optional[float]]:
    if signal is None:
        return {
            "FWHM": math.nan,
            "core_flank_ratio": math.nan,
            "centroid": math.nan,
            "skewness": math.nan,
            "norm_signal": None,
        }

    total = signal.sum()
    if not np.isfinite(total) or total <= 0:
        return {
            "FWHM": math.nan,
            "core_flank_ratio": math.nan,
            "centroid": math.nan,
            "skewness": math.nan,
            "norm_signal": None,
        }

    norm_signal = signal / total
    max_val = norm_signal.max() if norm_signal.size else 0.0
    if max_val <= 0 or not np.isfinite(max_val):
        fwhm = math.nan
    else:
        threshold = max_val * fwhm_threshold
        above = np.where(norm_signal >= threshold)[0]
        if above.size == 0:
            fwhm = math.nan
        else:
            fwhm = float(above[-1] - above[0] + 1)

    length = norm_signal.size
    centre = (length - 1) / 2.0
    positions = np.arange(length, dtype=float) - centre

    core_mask = np.abs(positions) <= core_half_width
    flank_min, flank_max = flank_range
    if flank_min > flank_max:
        flank_min, flank_max = flank_max, flank_min
    flank_mask = (np.abs(positions) >= flank_min) & (np.abs(positions) < flank_max)

    core_sum = float(norm_signal[core_mask].sum()) if np.any(core_mask) else 0.0
    flank_sum = float(norm_signal[flank_mask].sum()) if np.any(flank_mask) else 0.0
    ratio = core_sum / flank_sum if flank_sum > 0 else math.nan

    centroid = float(np.sum(positions * norm_signal))

    try:
        skewness = float(stats.skew(norm_signal, bias=False, nan_policy="omit"))
    except Exception:  # pragma: no cover - defensive
        skewness = math.nan

    return {
        "FWHM": fwhm,
        "core_flank_ratio": ratio,
        "centroid": centroid,
        "skewness": skewness,
        "norm_signal": norm_signal,
    }


def process_interval(interval: Interval, config: Dict[str, float]) -> Dict[str, object]:
    try:
        signal_a = fetch_signal("A", interval)
        signal_b = fetch_signal("B", interval)
        metrics_a = compute_metrics(
            signal_a,
            core_half_width=config["core_half_width"],
            flank_range=config["flank_range"],
            fwhm_threshold=config["fwhm_threshold"],
        )
        metrics_b = compute_metrics(
            signal_b,
            core_half_width=config["core_half_width"],
            flank_range=config["flank_range"],
            fwhm_threshold=config["fwhm_threshold"],
        )
    except Exception as exc:  # pragma: no cover - safety net
        LOGGER.exception("Failed to process %s: %s", interval.peak_id, exc)
        metrics_a = metrics_b = {
            "FWHM": math.nan,
            "core_flank_ratio": math.nan,
            "centroid": math.nan,
            "skewness": math.nan,
            "norm_signal": None,
        }

    return {
        "peak_id": interval.peak_id,
        "chrom": interval.chrom,
        "start": interval.start,
        "end": interval.end,
        "A_FWHM": metrics_a["FWHM"],
        "B_FWHM": metrics_b["FWHM"],
        "A_core_flank": metrics_a["core_flank_ratio"],
        "B_core_flank": metrics_b["core_flank_ratio"],
        "A_centroid": metrics_a["centroid"],
        "B_centroid": metrics_b["centroid"],
        "A_skewness": metrics_a["skewness"],
        "B_skewness": metrics_b["skewness"],
        "A_signal": metrics_a["norm_signal"],
        "B_signal": metrics_b["norm_signal"],
    }


def ensure_output_dirs(base_dir: str) -> Tuple[str, str]:
    table_dir = base_dir
    plots_dir = os.path.join(base_dir, "plots")
    os.makedirs(table_dir, exist_ok=True)
    os.makedirs(plots_dir, exist_ok=True)
    return table_dir, plots_dir


def create_histogram(data: pd.Series, title: str, xlabel: str, output_path: str) -> None:
    valid = data.replace([np.inf, -np.inf], np.nan).dropna()
    if valid.empty:
        LOGGER.warning("Skipping plot '%s' due to lack of valid data", title)
        return
    plt.figure(figsize=(8, 5))
    plt.hist(valid, bins=40, color="#4c72b0", alpha=0.8)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def create_scatter(data: pd.Series, title: str, ylabel: str, output_path: str) -> None:
    valid = data.replace([np.inf, -np.inf], np.nan).dropna()
    if valid.empty:
        LOGGER.warning("Skipping plot '%s' due to lack of valid data", title)
        return
    plt.figure(figsize=(8, 5))
    plt.scatter(range(len(valid)), valid, s=10, alpha=0.7, color="#dd8452")
    plt.title(title)
    plt.xlabel("Peak index")
    plt.ylabel(ylabel)
    plt.axhline(0, color="black", linewidth=0.8, linestyle="--")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def resample_signal(signal: np.ndarray, bins: int = 200) -> np.ndarray:
    if signal is None or signal.size == 0:
        return np.full(bins, np.nan)
    x_orig = np.linspace(0.0, 1.0, num=signal.size)
    x_new = np.linspace(0.0, 1.0, num=bins)
    return np.interp(x_new, x_orig, signal)


def build_outlier_heatmap(
    df: pd.DataFrame,
    plots_dir: str,
    max_peaks: int = 50,
    bins: int = 200,
) -> None:
    delta_cols = ["delta_FWHM", "delta_core_flank", "delta_centroid", "delta_skewness"]
    scores = np.zeros(len(df), dtype=float)
    for col in delta_cols:
        values = df[col].to_numpy(dtype=float)
        mask = np.isfinite(values)
        if not np.any(mask):
            continue
        subset = values[mask]
        mean = subset.mean()
        std = subset.std(ddof=0)
        if std == 0:
            continue
        z = np.zeros_like(values)
        z[mask] = (subset - mean) / std
        scores += np.abs(z)

    df = df.copy()
    df["outlier_score"] = scores
    df = df.sort_values("outlier_score", ascending=False)

    heatmap_rows: List[np.ndarray] = []
    labels: List[str] = []
    for _, row in df.head(max_peaks).iterrows():
        signal_a = row["A_signal"]
        signal_b = row["B_signal"]
        if signal_a is None or signal_b is None:
            continue
        resampled_a = resample_signal(signal_a, bins=bins)
        resampled_b = resample_signal(signal_b, bins=bins)
        diff = resampled_b - resampled_a
        heatmap_rows.append(diff)
        labels.append(row["peak_id"])

    if not heatmap_rows:
        LOGGER.warning("Skipping outlier heatmap (no valid signals)")
        return

    data = np.vstack(heatmap_rows)
    plt.figure(figsize=(10, max(4, len(labels) * 0.2)))
    masked = np.ma.masked_invalid(data)
    im = plt.imshow(masked, aspect="auto", cmap="RdBu_r", interpolation="nearest")
    plt.colorbar(im, label="Signal (B - A)")
    plt.yticks(range(len(labels)), labels, fontsize=6)
    plt.xlabel("Relative position (bins)")
    plt.ylabel("Peak")
    plt.title("Top outlier peaks (B - A signal difference)")
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "top_outlier_heatmap.png"), dpi=300)
    plt.close()


def main() -> None:
    args = parse_args()

    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
    LOGGER.info("Loading intervals from %s", args.bed)
    intervals = load_bed(args.bed)
    LOGGER.info("Loaded %d intervals", len(intervals))

    BIGWIG_PATHS["A"] = args.bigwig_a
    BIGWIG_PATHS["B"] = args.bigwig_b

    # Validate bigWig accessibility
    for sample, path in BIGWIG_PATHS.items():
        if not os.path.exists(path):
            raise FileNotFoundError(f"bigWig file for sample {sample} not found: {path}")
        with pyBigWig.open(path) as bw:
            if bw.chroms() is None:
                raise ValueError(f"Failed to read chromosome info from {path}")

    threads = max(1, int(args.threads))
    config = {
        "core_half_width": int(args.core),
        "flank_range": (int(args.flank[0]), int(args.flank[1])),
        "fwhm_threshold": float(args.fwhm_threshold),
    }

    LOGGER.info(
        "Processing peaks (threads=%d, core=±%dbp, flank=%s, FWHM threshold=%.2f)",
        threads,
        config["core_half_width"],
        config["flank_range"],
        config["fwhm_threshold"],
    )

    if threads == 1:
        results = [process_interval(interval, config) for interval in intervals]
    else:
        from concurrent.futures import ThreadPoolExecutor

        with ThreadPoolExecutor(max_workers=threads) as executor:
            results = list(executor.map(lambda iv: process_interval(iv, config), intervals))

    df = pd.DataFrame(results)
    df["delta_FWHM"] = df["B_FWHM"] - df["A_FWHM"]
    df["delta_core_flank"] = df["B_core_flank"] - df["A_core_flank"]
    df["delta_centroid"] = df["B_centroid"] - df["A_centroid"]
    df["delta_skewness"] = df["B_skewness"] - df["A_skewness"]

    output_dir, plots_dir = ensure_output_dirs(args.out)
    table_path = os.path.join(output_dir, "peak_shape.tsv")

    df_export = df.drop(columns=["A_signal", "B_signal"])
    df_export.to_csv(table_path, sep="\t", index=False, float_format="%.6f")
    LOGGER.info("Wrote metrics table to %s", table_path)

    create_histogram(
        df["delta_FWHM"],
        title="ΔFWHM distribution",
        xlabel="ΔFWHM (B - A, bp)",
        output_path=os.path.join(plots_dir, "delta_fwhm_hist.png"),
    )
    create_histogram(
        df["delta_core_flank"],
        title="Δcore:flank ratio distribution",
        xlabel="Δ(core:flank) (B - A)",
        output_path=os.path.join(plots_dir, "delta_core_flank_hist.png"),
    )
    create_scatter(
        df["delta_centroid"],
        title="Δcentroid scatter",
        ylabel="Δcentroid (B - A, bp)",
        output_path=os.path.join(plots_dir, "delta_centroid_scatter.png"),
    )

    build_outlier_heatmap(df, plots_dir=plots_dir, max_peaks=50, bins=200)

    LOGGER.info("Analysis complete")


if __name__ == "__main__":
    main()

