#!/usr/bin/env python3
"""chipdiff.py

End-to-end pipeline for CUT&Tag / ChIP-seq differential analysis.

The script expects a sample sheet describing the input BAM files and
(optional) pre-computed peak files.  The sample sheet must be a tab- or
comma-delimited text file with the following required columns:

    sample    Unique sample identifier (no spaces)
    condition    Experimental condition or group label
    bam    Path to the aligned reads in BAM format

Optional columns:

    peaks    Path to an existing peak file (summits.bed, narrowPeak, broadPeak)
    peak_type    One of {auto, narrow, broad, summit}.  ``auto`` (default)
                 attempts to infer the peak type from the file name.

Example TSV sample sheet::

    sample  condition   bam                 peaks               peak_type
    S1      treated     data/S1.bam         data/S1_summits.bed summit
    S2      treated     data/S2.bam         -                   -
    C1      control     data/C1.bam         data/C1_peaks.bed   narrow

The pipeline performs the following steps:

1. Peak calling with MACS2 (if required).
2. Construction of consensus peaks across samples.
3. Counting read overlaps per consensus peak using deepTools
   ``multiBamSummary``.
4. Differential analysis with PyDESeq2 (replicated designs) or the
   MARS method (no replicates).
5. Optional annotation against a GTF file and Enrichr enrichment via
   gseapy.
6. Plot generation (volcano, MA, correlation, heatmap) and metadata
   capture.

Dependencies: numpy, pandas, scipy, statsmodels, matplotlib, seaborn,
pyranges, gseapy, MACS2, deepTools.
"""
from __future__ import annotations

import argparse
import json
import logging
import math
import shutil
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

import numpy as np
import pandas as pd
import pyranges as pr
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from scipy import stats
from statsmodels.nonparametric.smoothers_lowess import lowess

try:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
except ImportError:  # pragma: no cover - optional dependency
    DeseqDataSet = None  # type: ignore[assignment]
    DeseqStats = None  # type: ignore[assignment]

try:
    import gseapy
except ImportError:  # pragma: no cover - optional dependency
    gseapy = None


# ---------------------------------------------------------------------------
# Data classes and utility helpers
# ---------------------------------------------------------------------------


@dataclass
class SampleEntry:
    """Representation of a single sample entry from the metadata sheet."""

    sample: str
    condition: str
    bam: Path
    peaks: Optional[Path] = None
    peak_type: str = "auto"

    def ensure_paths(self) -> None:
        if not self.bam.exists():
            raise FileNotFoundError(f"BAM file not found for sample {self.sample}: {self.bam}")
        if self.peaks is not None and not self.peaks.exists():
            raise FileNotFoundError(f"Peak file not found for sample {self.sample}: {self.peaks}")


# ---------------------------------------------------------------------------
# File and command helpers
# ---------------------------------------------------------------------------


def ensure_commands(commands: Sequence[str]) -> None:
    missing = [cmd for cmd in commands if shutil.which(cmd) is None]
    if missing:
        joined = ", ".join(sorted(missing))
        raise RuntimeError(
            "Missing required command(s): "
            f"{joined}. Please install them (e.g. via 'conda install -c bioconda macs2 deeptools')."
        )


def run_command(cmd: Sequence[str], *, workdir: Optional[Path] = None, log: bool = True) -> None:
    """Run a subprocess command with logging and error handling."""

    if log:
        logging.info("Running command: %s", " ".join(cmd))
    result = subprocess.run(cmd, cwd=str(workdir) if workdir else None, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed with exit code {result.returncode}: {' '.join(cmd)}")


def ensure_directory(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def read_table(path: Path) -> pd.DataFrame:
    """Read a delimited table inferring delimiter automatically."""

    try:
        df = pd.read_csv(path, sep=None, engine="python")
    except Exception as exc:  # pragma: no cover - passthrough error
        raise RuntimeError(f"Failed to read metadata file {path}: {exc}")
    return df


# ---------------------------------------------------------------------------
# Metadata parsing
# ---------------------------------------------------------------------------


def load_samples(metadata_path: Path) -> List[SampleEntry]:
    df = read_table(metadata_path)
    required = {"sample", "condition", "bam"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Metadata file missing required columns: {', '.join(sorted(missing))}")

    entries: List[SampleEntry] = []
    for _, row in df.iterrows():
        peaks_val = row.get("peaks")
        peaks = Path(str(peaks_val)) if isinstance(peaks_val, str) and peaks_val.strip() not in {"", "-", "NA", "None"} else None
        peak_type = str(row.get("peak_type", "auto")).lower()
        entry = SampleEntry(
            sample=str(row["sample"]),
            condition=str(row["condition"]),
            bam=Path(str(row["bam"])),
            peaks=peaks,
            peak_type=peak_type,
        )
        entry.ensure_paths()
        entries.append(entry)
    return entries


# ---------------------------------------------------------------------------
# Peak handling
# ---------------------------------------------------------------------------


def infer_peak_type(path: Path, declared: str, default: str) -> str:
    if declared and declared not in {"", "auto", "nan"}:
        if declared not in {"narrow", "broad", "summit"}:
            raise ValueError(f"Unknown peak_type '{declared}' for file {path}")
        return declared
    name = path.name.lower()
    if "summit" in name:
        return "summit"
    if name.endswith(".broadpeak"):
        return "broad"
    if name.endswith(".narrowpeak"):
        return "narrow"
    return default


def read_peak_file(path: Path, peak_type: str, summit_extension: int) -> pr.PyRanges:
    """Load peaks into a :class:`pyranges.PyRanges` object."""

    df = pd.read_csv(path, sep="\t", comment="#", header=None, dtype={0: str})
    if df.shape[1] < 3:
        raise ValueError(f"Peak file {path} must have at least 3 columns")
    df = df.iloc[:, :3]
    df.columns = ["Chromosome", "Start", "End"]

    if peak_type == "summit":
        center = (df["Start"].astype(int)).to_numpy()
        df["Start"] = np.maximum(center - summit_extension, 0)
        df["End"] = center + summit_extension
    elif peak_type == "narrow":
        start = df["Start"].astype(int) - summit_extension
        end = df["End"].astype(int) + summit_extension
        df["Start"] = np.maximum(start, 0)
        df["End"] = end
    else:  # broad
        df["Start"] = df["Start"].astype(int)
        df["End"] = df["End"].astype(int)

    return pr.PyRanges(df)


def call_macs2(sample: SampleEntry, *, output_dir: Path, macs2_genome: str, macs2_qval: float,
               peak_type: str, macs2_extra: Optional[List[str]] = None) -> Path:
    """Call MACS2 for a sample and return the resulting peak file path."""

    ensure_directory(output_dir)
    macs2_extra = macs2_extra or []
    name = sample.sample
    out_prefix = output_dir / name
    cmd = [
        "macs2",
        "callpeak",
        "-t",
        str(sample.bam),
        "-n",
        str(out_prefix),
        "-g",
        macs2_genome,
        "-q",
        str(macs2_qval),
    ]
    if peak_type == "broad":
        cmd.extend(["--broad"])
    cmd.extend(macs2_extra)
    run_command(cmd)

    if peak_type == "broad":
        peak_path = output_dir / f"{name}_peaks.broadPeak"
    elif peak_type == "summit":
        peak_path = output_dir / f"{name}_summits.bed"
    else:
        # Prefer summits if available, otherwise narrowPeak
        summit_path = output_dir / f"{name}_summits.bed"
        peak_path = summit_path if summit_path.exists() else output_dir / f"{name}_peaks.narrowPeak"
    if not peak_path.exists():
        raise FileNotFoundError(f"MACS2 output not found for sample {sample.sample}: {peak_path}")
    return peak_path


def load_all_peaks(samples: List[SampleEntry], *, summit_extension: int, default_peak_type: str,
                   macs2_params: Dict[str, str], peak_output_dir: Path) -> Dict[str, pr.PyRanges]:
    """Ensure every sample has peak calls and return PyRanges per sample."""

    peak_ranges: Dict[str, pr.PyRanges] = {}
    for sample in samples:
        if sample.peaks is None:
            peak_type = sample.peak_type if sample.peak_type != "auto" else default_peak_type
            logging.info("Calling MACS2 for sample %s (type=%s)", sample.sample, peak_type)
            peak_path = call_macs2(
                sample,
                output_dir=peak_output_dir,
                macs2_genome=macs2_params["genome"],
                macs2_qval=float(macs2_params["qvalue"]),
                peak_type=peak_type,
                macs2_extra=macs2_params.get("extra", []),
            )
        else:
            peak_type = infer_peak_type(sample.peaks, sample.peak_type, default_peak_type)
            peak_path = sample.peaks
            logging.info("Using provided peaks for sample %s (%s)", sample.sample, peak_type)

        pr_obj = read_peak_file(Path(peak_path), peak_type, summit_extension)
        df = pr_obj.df
        df["Sample"] = sample.sample
        pr_obj = pr.PyRanges(df)
        peak_ranges[sample.sample] = pr_obj
    return peak_ranges


def build_consensus(peak_ranges: Dict[str, pr.PyRanges], *, min_overlap: int) -> pr.PyRanges:
    """Build consensus peaks across samples with minimum overlap criteria."""

    logging.info("Building consensus peaks across %d samples", len(peak_ranges))
    combined = pr.concat(list(peak_ranges.values()))
    clustered = combined.cluster()
    df = clustered.df

    consensus_rows = []
    for cluster_id, group in df.groupby("Cluster"):
        samples = group["Sample"].unique()
        if len(samples) < min_overlap:
            continue
        chrom = group["Chromosome"].iloc[0]
        start = int(group["Start"].min())
        end = int(group["End"].max())
        consensus_rows.append((chrom, start, end, len(samples)))

    consensus_df = pd.DataFrame(consensus_rows, columns=["Chromosome", "Start", "End", "Support"])
    consensus_df.sort_values(["Chromosome", "Start", "End"], inplace=True)
    consensus_df.reset_index(drop=True, inplace=True)
    consensus_df["Name"] = [f"consensus_{i+1}" for i in range(len(consensus_df))]
    return pr.PyRanges(consensus_df[["Chromosome", "Start", "End", "Name", "Support"]])


# ---------------------------------------------------------------------------
# Counting with deepTools
# ---------------------------------------------------------------------------


def run_multibamsummary(consensus_bed: Path, samples: List[SampleEntry], output_dir: Path,
                        threads: int = 1) -> Path:
    ensure_directory(output_dir)
    out_npz = output_dir / "counts.npz"
    out_tsv = output_dir / "counts.tsv"

    cmd = [
        "multiBamSummary",
        "BED-file",
        "--BED",
        str(consensus_bed),
        "--bamfiles",
    ]
    cmd.extend(str(sample.bam) for sample in samples)
    cmd.extend([
        "--outFileName",
        str(out_npz),
        "--outRawCounts",
        str(out_tsv),
        "--numberOfProcessors",
        str(threads),
    ])
    run_command(cmd)
    if not out_tsv.exists():
        raise FileNotFoundError("multiBamSummary failed to produce counts TSV")
    return out_tsv


# ---------------------------------------------------------------------------
# Differential analysis utilities
# ---------------------------------------------------------------------------


def benjamini_hochberg(pvalues: pd.Series) -> pd.Series:
    pvals = pvalues.fillna(1.0).to_numpy(dtype=float, copy=True)
    n = len(pvals)
    if n == 0:
        return pd.Series(index=pvalues.index, dtype=float)

    order = np.argsort(pvals)
    ranks = np.arange(1, n + 1, dtype=float)
    adjusted_sorted = pvals[order] * n / ranks
    adjusted_sorted = np.minimum.accumulate(adjusted_sorted[::-1])[::-1]

    adjusted = np.empty_like(adjusted_sorted)
    adjusted[order] = adjusted_sorted
    adjusted = np.clip(adjusted, 0, 1)
    return pd.Series(adjusted, index=pvalues.index)


def pydeseq2_differential(counts: pd.DataFrame, conditions: pd.Series) -> pd.DataFrame:
    if DeseqDataSet is None or DeseqStats is None:
        raise ImportError("pydeseq2 is required for the DESeq2 workflow but is not installed")

    logging.info("Running PyDESeq2 differential analysis")
    samples = counts.columns.tolist()
    cond = conditions.loc[samples]
    if cond.nunique() != 2:
        raise ValueError("PyDESeq2 differential analysis requires exactly two conditions")

    condition_order = list(dict.fromkeys(cond.tolist()))
    reference = condition_order[0]
    contrast = condition_order[1]

    metadata = pd.DataFrame({"condition": cond.astype("category")}, index=samples)
    dds = DeseqDataSet(counts=counts.T.astype(int), metadata=metadata, design="~condition")
    dds.deseq2()

    stats = DeseqStats(dds, contrast=("condition", contrast, reference))
    stats.summary()
    res = stats.results_df.copy()
    res.index.name = "Peak"

    result = pd.DataFrame(index=res.index)
    if "baseMean" in res.columns:
        result["baseMean"] = res["baseMean"]
    result["log2FC"] = res["log2FoldChange"]
    if "lfcSE" in res.columns:
        result["lfcSE"] = res["lfcSE"]
    if "stat" in res.columns:
        result["waldStat"] = res["stat"]
    result["pvalue"] = res["pvalue"]
    result["padj"] = res["padj"].fillna(1.0)
    result["log2FC_shrunk"] = result["log2FC"]
    result["method"] = "pydeseq2"
    return result


def mars_differential(counts: pd.DataFrame, conditions: pd.Series) -> pd.DataFrame:
    """Implement the MARS method for designs without replicates."""

    logging.info("Running MARS differential analysis (no replicates)")
    unique_conditions = conditions.unique()
    if len(unique_conditions) != 2:
        raise ValueError("MARS method requires exactly two conditions")

    cond_a, cond_b = unique_conditions
    counts_a = counts.loc[:, conditions[conditions == cond_a].index].sum(axis=1)
    counts_b = counts.loc[:, conditions[conditions == cond_b].index].sum(axis=1)

    x = counts_a.to_numpy() + 1.0
    y = counts_b.to_numpy() + 1.0
    M = np.log2(y / x)
    A = 0.5 * np.log2(x * y)

    fitted = lowess(M, A, frac=0.3, return_sorted=False)
    residuals = M - fitted

    # Estimate variance as a smooth function of A using binning
    bins = max(10, int(np.sqrt(len(A))))
    df = pd.DataFrame({"A": A, "residual": residuals})
    df["bin"] = pd.cut(df["A"], bins, duplicates="drop")
    var_by_bin = df.groupby("bin")["residual"].var().fillna(df["residual"].var())
    bin_centers = df.groupby("bin")["A"].mean()
    # Interpolate variance for each observation
    interp_var = np.interp(A, bin_centers.fillna(0).to_numpy(), var_by_bin.to_numpy(), left=var_by_bin.iloc[0], right=var_by_bin.iloc[-1])
    z_scores = residuals / np.sqrt(interp_var + 1e-6)
    pvals = 2 * stats.norm.sf(np.abs(z_scores))

    res_df = pd.DataFrame({
        "Peak": counts.index,
        "log2FC": M,
        "A": A,
        "pvalue": pvals,
        "log2FC_shrunk": M * (interp_var / (interp_var + 1.0)),
    }).set_index("Peak")
    res_df["padj"] = benjamini_hochberg(res_df["pvalue"]).fillna(1.0)
    res_df["method"] = "mars"
    return res_df


# ---------------------------------------------------------------------------
# Differential orchestration helpers
# ---------------------------------------------------------------------------


def call_differential_analysis(counts: pd.DataFrame, conditions: pd.Series) -> pd.DataFrame:
    """Select and run the appropriate differential analysis workflow."""

    if counts.empty:
        raise ValueError("Counts matrix is empty; cannot perform differential analysis")
    if conditions.nunique() != 2:
        raise ValueError("Differential analysis requires exactly two experimental conditions")

    replicates_per_condition = conditions.value_counts().min()
    if replicates_per_condition >= 2:
        logging.info("Detected replicates per condition; using DESeq2 analysis via PyDESeq2")
        return pydeseq2_differential(counts, conditions)

    logging.info("No replicates detected; using MARS analysis")
    return mars_differential(counts, conditions)


# ---------------------------------------------------------------------------
# Annotation and enrichment
# ---------------------------------------------------------------------------


def annotate_peaks(consensus: pr.PyRanges, gtf_path: Path) -> pd.DataFrame:
    logging.info("Annotating peaks using GTF: %s", gtf_path)
    gtf = pr.read_gtf(gtf_path)
    genes = gtf[gtf.Feature == "gene"]
    nearest = consensus.nearest(genes, how="nearest")
    df = nearest.df
    gene_col = "gene_name" if "gene_name" in df.columns else "Name_b"
    distance_col = "Distance" if "Distance" in df.columns else df.filter(like="Distance").columns[0]
    df = df.rename(columns={gene_col: "NearestGene", distance_col: "Distance"})
    keep_cols = [col for col in ["Chromosome", "Start", "End", "Name", "Support", "NearestGene", "Distance"] if col in df.columns]
    return df[keep_cols]


def run_enrichr(genes: Sequence[str], out_dir: Path, description: str = "top_genes") -> Optional[Path]:
    if gseapy is None:
        logging.warning("gseapy not available; skipping enrichment analysis")
        return None
    ensure_directory(out_dir)
    try:
        enr = gseapy.enrichr(
            gene_list=list(genes),
            gene_sets=["GO_Biological_Process_2021"],
            description=description,
            outdir=str(out_dir),
            cutoff=0.5,
        )
    except Exception as exc:  # pragma: no cover - network dependent
        logging.warning("Enrichr analysis failed: %s", exc)
        return None
    report = Path(enr.res2d_path) if hasattr(enr, "res2d_path") else None
    return report


# ---------------------------------------------------------------------------
# Plotting utilities
# ---------------------------------------------------------------------------


def save_plot(fig: plt.Figure, path: Path) -> None:
    ensure_directory(path.parent)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_volcano(results: pd.DataFrame, output: Path, padj_threshold: float = 0.05,
                 lfc_threshold: float = 1.0) -> None:
    fig, ax = plt.subplots(figsize=(6, 6))
    res = results.copy()
    padj_nonzero = res.loc[res["padj"] > 0, "padj"]
    min_padj = padj_nonzero.min() if not padj_nonzero.empty else 1e-6
    res["-log10(padj)"] = -np.log10(res["padj"].replace(0, min_padj / 10))
    sns.scatterplot(data=res, x="log2FC", y="-log10(padj)", ax=ax, hue=res["padj"] < padj_threshold,
                    palette={True: "red", False: "grey"}, legend=False)
    ax.axvline(lfc_threshold, color="black", linestyle="--", linewidth=0.8)
    ax.axvline(-lfc_threshold, color="black", linestyle="--", linewidth=0.8)
    ax.axhline(-math.log10(padj_threshold), color="black", linestyle="--", linewidth=0.8)
    ax.set_title("Volcano plot")
    ax.set_xlabel("log2 Fold Change")
    ax.set_ylabel("-log10 adjusted p-value")
    save_plot(fig, output)


def plot_ma(results: pd.DataFrame, counts: pd.DataFrame, output: Path) -> None:
    base_mean = results.get("baseMean")
    if base_mean is None:
        base_mean = counts.mean(axis=1)
    fig, ax = plt.subplots(figsize=(6, 6))
    A = np.log2(base_mean + 1e-6)
    M = results["log2FC"]
    sns.scatterplot(x=A, y=M, ax=ax, s=10, color="steelblue")
    ax.axhline(0, color="black", linestyle="--", linewidth=0.8)
    ax.set_xlabel("Average log2 expression")
    ax.set_ylabel("log2 Fold Change")
    ax.set_title("MA plot")
    save_plot(fig, output)


def plot_sample_correlation(counts: pd.DataFrame, output: Path) -> None:
    corr = counts.apply(lambda x: np.log1p(x)).corr()
    fig, ax = plt.subplots(figsize=(6, 6))
    sns.heatmap(corr, cmap="viridis", ax=ax, annot=True, fmt=".2f")
    ax.set_title("Sample correlation")
    save_plot(fig, output)


def plot_top_heatmap(counts: pd.DataFrame, results: pd.DataFrame, output: Path, top_n: int = 50) -> None:
    top = results.sort_values("padj").head(top_n).index
    data = counts.loc[top]
    log_data = np.log2(data + 1)
    norm = log_data.sub(log_data.mean(axis=1), axis=0)
    fig = plt.figure(figsize=(8, max(4, len(top) * 0.2)))
    sns.heatmap(norm, cmap="RdBu_r", center=0)
    plt.title(f"Top {top_n} differential peaks")
    plt.ylabel("Peaks")
    plt.xlabel("Samples")
    save_plot(fig, output)


def plot_differential_summary(results: pd.DataFrame, output: Path, *, counts: Optional[pd.DataFrame] = None,
                              top_n: int = 20) -> None:
    """Create an overview scatter plot inspired by clusterProfiler dotplots."""

    required = {"padj", "log2FC"}
    missing = [col for col in required if col not in results.columns]
    if missing:
        logging.warning("Cannot draw differential summary plot; missing columns: %s", ", ".join(missing))
        return

    df = results.replace([np.inf, -np.inf], np.nan).dropna(subset=["padj", "log2FC"])
    if df.empty:
        logging.warning("No finite differential results available for summary plot")
        return

    df = df.sort_values("padj").head(top_n).copy()
    if df.empty:
        logging.warning("Differential results contain no entries after filtering for top peaks")
        return

    padj_nonzero = df.loc[df["padj"] > 0, "padj"]
    min_nonzero = padj_nonzero.min() if not padj_nonzero.empty else 1e-6
    df["padj"] = df["padj"].replace(0, min_nonzero / 10)
    df["Peak"] = df.index.astype(str)

    if "baseMean" in df.columns:
        base_mean = df["baseMean"].astype(float)
    elif counts is not None and not counts.empty:
        base_mean = counts.mean(axis=1).reindex(df.index)
    else:
        base_mean = pd.Series(1.0, index=df.index)

    if base_mean.isna().all():
        base_mean = pd.Series(1.0, index=df.index)
    base_mean = base_mean.fillna(base_mean.median() if base_mean.notna().any() else 1.0)
    base_mean = base_mean.clip(lower=1e-3)
    df["MeanCount"] = base_mean
    df["MeanCountDisplay"] = np.log10(df["MeanCount"] + 1.0)
    df["neg_log10_padj"] = -np.log10(df["padj"])

    fc_values = df["log2FC"].astype(float)
    fc_max = np.nanmax(np.abs(fc_values.to_numpy()))
    if not np.isfinite(fc_max) or fc_max == 0:
        fc_max = 1.0
    hue_norm = Normalize(vmin=-fc_max, vmax=fc_max)

    fig, ax = plt.subplots(figsize=(8, max(4.0, len(df) * 0.45)))
    size_range = (80, 700)
    scatter = sns.scatterplot(
        data=df,
        x="neg_log10_padj",
        y="Peak",
        size="MeanCountDisplay",
        hue="log2FC",
        palette="RdBu_r",
        sizes=size_range,
        hue_norm=hue_norm,
        ax=ax,
        edgecolor="black",
        linewidth=0.5,
    )
    ax.set_xlabel("-log10 adjusted p-value")
    ax.set_ylabel("Peak")
    ax.set_title("Differential peak landscape")
    ax.grid(axis="x", linestyle="--", linewidth=0.5, alpha=0.4)

    sm = plt.cm.ScalarMappable(cmap="RdBu_r", norm=hue_norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label("log2 Fold Change")

    min_display, max_display = df["MeanCount"].min(), df["MeanCount"].max()
    if math.isclose(min_display, max_display):
        legend_handles = [Line2D([0], [0], marker="o", linestyle="", color="black",
                                  markersize=math.sqrt(np.mean(size_range)), markerfacecolor="none")]
        legend_labels = [f"{min_display:.1f}"]
    else:
        legend_values = np.linspace(min_display, max_display, num=3)
        legend_values = np.unique(np.round(legend_values, 2))
        display_min, display_max = df["MeanCountDisplay"].min(), df["MeanCountDisplay"].max()

        def _size_for(val: float) -> float:
            display_val = math.log10(val + 1.0)
            if math.isclose(display_min, display_max):
                return float(np.mean(size_range))
            return float(np.interp(display_val, [display_min, display_max], size_range))

        legend_handles = [
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="",
                color="black",
                markersize=math.sqrt(_size_for(val)),
                markerfacecolor="none",
            )
            for val in legend_values
        ]
        legend_labels = [f"{val:.1f}" for val in legend_values]

    size_legend = ax.legend(
        legend_handles,
        legend_labels,
        title="Mean count",
        loc="lower right",
        frameon=False,
    )
    ax.add_artist(size_legend)

    save_plot(fig, output)


def generate_differential_plots(results: pd.DataFrame, counts: pd.DataFrame, output_dir: Path) -> Dict[str, Path]:
    """Produce all differential analysis visualisations."""

    ensure_directory(output_dir)
    outputs = {
        "sample_correlation": output_dir / "sample_correlation.png",
        "ma": output_dir / "ma_plot.png",
        "volcano": output_dir / "volcano.png",
        "top_heatmap": output_dir / "top_peaks_heatmap.png",
        "summary": output_dir / "differential_summary.png",
    }

    plot_sample_correlation(counts, outputs["sample_correlation"])
    plot_ma(results, counts, outputs["ma"])
    plot_volcano(results, outputs["volcano"])
    plot_top_heatmap(counts, results, outputs["top_heatmap"])
    plot_differential_summary(results, outputs["summary"], counts=counts)

    return outputs


# ---------------------------------------------------------------------------
# Metadata persistence
# ---------------------------------------------------------------------------


def save_metadata(metadata: Dict, output_path: Path) -> None:
    ensure_directory(output_path.parent)
    with output_path.open("w") as fh:
        json.dump(metadata, fh, indent=2)


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------


def run_pipeline(args: argparse.Namespace) -> None:
    samples = load_samples(Path(args.metadata))
    conditions = pd.Series({s.sample: s.condition for s in samples})

    required_cmds = ["multiBamSummary"]
    needs_peak_calling = any(sample.peaks is None for sample in samples)
    if needs_peak_calling:
        required_cmds.append("macs2")
    ensure_commands(required_cmds)

    macs2_params = {
        "genome": args.macs2_genome,
        "qvalue": args.macs2_qvalue,
        "extra": args.macs2_extra,
    }
    peak_ranges = load_all_peaks(
        samples,
        summit_extension=args.summit_extension,
        default_peak_type=args.peak_type,
        macs2_params=macs2_params,
        peak_output_dir=Path(args.peak_dir),
    )

    consensus = build_consensus(peak_ranges, min_overlap=args.min_overlap)
    results_dir = ensure_directory(Path(args.output_dir))
    consensus_bed = results_dir / "consensus_peaks.bed"
    consensus.df[["Chromosome", "Start", "End", "Name"]].to_csv(consensus_bed, sep="\t", header=False, index=False)

    counts_tsv = run_multibamsummary(consensus_bed, samples, results_dir / "counts", threads=args.threads)
    raw_counts = pd.read_csv(counts_tsv, sep="\t")

    # deepTools 3.5 switched to quoting header labels in TSV output.  Clean up
    # the column names so downstream logic can match against the expected
    # ``chrom`` / ``start`` / ``end`` headers and sample BAM names regardless
    # of whether they were quoted or prefixed with ``#``.
    def _normalise_header(value: str) -> str:
        cleaned = value.strip()
        cleaned = cleaned.lstrip("#")
        cleaned = cleaned.strip("'\"")
        return cleaned or value

    raw_counts.rename(columns={col: _normalise_header(col) for col in raw_counts.columns}, inplace=True)

    # deepTools historically used ``#chr`` for the chromosome column but some
    # versions emit ``#chrom`` or even plain ``chrom``.  Normalise these headers
    # so downstream joins can rely on a canonical ``Chromosome`` column.
    chromosome_aliases = {
        "#chrom",
        "#chr",
        "chrom",
        "chr",
        "#Chromosome",
        "Chromosome",
    }
    matched_chrom = next((col for col in raw_counts.columns if col in chromosome_aliases), None)
    if matched_chrom is None:
        matched_chrom = next(
            (
                col
                for col in raw_counts.columns
                if col.lower().lstrip("#") in {"chrom", "chr", "chromosome"}
            ),
            None,
        )
    if matched_chrom is None:
        raise ValueError(
            "Counts matrix is missing a chromosome column (expected one of #chr, #chrom, chrom)."
        )
    raw_counts.rename(columns={matched_chrom: "Chromosome"}, inplace=True)

    consensus_df = consensus.df.rename(columns={"Start": "start", "End": "end"})
    merged = raw_counts.merge(consensus_df[["Chromosome", "start", "end", "Name"]],
                              on=["Chromosome", "start", "end"], how="left")
    merged["Peak"] = merged["Name"].fillna(
        merged["Chromosome"].astype(str)
        + ":"
        + merged["start"].astype(int).astype(str)
        + "-"
        + merged["end"].astype(int).astype(str)
    )
    count_cols = [col for col in merged.columns if col not in {"Chromosome", "start", "end", "Peak", "Name", "Support"}]
    counts_df = merged.set_index("Peak")[count_cols]
    column_map: Dict[str, str] = {}
    for sample in samples:
        bam_path = Path(sample.bam)
        column_map[str(sample.bam)] = sample.sample
        column_map[bam_path.name] = sample.sample
        column_map[bam_path.stem] = sample.sample
    counts_df = counts_df.rename(columns=column_map)
    missing_cols = [s.sample for s in samples if s.sample not in counts_df.columns]
    if missing_cols:
        raise ValueError(f"Counts matrix missing columns for samples: {', '.join(missing_cols)}")
    counts_df = counts_df[[s.sample for s in samples]]

    diff_res = call_differential_analysis(counts_df, conditions)

    diff_path = results_dir / "differential_results.tsv"
    diff_res.to_csv(diff_path, sep="\t")

    plot_dir = results_dir / "plots"
    plot_paths = generate_differential_plots(diff_res, counts_df, plot_dir)

    annotation_df = None
    annotation_path = None
    if args.gtf:
        annotation_df = annotate_peaks(consensus, Path(args.gtf))
        annotation_path = results_dir / "consensus_annotation.tsv"
        annotation_df.to_csv(annotation_path, sep="\t", index=False)

    enrichr_path = None
    if args.enrichr:
        if gseapy is None:
            logging.warning("gseapy not installed; skipping Enrichr analysis")
        elif annotation_df is None or "NearestGene" not in annotation_df.columns:
            logging.warning("Annotation required for Enrichr; provide --gtf to map peaks to genes")
        else:
            top_peaks = diff_res.sort_values("padj").head(args.enrichr_top).index
            gene_map = annotation_df.set_index("Name")["NearestGene"].dropna()
            top_genes = gene_map.reindex(top_peaks).dropna().unique().tolist()
            if not top_genes:
                logging.warning("No genes available for Enrichr after annotation")
            else:
                enrichr_path = run_enrichr(top_genes, results_dir / "enrichr", description="chipdiff")

    sample_metadata = []
    for sample in samples:
        sample_metadata.append(
            {
                "sample": sample.sample,
                "condition": sample.condition,
                "bam": str(sample.bam),
                "peaks": str(sample.peaks) if sample.peaks is not None else None,
                "peak_type": sample.peak_type,
            }
        )

    metadata = {
        "timestamp": datetime.utcnow().isoformat(),
        "args": vars(args),
        "samples": sample_metadata,
        "counts_matrix": str(counts_tsv),
        "differential_results": str(diff_path),
        "plots": {key: str(path) for key, path in plot_paths.items()},
        "annotation": str(annotation_path) if annotation_path else None,
        "enrichr": str(enrichr_path) if enrichr_path else None,
    }
    save_metadata(metadata, results_dir / "metadata.json")


# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="CUT&Tag / ChIP-seq differential analysis pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--metadata", required=True, help="Sample sheet (TSV/CSV)")
    parser.add_argument("--output-dir", default="results", help="Output directory")
    parser.add_argument("--peak-dir", default="peaks", help="Directory for peak calls")
    parser.add_argument("--peak-type", default="narrow", choices=["narrow", "broad", "summit"],
                        help="Default peak type when calling MACS2")
    parser.add_argument("--summit-extension", type=int, default=250, help="Extension for summits/narrow peaks (bp)")
    parser.add_argument("--min-overlap", type=int, default=2, help="Minimum samples required for consensus peak")
    parser.add_argument("--macs2-genome", default="hs", help="MACS2 genome size (e.g. hs, mm, 2.7e9)")
    parser.add_argument("--macs2-qvalue", type=float, default=0.01, help="MACS2 q-value cutoff")
    parser.add_argument("--macs2-extra", nargs=argparse.REMAINDER, default=[], help="Additional arguments for MACS2")
    parser.add_argument(
        "--threads",
        type=int,
        default=16,
        help="Threads for multiBamSummary (--numberOfProcessors)",
    )
    parser.add_argument("--gtf", help="Optional GTF file for annotation")
    parser.add_argument("--enrichr", action="store_true", help="Run Enrichr GO Biological Process analysis")
    parser.add_argument("--enrichr-top", type=int, default=200, help="Number of top peaks for enrichment")
    parser.add_argument("--log-level", default="INFO", help="Logging level")
    return parser


# ---------------------------------------------------------------------------
# Entrypoint
# ---------------------------------------------------------------------------


def main(argv: Optional[Sequence[str]] = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level.upper(), logging.INFO),
                        format="[%(asctime)s] %(levelname)s: %(message)s")

    try:
        run_pipeline(args)
    except Exception as exc:  # pragma: no cover - CLI exception reporting
        logging.error("Pipeline failed: %s", exc)
        sys.exit(1)


if __name__ == "__main__":
    main()
