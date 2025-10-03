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
4. Differential analysis with either a DESeq2-like NB-GLM workflow
   (replicated designs) or the MARS method (no replicates).
5. Optional annotation against a GTF file and Enrichr enrichment via
   gseapy.
6. Plot generation (volcano, MA, correlation, heatmap) and metadata
   capture.

Dependencies: numpy, pandas, scipy, statsmodels, matplotlib, seaborn,
pyranges, gseapy, MACS2, deepTools.
"""
from __future__ import annotations

import argparse
import dataclasses
import json
import logging
import math
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import pyranges as pr
import seaborn as sns
from matplotlib import pyplot as plt
from scipy import stats
from statsmodels.nonparametric.smoothers_lowess import lowess
import statsmodels.api as sm

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


def median_ratio_normalization(counts: pd.DataFrame) -> Tuple[pd.DataFrame, np.ndarray]:
    """Perform median ratio normalization (DESeq2)."""

    logging.info("Performing median-ratio normalization")
    counts = counts.astype(float)
    values = counts.to_numpy()
    with np.errstate(divide="ignore"):
        log_counts = np.where(values > 0, np.log(values), np.nan)
    geometric_means = np.exp(np.nanmean(log_counts, axis=1))
    geometric_means[~np.isfinite(geometric_means)] = 1.0
    ratios = counts.divide(geometric_means, axis=0).replace([np.inf, -np.inf], np.nan)
    size_factors = ratios.median(axis=0, skipna=True)
    size_factors = size_factors.replace(0, np.nan).fillna(1.0).to_numpy()
    normalized = counts.divide(size_factors, axis=1)
    return normalized, size_factors


def estimate_dispersions(normalized: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Estimate raw dispersions, fitted trend, and shrunk dispersions."""

    means = normalized.mean(axis=1)
    variances = normalized.var(axis=1, ddof=1)
    raw_disp = (variances - means) / (means ** 2)
    raw_disp = raw_disp.clip(lower=1e-8)

    log_means = np.log(means + 1e-8)
    log_disp = np.log(raw_disp)
    # Fit LOWESS trend
    fitted = lowess(log_disp, log_means, frac=0.3, return_sorted=False)
    trend = np.exp(fitted)
    prior_weight = 10.0
    shrunk = (raw_disp + prior_weight * trend) / (1.0 + prior_weight)
    return raw_disp.to_numpy(), trend, shrunk


def wald_test_deseq(counts: pd.DataFrame, conditions: pd.Series,
                     size_factors: np.ndarray, dispersions: np.ndarray) -> pd.DataFrame:
    samples = counts.columns.tolist()
    design = pd.get_dummies(conditions.loc[samples], drop_first=True)
    if design.shape[1] != 1:
        raise ValueError("DESeq2-like workflow currently supports exactly two conditions with replicates")
    cond_col = design.columns[0]
    X = sm.add_constant(design.iloc[:, 0].to_numpy())

    results = []
    offsets = np.log(size_factors)
    for idx, (peak, row) in enumerate(counts.iterrows()):
        y = row.to_numpy()
        disp = float(dispersions[idx])
        try:
            model = sm.GLM(y, X, family=sm.families.NegativeBinomial(alpha=disp), offset=offsets)
            fit = model.fit(maxiter=200, disp=0)
            beta = fit.params[1]
            se = fit.bse[1]
            wald = fit.tvalues[1]
            pval = 2 * stats.norm.sf(abs(wald))
        except Exception as exc:
            logging.debug("GLM failed for %s (%s), falling back to Wald on means", peak, exc)
            group_a = design.index[design.iloc[:, 0] == 0]
            group_b = design.index[design.iloc[:, 0] == 1]
            mean_a = row[group_a].mean() + 1e-6
            mean_b = row[group_b].mean() + 1e-6
            beta = math.log(mean_b / mean_a)
            se = math.sqrt(1.0 / mean_a + 1.0 / mean_b)
            wald = beta / se
            pval = 2 * stats.norm.sf(abs(wald))
        log2fc = beta / math.log(2)
        se_log2 = se / math.log(2)
        results.append((peak, log2fc, se_log2, pval))

    res_df = pd.DataFrame(results, columns=["Peak", "log2FC", "SE_log2FC", "pvalue"])
    res_df.set_index("Peak", inplace=True)
    return res_df


def shrink_log2fc(log2fc: pd.Series, se: pd.Series) -> pd.Series:
    prior_var = np.nanmedian(log2fc**2)
    prior_var = prior_var if np.isfinite(prior_var) and prior_var > 0 else 0.5
    shrinkage = prior_var / (prior_var + se**2)
    return log2fc * shrinkage


def benjamini_hochberg(pvalues: pd.Series) -> pd.Series:
    pvals = pvalues.fillna(1.0).to_numpy()
    n = len(pvals)
    order = np.argsort(pvals)
    ranks = np.empty(n, dtype=float)
    ranks[order] = np.arange(1, n + 1)
    adjusted = pvals * n / ranks
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    return pd.Series(adjusted, index=pvalues.index)


def deseqlike_differential(counts: pd.DataFrame, conditions: pd.Series) -> pd.DataFrame:
    normalized, size_factors = median_ratio_normalization(counts)
    raw_disp, trend, shrunk_disp = estimate_dispersions(normalized)
    res = wald_test_deseq(counts, conditions, size_factors, shrunk_disp)
    res["baseMean"] = normalized.mean(axis=1)
    res["lfcSE"] = res["SE_log2FC"]
    res["padj"] = benjamini_hochberg(res["pvalue"]).fillna(1.0)
    res["log2FC_shrunk"] = shrink_log2fc(res["log2FC"], res["SE_log2FC"])
    res["method"] = "deseq_like"
    return res


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
    raw_counts.rename(columns={"#chrom": "Chromosome"}, inplace=True)
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

    has_replicates = conditions.value_counts().min() >= 2
    if has_replicates:
        logging.info("Detected replicates per condition; using DESeq2-like analysis")
        diff_res = deseqlike_differential(counts_df, conditions)
    else:
        logging.info("No replicates detected; using MARS analysis")
        diff_res = mars_differential(counts_df, conditions)

    diff_path = results_dir / "differential_results.tsv"
    diff_res.to_csv(diff_path, sep="\t")

    plot_dir = results_dir / "plots"
    plot_sample_correlation(counts_df, plot_dir / "sample_correlation.png")
    plot_ma(diff_res, counts_df, plot_dir / "ma_plot.png")
    plot_volcano(diff_res, plot_dir / "volcano.png")
    plot_top_heatmap(counts_df, diff_res, plot_dir / "top_peaks_heatmap.png")

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

    metadata = {
        "timestamp": datetime.utcnow().isoformat(),
        "args": vars(args),
        "samples": [dataclasses.asdict(s) for s in samples],
        "counts_matrix": str(counts_tsv),
        "differential_results": str(diff_path),
        "plots": {
            "volcano": str((plot_dir / "volcano.png")),
            "ma": str((plot_dir / "ma_plot.png")),
            "sample_correlation": str((plot_dir / "sample_correlation.png")),
            "top_heatmap": str((plot_dir / "top_peaks_heatmap.png")),
        },
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
    parser.add_argument("--threads", type=int, default=4, help="Threads for multiBamSummary")
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
