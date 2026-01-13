#!/usr/bin/env python3
"""chipdiff.py

End-to-end pipeline for CUT&Tag / ChIP-seq differential analysis.
"""
from __future__ import annotations

import json
import logging
import math
import shutil
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pyranges as pr
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from scipy import stats

import typer
from typing_extensions import Annotated
from rich.logging import RichHandler

try:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
except ImportError:
    DeseqDataSet = None
    DeseqStats = None

try:
    import gseapy
except ImportError:
    gseapy = None

import peak_shape
from io_utils import ensure_integer_columns, read_bed_frame
from prior_utils import PriorRegistry, load_prior_manifest

# Initialize Typer app
app = typer.Typer(
    name="peakforge",
    help="PeakForge: A modern ChIP-seq/CUT&Tag differential analysis toolkit.",
    add_completion=False,
    no_args_is_help=True,
    context_settings={"help_option_names": ["-h", "--help"]},
)


def _detect_macs_command() -> str:
    for candidate in ("macs2", "macs3"):
        if shutil.which(candidate):
            return candidate
    raise RuntimeError(
        "Missing required command(s): macs2, macs3. Install MACS via 'pip install macs3'"
    )


MACS_COMMAND: Optional[str] = None


# ---------------------------------------------------------------------------
# Data classes and utility helpers
# ---------------------------------------------------------------------------


@dataclass
class SampleEntry:
    sample: str
    condition: str
    bam: Path
    peaks: Optional[Path] = None
    peak_type: str = "auto"
    is_paired: Optional[bool] = None

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
            f"{joined}. Install MACS via 'pip install macs3' (or macs2) and deepTools via "
            "'pip install deeptools' "
            "and samtools via 'conda install -c bioconda samtools'."
        )


def get_macs_command() -> str:
    global MACS_COMMAND
    if MACS_COMMAND is None:
        MACS_COMMAND = _detect_macs_command()
    return MACS_COMMAND


def ensure_python_version(min_version: tuple[int, int] = (3, 10)) -> None:
    if sys.version_info < min_version:
        formatted = ".".join(str(part) for part in min_version)
        raise RuntimeError(
            f"PeakForge requires Python {formatted} or newer; detected {sys.version.split()[0]}"
        )


def run_command(cmd: Sequence[str], *, workdir: Optional[Path] = None, log: bool = True) -> None:
    if log:
        logging.info("Running command: %s", " ".join(cmd))
    result = subprocess.run(cmd, cwd=str(workdir) if workdir else None, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed with exit code {result.returncode}: {' '.join(cmd)}")


def ensure_directory(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def _bam_index_candidates(bam: Path) -> List[Path]:
    candidates: List[Path] = []
    candidates.append(Path(f"{bam}.bai"))
    if bam.suffix:
        candidates.append(bam.with_suffix(".bai"))
    seen: Set[Path] = set()
    unique: List[Path] = []
    for candidate in candidates:
        if candidate not in seen:
            unique.append(candidate)
            seen.add(candidate)
    return unique


def ensure_bam_index(bam: Path, samtools_path: str, threads: int = 1) -> None:
    for candidate in _bam_index_candidates(bam):
        if candidate.exists():
            return
    logging.info("Indexing BAM for library size estimation: %s", bam)
    cmd = [samtools_path, "index"]
    if threads > 1:
        cmd.extend(["-@", str(threads)])
    cmd.append(str(bam))
    run_command(cmd)


def detect_paired_end_bam(bam: Path, samtools_path: str, threads: int = 1) -> bool:
    cmd = [samtools_path, "view", "-c", "-f", "1"]
    if threads > 1:
        cmd.extend(["-@", str(threads)])
    cmd.append(str(bam))

    result = subprocess.run(cmd, check=False, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"samtools view failed for {bam} with exit code {result.returncode}: {result.stderr.strip()}"
        )
    try:
        count = int(result.stdout.strip() or 0)
    except ValueError as exc:
        raise RuntimeError(f"Unable to parse samtools view output for {bam}: {result.stdout!r}") from exc

    paired = count > 0
    logging.info("Detected %s BAM for %s", "paired-end" if paired else "single-end", bam)
    return paired


def bam_total_mapped_reads(bam: Path, samtools_path: str, threads: int = 1) -> int:
    ensure_bam_index(bam, samtools_path, threads)
    result = subprocess.run(
        [samtools_path, "idxstats", str(bam)],
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"samtools idxstats failed for {bam} with exit code {result.returncode}: {result.stderr.strip()}"
        )

    total = 0
    for line in result.stdout.splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 3:
            continue
        try:
            total += int(fields[2])
        except ValueError:
            continue

    if total <= 0:
        raise ValueError(f"Unable to determine mapped reads for BAM {bam}")
    return total


def compute_library_sizes(samples: Sequence[SampleEntry], samtools_path: str, threads: int = 1) -> pd.Series:
    sizes = {}
    for sample in samples:
        logging.info("Estimating library size for sample %s", sample.sample)
        sizes[sample.sample] = bam_total_mapped_reads(sample.bam, samtools_path, threads)
    return pd.Series(sizes, dtype=float)


def read_table(path: Path) -> pd.DataFrame:
    try:
        df = pd.read_csv(path, sep=None, engine="python")
    except Exception as exc:
        raise RuntimeError(f"Failed to read metadata file {path}: {exc}")
    return df


# ---------------------------------------------------------------------------
# Metadata parsing
# ---------------------------------------------------------------------------


def _normalise_optional_path(value: object) -> Optional[Path]:
    if isinstance(value, Path):
        return value
    if isinstance(value, str):
        stripped = value.strip()
        if stripped and stripped not in {"-", "NA", "None", "nan"}:
            return Path(stripped)
    return None


def load_samples(metadata_path: Path) -> List[SampleEntry]:
    df = read_table(metadata_path)
    required = {"sample", "condition", "bam"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Metadata file missing required columns: {', '.join(sorted(missing))}")

    entries: List[SampleEntry] = []
    for row in df.itertuples(index=False):
        sample = getattr(row, "sample")
        condition = getattr(row, "condition")
        bam = getattr(row, "bam")
        peaks = _normalise_optional_path(getattr(row, "peaks", None))
        peak_type_raw = getattr(row, "peak_type", "auto")
        peak_type = str(peak_type_raw).lower() if peak_type_raw is not None else "auto"
        entry = SampleEntry(
            sample=str(sample),
            condition=str(condition),
            bam=Path(str(bam)),
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
        if declared not in {"narrow", "broad"}:
            raise ValueError(
                f"Unknown peak_type '{declared}' for file {path}; only 'narrow' and 'broad' are supported"
            )
        return declared
    name = path.name.lower()
    if name.endswith(".broadpeak"):
        return "broad"
    if name.endswith(".narrowpeak"):
        return "narrow"
    if "summit" in name:
        raise ValueError(
            "Summit-only BED files are no longer supported; please provide a narrowPeak or broadPeak file"
        )
    return default


def read_peak_file(path: Path, peak_type: str, peak_extension: int) -> pr.PyRanges:
    frame = read_bed_frame(path)
    frame = ensure_integer_columns(frame, ("Start", "End"))

    if peak_type == "narrow":
        start = frame["Start"].to_numpy() - peak_extension
        end = frame["End"].to_numpy() + peak_extension
        frame["Start"] = np.maximum(start, 0)
        frame["End"] = end

    return pr.PyRanges(frame)


def _macs2_command(
    sample: SampleEntry,
    *,
    output_dir: Path,
    macs2_genome: str,
    macs2_qval: float,
    peak_type: str,
    macs2_extra: Optional[List[str]] = None,
) -> tuple[list[str], Path]:
    ensure_directory(output_dir)
    macs2_extra = macs2_extra or []
    name = sample.sample
    out_prefix = output_dir / name
    cmd = [
        get_macs_command(),
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
    if sample.is_paired:
        cmd.extend(["-f", "BAMPE"])

    if peak_type == "broad":
        cmd.extend(["--broad"])
    cmd.extend(macs2_extra)

    if peak_type == "broad":
        peak_path = output_dir / f"{name}_peaks.broadPeak"
    else:
        peak_path = output_dir / f"{name}_peaks.narrowPeak"
    return cmd, peak_path


def call_macs2(
    sample: SampleEntry,
    *,
    output_dir: Path,
    macs2_genome: str,
    macs2_qval: float,
    peak_type: str,
    macs2_extra: Optional[List[str]] = None,
) -> Path:
    cmd, peak_path = _macs2_command(
        sample,
        output_dir=output_dir,
        macs2_genome=macs2_genome,
        macs2_qval=macs2_qval,
        peak_type=peak_type,
        macs2_extra=macs2_extra,
    )
    run_command(cmd)

    if not peak_path.exists():
        raise FileNotFoundError(f"MACS2 output not found for sample {sample.sample}: {peak_path}")
    return peak_path


@dataclass
class Macs2Job:
    sample: SampleEntry
    peak_type: str
    peak_path: Path
    process: subprocess.Popen[str]


def load_all_peaks(
    samples: List[SampleEntry],
    *,
    peak_extension: int,
    default_peak_type: str,
    macs2_params: Dict[str, str],
    peak_output_dir: Path,
    prior_registry: Optional[PriorRegistry] = None,
) -> Dict[str, pr.PyRanges]:
    peak_ranges: Dict[str, pr.PyRanges] = {}
    macs2_jobs: List[Macs2Job] = []

    for sample in samples:
        if sample.peaks is None:
            peak_type = sample.peak_type if sample.peak_type != "auto" else default_peak_type
            logging.info("Launching MACS2 for sample %s (type=%s)", sample.sample, peak_type)
            cmd, peak_path = _macs2_command(
                sample,
                output_dir=peak_output_dir,
                macs2_genome=macs2_params["genome"],
                macs2_qval=float(macs2_params["qvalue"]),
                peak_type=peak_type,
                macs2_extra=macs2_params.get("extra", []),
            )
            process = subprocess.Popen(cmd)
            macs2_jobs.append(
                Macs2Job(
                    sample=sample,
                    peak_type=peak_type,
                    peak_path=peak_path,
                    process=process,
                )
            )
        else:
            peak_type = infer_peak_type(sample.peaks, sample.peak_type, default_peak_type)
            peak_path = sample.peaks
            logging.info("Using provided peaks for sample %s (%s)", sample.sample, peak_type)

            pr_obj = read_peak_file(Path(peak_path), peak_type, peak_extension)
            df = pr_obj.df
            df["Sample"] = sample.sample
            pr_obj = pr.PyRanges(df)
            peak_ranges[sample.sample] = pr_obj
            if prior_registry is not None and prior_registry.enabled:
                prior_registry.record_sample(sample.sample, pr_obj.df)

    macs2_results: List[tuple[Macs2Job, int]] = []
    for job in macs2_jobs:
        returncode = job.process.wait()
        macs2_results.append((job, returncode))

    failed = [job for job, code in macs2_results if code != 0]
    if failed:
        errors = ", ".join(f"{job.sample.sample} (exit {job.process.returncode})" for job in failed)
        raise RuntimeError(f"MACS2 failed for sample(s): {errors}")

    for job, _ in macs2_results:
        if not job.peak_path.exists():
            raise FileNotFoundError(
                f"MACS2 output not found for sample {job.sample.sample}: {job.peak_path}"
            )

        pr_obj = read_peak_file(Path(job.peak_path), job.peak_type, peak_extension)
        df = pr_obj.df
        df["Sample"] = job.sample.sample
        pr_obj = pr.PyRanges(df)
        peak_ranges[job.sample.sample] = pr_obj
        if prior_registry is not None and prior_registry.enabled:
            prior_registry.record_sample(job.sample.sample, pr_obj.df)
    return peak_ranges


def build_consensus(peak_ranges: Dict[str, pr.PyRanges], *, min_overlap: int) -> pr.PyRanges:
    logging.info("Building consensus peaks across %d samples", len(peak_ranges))
    if not peak_ranges:
        return pr.PyRanges()
    combined = pr.concat(list(peak_ranges.values()))
    clustered = combined.cluster()
    df = clustered.df

    grouped = (
        df.groupby("Cluster")
        .agg(
            Chromosome=("Chromosome", "first"),
            Start=("Start", "min"),
            End=("End", "max"),
            Support=("Sample", pd.Series.nunique),
        )
        .reset_index(drop=True)
    )

    consensus_df = grouped[grouped["Support"] >= max(1, min_overlap)].copy()
    consensus_df.sort_values(["Chromosome", "Start", "End"], inplace=True)
    consensus_df.reset_index(drop=True, inplace=True)
    consensus_df["Name"] = [f"consensus_{i + 1}" for i in range(len(consensus_df))]
    return pr.PyRanges(consensus_df[["Chromosome", "Start", "End", "Name", "Support"]])


def load_consensus_bed(path: Path) -> pr.PyRanges:
    if not path.exists():
        raise FileNotFoundError(f"Consensus BED file not found: {path}")

    df = pd.read_csv(path, sep="\t", comment="#", header=None)
    if df.shape[1] < 3:
        raise ValueError(
            f"Consensus BED {path} must have at least three columns (chrom, start, end)"
        )

    base = df.iloc[:, :3].copy()
    base.columns = ["Chromosome", "Start", "End"]
    base = ensure_integer_columns(base, ("Start", "End"))

    names: List[str] = []
    provided_names = df.iloc[:, 3] if df.shape[1] >= 4 else None
    seen: Set[str] = set()
    for idx in range(len(base)):
        value: Optional[str] = None
        if provided_names is not None:
            raw = provided_names.iloc[idx]
            if pd.notna(raw):
                raw_str = str(raw).strip()
                if raw_str:
                    value = raw_str
        if not value:
            value = f"consensus_{idx + 1}"
        candidate = value
        suffix = 1
        while candidate in seen:
            suffix += 1
            candidate = f"{value}_{suffix}"
        seen.add(candidate)
        names.append(candidate)

    base["Name"] = names

    support = pd.Series([pd.NA] * len(base))
    if df.shape[1] >= 5:
        support = pd.to_numeric(df.iloc[:, 4], errors="coerce")
    base["Support"] = support

    return pr.PyRanges(base[["Chromosome", "Start", "End", "Name", "Support"]])


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


def compute_prior_score(df):
    beta1, beta2, beta3 = 2.0, 1.0, 1.0

    def _extract(column: str) -> np.ndarray:
        if column in df:
            series = pd.to_numeric(df[column], errors="coerce").fillna(0.0)
            return series.to_numpy(dtype=float)
        return np.zeros(len(df), dtype=float)

    overlap = _extract("PriorOverlap")
    int_z = np.clip(np.abs(_extract("IntZ")), 0, 4)
    shape_z = np.clip(np.abs(_extract("ShapeZ")), 0, 4)

    return beta1 * overlap - beta2 * int_z - beta3 * shape_z


def compute_prior_weights(df, gamma=0.5):
    scores = compute_prior_score(df)
    if scores.size == 0:
        return scores

    median = np.median(scores)
    mad = np.median(np.abs(scores - median))
    if mad == 0:
        normalized = np.zeros_like(scores)
    else:
        normalized = (scores - median) / mad

    w_raw = np.exp(gamma * normalized)
    mean_weight = w_raw.mean()
    if mean_weight == 0:
        return np.ones_like(w_raw)
    return w_raw / mean_weight


def weighted_bh(p_values, weights, alpha=0.05):
    pvals = np.asarray(p_values, dtype=float)
    w = np.asarray(weights, dtype=float)

    if pvals.shape != w.shape:
        raise ValueError("p_values and weights must have the same shape")

    if pvals.size == 0:
        return pvals, pvals, np.asarray([], dtype=bool)

    pvals = np.clip(np.nan_to_num(pvals, nan=1.0, posinf=1.0, neginf=0.0), 0.0, 1.0)
    mean_weight = w.mean()
    if mean_weight <= 0:
        w_norm = np.ones_like(w)
    else:
        w_norm = w / mean_weight

    with np.errstate(divide="ignore", invalid="ignore"):
        p_weighted = np.divide(pvals, w_norm, out=np.ones_like(pvals), where=w_norm != 0)
    p_weighted = np.clip(p_weighted, 0.0, 1.0)

    order = np.argsort(p_weighted)
    ranks = np.arange(1, p_weighted.size + 1, dtype=float)
    adjusted_sorted = p_weighted[order] * p_weighted.size / ranks
    adjusted_sorted = np.minimum.accumulate(adjusted_sorted[::-1])[::-1]

    q_weighted = np.empty_like(adjusted_sorted)
    q_weighted[order] = adjusted_sorted
    q_weighted = np.clip(q_weighted, 0.0, 1.0)

    significant = q_weighted <= alpha

    return p_weighted, q_weighted, significant


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


def mars_differential(
    counts: pd.DataFrame, conditions: pd.Series, library_sizes: pd.Series
) -> pd.DataFrame:
    logging.info("Running MARS differential analysis (no replicates)")
    samples = counts.columns.tolist()
    cond_series = conditions.loc[samples]
    unique_conditions = list(dict.fromkeys(cond_series.tolist()))
    if len(unique_conditions) != 2:
        raise ValueError("MARS method requires exactly two conditions")

    reference, contrast = unique_conditions
    contrast_cols = cond_series[cond_series == contrast].index
    reference_cols = cond_series[cond_series == reference].index

    if contrast_cols.empty or reference_cols.empty:
        raise ValueError("Each condition must contribute at least one sample for MARS analysis")

    contrast_counts = counts.loc[:, contrast_cols].sum(axis=1)
    reference_counts = counts.loc[:, reference_cols].sum(axis=1)

    if not isinstance(library_sizes, pd.Series):
        library_sizes = pd.Series(library_sizes)

    library_sizes = library_sizes.reindex(samples)
    if library_sizes.isna().any():
        missing = library_sizes[library_sizes.isna()].index.tolist()
        raise ValueError(
            "Library size information missing for samples: " + ", ".join(missing)
        )

    contrast_total = float(library_sizes.loc[list(contrast_cols)].sum())
    reference_total = float(library_sizes.loc[list(reference_cols)].sum())

    if contrast_total <= 0 or reference_total <= 0:
        raise ValueError("Total read counts must be positive for both conditions in MARS analysis")

    c1 = contrast_counts.to_numpy(dtype=float)
    c2 = reference_counts.to_numpy(dtype=float)

    with np.errstate(divide="ignore"):
        log2_c1 = np.log2(c1)
        log2_c2 = np.log2(c2)

    valid_mask = np.isfinite(log2_c1) & np.isfinite(log2_c2)

    M = np.full_like(c1, np.nan, dtype=float)
    A = np.full_like(c1, np.nan, dtype=float)
    M[valid_mask] = log2_c1[valid_mask] - log2_c2[valid_mask]
    A[valid_mask] = 0.5 * (log2_c1[valid_mask] + log2_c2[valid_mask])

    sqrt_total = math.sqrt(contrast_total * reference_total)
    p = np.full_like(M, np.nan, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        p[valid_mask] = np.exp2(A[valid_mask]) / sqrt_total

    epsilon = 1e-12
    if np.any(valid_mask):
        p_valid = p[valid_mask]
        p_valid = np.clip(p_valid, epsilon, 1 - epsilon)
        p[valid_mask] = p_valid

    log_factor = math.log(2.0)
    denom = (contrast_total + reference_total) * p
    with np.errstate(divide="ignore", invalid="ignore"):
        variance = 4.0 * (1.0 - p) / (denom * (log_factor ** 2))
    sd = np.sqrt(variance)

    with np.errstate(divide="ignore", invalid="ignore"):
        mean = (np.log(contrast_total * p) - np.log(reference_total * p)) / log_factor

    z_scores = np.full_like(M, np.nan, dtype=float)
    valid_z = valid_mask & np.isfinite(sd) & (sd > 0)
    z_scores[valid_z] = (M[valid_z] - mean[valid_z]) / sd[valid_z]

    pvals = np.ones_like(M, dtype=float)
    finite_z = np.isfinite(z_scores)
    pvals[finite_z] = 2.0 * stats.norm.sf(np.abs(z_scores[finite_z]))

    log2fc_output = np.log2((c1 + 0.5) / (c2 + 0.5))

    res_df = pd.DataFrame(
        {
            "Peak": counts.index,
            "log2FC": log2fc_output,
            "A": A,
            "M": M,
            "pvalue": pvals,
            "log2FC_shrunk": log2fc_output,
        }
    ).set_index("Peak")
    res_df["padj"] = benjamini_hochberg(res_df["pvalue"]).fillna(1.0)
    res_df["method"] = "mars"
    return res_df


# ---------------------------------------------------------------------------
# Differential orchestration helpers
# ---------------------------------------------------------------------------


def call_differential_analysis(
    counts: pd.DataFrame, conditions: pd.Series, library_sizes: pd.Series
) -> pd.DataFrame:
    if counts.empty:
        raise ValueError("Counts matrix is empty; cannot perform differential analysis")
    if conditions.nunique() != 2:
        raise ValueError("Differential analysis requires exactly two experimental conditions")

    replicates_per_condition = conditions.value_counts().min()
    if replicates_per_condition >= 2:
        logging.info("Detected replicates per condition; using DESeq2 analysis via PyDESeq2")
        return pydeseq2_differential(counts, conditions)

    logging.info("No replicates detected; using MARS analysis")
    return mars_differential(counts, conditions, library_sizes)


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
    except Exception as exc:
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


def run_pipeline(
    args: SimpleNamespace,
    *,
    samples: Optional[List[SampleEntry]] = None,
    metadata_path: Optional[Path] = None,
) -> None:
    if samples is None:
        metadata_value = getattr(args, "metadata", None)
        if metadata_value is None:
            raise ValueError("Metadata file must be provided when samples are not supplied explicitly")
        metadata_path = Path(metadata_value)
        samples = load_samples(metadata_path)
    else:
        for sample in samples:
            sample.ensure_paths()

    conditions = pd.Series({s.sample: s.condition for s in samples})

    consensus_arg = getattr(args, "consensus_peaks", None)
    consensus_path = Path(consensus_arg) if consensus_arg else None

    required_cmds = ["multiBamSummary", "samtools"]
    needs_peak_calling = consensus_path is None and any(sample.peaks is None for sample in samples)
    if needs_peak_calling:
        required_cmds.append(get_macs_command())
    ensure_commands(required_cmds)

    samtools_path = shutil.which("samtools")
    if samtools_path is None:
        raise RuntimeError("samtools not found on PATH after validation")

    for sample in samples:
        sample.is_paired = detect_paired_end_bam(sample.bam, samtools_path, args.threads)

    library_sizes = compute_library_sizes(samples, samtools_path, args.threads)

    results_dir = ensure_directory(Path(args.output_dir))

    prior_manifest_value = getattr(args, "prior_manifest", None)
    manifest_data: Dict[str, object] = {}
    manifest_dir: Optional[Path] = None
    if prior_manifest_value:
        manifest_path = Path(prior_manifest_value)
        manifest_data = load_prior_manifest(manifest_path)
        manifest_dir = manifest_path.parent

    def _manifest_lookup(*keys: str) -> Optional[object]:
        for key in keys:
            if key in manifest_data:
                return manifest_data[key]
        return None

    def _normalise_path(value: object) -> Optional[Path]:
        if value is None:
            return None
        text = str(value).strip()
        if not text or text in {"-", "None", "none", "NA", "na"}:
            return None
        path = Path(text).expanduser()
        if manifest_dir and not path.is_absolute():
            path = (manifest_dir / path).expanduser()
        return path

    prior_bed_path = _normalise_path(getattr(args, "prior_bed", None))
    if prior_bed_path is None:
        prior_bed_path = _normalise_path(_manifest_lookup("prior_bed", "bed", "peaks"))

    prior_bigwig_path = _normalise_path(getattr(args, "prior_bigwig", None))
    if prior_bigwig_path is None:
        prior_bigwig_path = _normalise_path(_manifest_lookup("prior_bigwig", "bigwig", "bw"))

    prior_stats_path = _normalise_path(_manifest_lookup("prior_stats", "stats", "shape", "shape_stats"))

    manifest_weight = _manifest_lookup("prior_weight", "weight")
    if manifest_weight is not None:
        try:
            prior_weight = float(manifest_weight)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"Prior manifest weight must be numeric, received {manifest_weight!r}"
            ) from exc
    else:
        prior_weight = float(getattr(args, "prior_weight", 0.3))
    if prior_weight < 0:
        raise ValueError("--prior-weight must be non-negative")

    prior_registry = PriorRegistry(
        prior_bed=prior_bed_path,
        prior_bigwig=prior_bigwig_path,
        prior_stats=prior_stats_path,
        weight=prior_weight,
    )
    prior_registry.load_prior_distributions()
    if prior_registry.enabled:
        logging.info(
            "Prior integration enabled (weight=%.2f, bed=%s, bigwig=%s, stats=%s)",
            prior_registry.weight,
            prior_bed_path,
            prior_bigwig_path,
            prior_stats_path,
        )

    consensus: pr.PyRanges
    consensus_bed = results_dir / "consensus_peaks.bed"
    consensus_metadata: Dict[str, Optional[str]] = {
        "source": "generated" if consensus_path is None else "provided",
        "input": str(consensus_path) if consensus_path else None,
        "path": str(consensus_bed),
    }

    if consensus_path is None:
        macs2_params = {
            "genome": args.macs2_genome,
            "qvalue": args.macs2_qvalue,
            "extra": args.macs2_extra,
        }
        peak_ranges = load_all_peaks(
            samples,
            peak_extension=args.peak_extension,
            default_peak_type=args.peak_type,
            macs2_params=macs2_params,
            peak_output_dir=Path(args.peak_dir),
            prior_registry=prior_registry,
        )

        consensus = build_consensus(peak_ranges, min_overlap=args.min_overlap)
        consensus.df[["Chromosome", "Start", "End", "Name"]].to_csv(
            consensus_bed, sep="\t", header=False, index=False
        )
    else:
        logging.info("Using provided consensus peaks: %s", consensus_path)
        consensus = load_consensus_bed(consensus_path)
        ensure_directory(consensus_bed.parent)
        if consensus_path.resolve() != consensus_bed.resolve():
            shutil.copyfile(consensus_path, consensus_bed)

    if prior_registry.enabled:
        prior_registry.record_consensus(consensus.df)
        prior_registry.write_peak_tables(results_dir)
        summary = prior_registry.summarize()
        consensus_metadata["prior_overlap_fraction"] = summary.get("consensus_peaks_fraction_overlap")
        consensus_metadata["prior_overlap_count"] = summary.get("consensus_peaks_overlap")
        consensus_metadata["prior_average_weight"] = summary.get("consensus_average_weight")

    if len(consensus) == 0:
        raise ValueError("Consensus peak set is empty; cannot proceed with counting")

    counts_tsv = run_multibamsummary(consensus_bed, samples, results_dir / "counts", threads=args.threads)
    raw_counts = pd.read_csv(counts_tsv, sep="\t")

    def _normalise_header(value: str) -> str:
        cleaned = value.strip()
        cleaned = cleaned.lstrip("#")
        cleaned = cleaned.strip("'\"")
        return cleaned or value

    raw_counts.rename(columns={col: _normalise_header(col) for col in raw_counts.columns}, inplace=True)

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

    if prior_registry.enabled:
        consensus_prior_weights = prior_registry.get_consensus_weights(counts_df.index)
    else:
        consensus_prior_weights = pd.Series(0.0, index=counts_df.index, dtype=float)

    diff_res = call_differential_analysis(counts_df, conditions, library_sizes)

    prior_adjusted_path: Optional[Path] = None
    prior_adjusted: Optional[pd.DataFrame] = None
    if not diff_res.empty:
        prior_adjusted = diff_res.copy()
        if "PriorOverlap" not in prior_adjusted.columns:
            prior_adjusted["PriorOverlap"] = consensus_prior_weights > 0
        if "IntZ" not in prior_adjusted.columns:
            prior_adjusted["IntZ"] = 0.0
        if "ShapeZ" not in prior_adjusted.columns:
            prior_adjusted["ShapeZ"] = 0.0

        weights = compute_prior_weights(prior_adjusted)
        p_weighted, q_weighted, sig = weighted_bh(
            prior_adjusted["pvalue"].to_numpy(), weights, alpha=0.05
        )

        prior_adjusted["prior_weight"] = weights
        prior_adjusted["p_weighted"] = p_weighted
        prior_adjusted["q_weighted"] = q_weighted
        prior_adjusted["significant"] = sig

        diff_res["prior_weight"] = prior_adjusted["prior_weight"]
        diff_res["p_weighted"] = prior_adjusted["p_weighted"]
        diff_res["q_weighted"] = prior_adjusted["q_weighted"]
        diff_res["significant"] = prior_adjusted["significant"]

    diff_path = results_dir / "differential_results.tsv"
    diff_res.to_csv(diff_path, sep="\t")

    if prior_registry.enabled and prior_adjusted is not None:
        prior_adjusted_path = results_dir / "differential_results_prior.tsv"
        prior_adjusted.to_csv(prior_adjusted_path, sep="\t")

    plot_dir = results_dir / "plots"
    plot_paths = generate_differential_plots(diff_res, counts_df, plot_dir)

    if prior_registry.enabled:
        observed_metrics = {
            "width": (consensus.df["End"] - consensus.df["Start"]).astype(float),
            "intensity": counts_df.mean(axis=1),
        }
        prior_plot = prior_registry.plot_prior_vs_observed(
            observed_metrics, plot_dir / "prior_vs_observed.png"
        )
        if prior_plot is not None:
            plot_paths["prior_vs_observed"] = prior_plot
        prior_registry.save_distributions(results_dir / "metadata" / "prior_distributions.json")

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
                "library_size": float(library_sizes.loc[sample.sample]),
            }
        )

    args_dict = {key: value for key, value in vars(args).items() if key != "samples"}
    metadata = {
        "timestamp": datetime.utcnow().isoformat(),
        "args": args_dict,
        "metadata_sheet": str(metadata_path) if metadata_path else None,
        "samples": sample_metadata,
        "counts_matrix": str(counts_tsv),
        "differential_results": str(diff_path),
        "differential_results_prior": str(prior_adjusted_path) if prior_adjusted_path else None,
        "plots": {key: str(path) for key, path in plot_paths.items()},
        "annotation": str(annotation_path) if annotation_path else None,
        "enrichr": str(enrichr_path) if enrichr_path else None,
        "library_sizes": {name: float(value) for name, value in library_sizes.to_dict().items()},
        "consensus": consensus_metadata,
        "prior": prior_registry.metadata_entry(),
    }
    save_metadata(metadata, results_dir / "metadata.json")


# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------

def _build_condition_samples(
    condition: str,
    bam_files: Sequence[str],
    peak_files: Optional[Sequence[str]],
    used_names: Set[str],
) -> List[SampleEntry]:
    if not bam_files:
        raise ValueError(f"No BAM files supplied for condition {condition}")

    if peak_files is not None and len(peak_files) not in {0, len(bam_files)}:
        raise ValueError(
            f"Number of peak files for condition {condition} must match the number of BAM files"
        )

    entries: List[SampleEntry] = []
    for idx, bam in enumerate(bam_files):
        bam_path = Path(bam)
        base_name = bam_path.stem or f"{condition}_rep{idx + 1}"
        candidate = base_name
        suffix = 1
        while candidate in used_names:
            suffix += 1
            candidate = f"{base_name}_{suffix}"
        used_names.add(candidate)

        peak_path = None
        if peak_files:
            peak_value = peak_files[idx]
            if peak_value not in {"", "-", "None", "none", "NA", "na"}:
                peak_path = Path(peak_value)

        entries.append(
            SampleEntry(
                sample=candidate,
                condition=condition,
                bam=bam_path,
                peaks=peak_path,
                peak_type="auto",
            )
        )

    return entries


def configure_logging(level: str = "INFO"):
    logging.basicConfig(
        level=level,
        format="%(message)s",
        datefmt="[%X]",
        handlers=[RichHandler(rich_tracebacks=True, markup=True)]
    )


@app.command(name="run", help="Run the pipeline using a sample sheet (TSV/CSV).")
def run_command(
    metadata: Annotated[Path, typer.Argument(help="Path to metadata TSV/CSV file.", exists=True)],
    output_dir: Annotated[Path, typer.Option(help="Output directory")] = Path("results"),
    peak_dir: Annotated[Path, typer.Option(help="Directory for peak calls")] = Path("peaks"),
    consensus_peaks: Annotated[Optional[Path], typer.Option(help="Use existing consensus BED")] = None,
    peak_type: Annotated[str, typer.Option(help="Default peak type for MACS2 (narrow or broad)")] = "narrow",
    peak_extension: Annotated[int, typer.Option(help="Extension for narrow peaks (bp)")] = 250,
    min_overlap: Annotated[int, typer.Option(help="Minimum samples required for consensus")] = 2,
    macs2_genome: Annotated[str, typer.Option(help="MACS2 genome size (hs, mm, etc)")] = "hs",
    macs2_qvalue: Annotated[float, typer.Option(help="MACS2 q-value cutoff")] = 0.01,
    macs2_extra: Annotated[Optional[List[str]], typer.Option(help="Extra arguments for MACS2")] = None,
    prior_bed: Annotated[Optional[Path], typer.Option(help="BED file describing prior peaks")] = None,
    prior_bigwig: Annotated[Optional[Path], typer.Option(help="bigWig for prior signal")] = None,
    prior_manifest: Annotated[Optional[Path], typer.Option(help="Manifest JSON for priors")] = None,
    prior_weight: Annotated[float, typer.Option(help="Prior regularization weight")] = 0.3,
    threads: Annotated[int, typer.Option(help="Number of threads")] = 16,
    gtf: Annotated[Optional[Path], typer.Option(help="GTF file for annotation")] = None,
    enrichr: Annotated[bool, typer.Option(help="Run Enrichr analysis")] = False,
    enrichr_top: Annotated[int, typer.Option(help="Top peaks for enrichment")] = 200,
    log_level: Annotated[str, typer.Option(help="Logging level")] = "INFO",
):
    configure_logging(log_level)
    args = SimpleNamespace(
        metadata=metadata,
        output_dir=output_dir,
        peak_dir=peak_dir,
        consensus_peaks=consensus_peaks,
        peak_type=peak_type,
        peak_extension=peak_extension,
        min_overlap=min_overlap,
        macs2_genome=macs2_genome,
        macs2_qvalue=macs2_qvalue,
        macs2_extra=macs2_extra or [],
        prior_bed=prior_bed,
        prior_bigwig=prior_bigwig,
        prior_manifest=prior_manifest,
        prior_weight=prior_weight,
        threads=threads,
        gtf=gtf,
        enrichr=enrichr,
        enrichr_top=enrichr_top,
        log_level=log_level,
    )
    run_pipeline(args, metadata_path=metadata)


@app.command(name="diff", help="Run differential analysis with explicit inputs (no sample sheet).")
def diff_command(
    condition_a: Annotated[str, typer.Option(help="Name of condition A")],
    a_bams: Annotated[List[str], typer.Option(help="BAM files for condition A")],
    condition_b: Annotated[str, typer.Option(help="Name of condition B")],
    b_bams: Annotated[List[str], typer.Option(help="BAM files for condition B")],
    a_peaks: Annotated[Optional[List[str]], typer.Option(help="Peak files for condition A")] = None,
    b_peaks: Annotated[Optional[List[str]], typer.Option(help="Peak files for condition B")] = None,
    output_dir: Annotated[Path, typer.Option(help="Output directory")] = Path("results"),
    peak_dir: Annotated[Path, typer.Option(help="Directory for peak calls")] = Path("peaks"),
    consensus_peaks: Annotated[Optional[Path], typer.Option(help="Use existing consensus BED")] = None,
    peak_type: Annotated[str, typer.Option(help="Default peak type for MACS2 (narrow or broad)")] = "narrow",
    peak_extension: Annotated[int, typer.Option(help="Extension for narrow peaks (bp)")] = 250,
    min_overlap: Annotated[int, typer.Option(help="Minimum samples required for consensus")] = 2,
    macs2_genome: Annotated[str, typer.Option(help="MACS2 genome size (hs, mm, etc)")] = "hs",
    macs2_qvalue: Annotated[float, typer.Option(help="MACS2 q-value cutoff")] = 0.01,
    macs2_extra: Annotated[Optional[List[str]], typer.Option(help="Extra arguments for MACS2")] = None,
    prior_bed: Annotated[Optional[Path], typer.Option(help="BED file describing prior peaks")] = None,
    prior_bigwig: Annotated[Optional[Path], typer.Option(help="bigWig for prior signal")] = None,
    prior_manifest: Annotated[Optional[Path], typer.Option(help="Manifest JSON for priors")] = None,
    prior_weight: Annotated[float, typer.Option(help="Prior regularization weight")] = 0.3,
    threads: Annotated[int, typer.Option(help="Number of threads")] = 16,
    gtf: Annotated[Optional[Path], typer.Option(help="GTF file for annotation")] = None,
    enrichr: Annotated[bool, typer.Option(help="Run Enrichr analysis")] = False,
    enrichr_top: Annotated[int, typer.Option(help="Top peaks for enrichment")] = 200,
    log_level: Annotated[str, typer.Option(help="Logging level")] = "INFO",
):
    configure_logging(log_level)
    # Construct a mock args object for build_runmode_samples
    args = SimpleNamespace(
        condition_a=condition_a,
        a_bams=a_bams,
        a_peaks=a_peaks,
        condition_b=condition_b,
        b_bams=b_bams,
        b_peaks=b_peaks,
        output_dir=output_dir,
        peak_dir=peak_dir,
        consensus_peaks=consensus_peaks,
        peak_type=peak_type,
        peak_extension=peak_extension,
        min_overlap=min_overlap,
        macs2_genome=macs2_genome,
        macs2_qvalue=macs2_qvalue,
        macs2_extra=macs2_extra or [],
        prior_bed=prior_bed,
        prior_bigwig=prior_bigwig,
        prior_manifest=prior_manifest,
        prior_weight=prior_weight,
        threads=threads,
        gtf=gtf,
        enrichr=enrichr,
        enrichr_top=enrichr_top,
        log_level=log_level,
    )
    samples = build_runmode_samples(args)
    run_pipeline(args, samples=samples)


@app.command(name="shape", help="Analyze peak shapes.")
def shape_command(
    bed: Annotated[str, typer.Option(help="BED file with regions")],
    bigwig_a: Annotated[Optional[str], typer.Option(help="Path to sample A bigWig")] = None,
    bam_a: Annotated[Optional[str], typer.Option(help="Path to sample A BAM")] = None,
    bigwig_b: Annotated[Optional[str], typer.Option(help="Path to sample B bigWig")] = None,
    bam_b: Annotated[Optional[str], typer.Option(help="Path to sample B BAM")] = None,
    core: Annotated[int, typer.Option(help="Half-width of core window (bp)")] = 500,
    flank: Annotated[Tuple[int, int], typer.Option(help="Flank window range (bp)")] = (1000, 3000),
    fwhm_threshold: Annotated[float, typer.Option(help="Fraction of max for FWHM")] = 0.5,
    threads: Annotated[int, typer.Option(help="Number of threads")] = 1,
    out: Annotated[str, typer.Option(help="Output directory")] = "results/shape",
    prior_shape: Annotated[Optional[str], typer.Option(help="Prior shape stats file")] = None,
    prior_weight: Annotated[float, typer.Option(help="Prior weight")] = 0.3,
    bamcoverage_bin_size: Annotated[int, typer.Option(help="Bin size for bamCoverage")] = 10,
    bamcoverage_normalization: Annotated[str, typer.Option(help="Normalization for bamCoverage")] = "RPKM",
    bamcoverage_extra_args: Annotated[str, typer.Option(help="Extra args for bamCoverage")] = "",
    log_level: Annotated[str, typer.Option(help="Logging level")] = "INFO",
):
    configure_logging(log_level)
    args = SimpleNamespace(
        bed=bed,
        bigwig_a=bigwig_a,
        bam_a=bam_a,
        bigwig_b=bigwig_b,
        bam_b=bam_b,
        core=core,
        flank=list(flank),
        fwhm_threshold=fwhm_threshold,
        threads=threads,
        out=out,
        prior_shape=prior_shape,
        prior_weight=prior_weight,
        bamcoverage_bin_size=bamcoverage_bin_size,
        bamcoverage_normalization=bamcoverage_normalization,
        bamcoverage_extra_args=bamcoverage_extra_args,
        log_level=log_level
    )
    peak_shape.run_peak_shape(args)

# ---------------------------------------------------------------------------
# Entrypoint
# ---------------------------------------------------------------------------


def main() -> None:
    ensure_python_version()
    try:
        app()
    except Exception as exc:
        logging.error("Pipeline failed: %s", exc)
        sys.exit(1)


if __name__ == "__main__":
    main()
