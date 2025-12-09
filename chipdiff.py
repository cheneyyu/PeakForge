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

    peaks    Path to an existing peak file (narrowPeak or broadPeak)
    peak_type    One of {auto, narrow, broad}.  ``auto`` (default)
                 attempts to infer the peak type from the file name.

Example TSV sample sheet::

    sample  condition   bam                 peaks               peak_type
    S1      treated     data/S1.bam         data/S1_peaks.bed   narrow
    S2      treated     data/S2.bam         -                   -
    C1      control     data/C1.bam         data/C1_broad.bed   broad

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
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd
import pyranges as pr
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from scipy import stats

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

import peak_shape
from io_utils import ensure_integer_columns, read_bed_frame
from prior_utils import PriorRegistry, load_prior_manifest


def _detect_macs_command() -> str:
    """Resolve the available MACS executable.

    Preference is given to ``macs2`` when present; otherwise ``macs3`` is used.
    A runtime error is raised when neither executable is found on ``PATH``.
    """

    for candidate in ("macs2", "macs3"):
        if shutil.which(candidate):
            return candidate
    raise RuntimeError(
        "Missing required command(s): macs2, macs3. Install MACS via 'pip install macs3'"
    )


MACS_COMMAND: Optional[str] = None
"""Name of the resolved MACS executable (``macs2`` preferred, ``macs3`` fallback)."""



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
    is_paired: Optional[bool] = None
    parent_sample: Optional[str] = None

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
    """Return the available MACS executable, preferring ``macs2``.

    The resolved command is cached for subsequent calls.
    """

    global MACS_COMMAND
    if MACS_COMMAND is None:
        MACS_COMMAND = _detect_macs_command()
    return MACS_COMMAND


def ensure_python_version(min_version: tuple[int, int] = (3, 10)) -> None:
    """Guard against unsupported Python interpreters."""

    if sys.version_info < min_version:
        formatted = ".".join(str(part) for part in min_version)
        raise RuntimeError(
            f"PeakForge requires Python {formatted} or newer; detected {sys.version.split()[0]}"
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


def _maybe_threaded_cmd(cmd: List[str], threads: int) -> List[str]:
    if threads <= 1:
        return cmd

    threaded = list(cmd)
    # Insert thread flag immediately after the subcommand (e.g. 'view').
    if len(threaded) >= 2:
        threaded.insert(2, "-@")
        threaded.insert(3, str(threads))
    return threaded


def _bam_index_candidates(bam: Path) -> List[Path]:
    candidates: List[Path] = []
    candidates.append(Path(f"{bam}.bai"))
    if bam.suffix:
        candidates.append(bam.with_suffix(".bai"))
    # Remove duplicates while preserving order
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
    """Return ``True`` if the BAM contains paired-end reads."""

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
    except ValueError as exc:  # pragma: no cover - defensive
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
    """Read a delimited table inferring delimiter automatically."""

    try:
        df = pd.read_csv(path, sep=None, engine="python")
    except Exception as exc:  # pragma: no cover - passthrough error
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


def _sequential_split_bam(
    sample: SampleEntry,
    output_dir: Path,
    samtools_path: str,
    replicates: int,
    threads: int = 1,
) -> List[Path]:
    if replicates < 1:
        raise ValueError("replicates must be a positive integer")

    ensure_directory(output_dir)

    try:
        total_reads = int(
            subprocess.check_output(
                _maybe_threaded_cmd([samtools_path, "view", "-c", str(sample.bam)], threads),
                text=True,
            ).strip()
        )
    except subprocess.CalledProcessError as exc:  # pragma: no cover - external command
        raise RuntimeError(f"Failed to count reads for {sample.bam}: {exc}") from exc
    except ValueError as exc:
        raise RuntimeError(
            f"Non-numeric read count returned by samtools for {sample.bam}"
        ) from exc

    if total_reads <= 0:
        raise ValueError(f"No alignments available in BAM {sample.bam}")

    chunk_size = math.ceil(total_reads / replicates)
    header_proc = subprocess.run(
        _maybe_threaded_cmd([samtools_path, "view", "-H", str(sample.bam)], threads),
        capture_output=True,
        text=True,
        check=True,
    )
    header_lines = header_proc.stdout.splitlines(keepends=True)

    procs: List[subprocess.Popen[str]] = []
    outputs: List[Path] = []
    for idx in range(replicates):
        out_path = output_dir / f"{sample.sample}_pseudo{idx + 1}.bam"
        proc = subprocess.Popen(
            _maybe_threaded_cmd(
                [samtools_path, "view", "-b", "-o", str(out_path), "-"], threads
            ),
            stdin=subprocess.PIPE,
            text=True,
        )
        procs.append(proc)
        outputs.append(out_path)

    stream_proc = subprocess.Popen(
        _maybe_threaded_cmd([samtools_path, "view", str(sample.bam)], threads),
        stdout=subprocess.PIPE,
        text=True,
    )

    processed = 0
    target_idx = 0
    try:
        for line in header_lines:
            for proc in procs:
                if proc.stdin:
                    proc.stdin.write(line)

        if stream_proc.stdout is None:  # pragma: no cover - defensive guard
            raise RuntimeError("samtools stream missing stdout handle")

        for line in stream_proc.stdout:
            target_proc = procs[target_idx]
            if target_proc.stdin:
                target_proc.stdin.write(line)
            processed += 1
            if processed % chunk_size == 0 and target_idx < replicates - 1:
                target_idx += 1
    finally:
        if stream_proc.stdout:
            stream_proc.stdout.close()

    stream_proc.wait()
    for proc in procs:
        if proc.stdin:
            proc.stdin.close()
    return_codes = [proc.wait() for proc in procs]

    if stream_proc.returncode not in (0, None):
        raise RuntimeError(f"samtools view failed while splitting {sample.bam}")
    for idx, code in enumerate(return_codes):
        if code != 0:
            raise RuntimeError(
                f"samtools view failed for pseudo-replicate {idx + 1} of {sample.bam}"
            )

    return outputs


def generate_pseudo_replicates(
    samples: Sequence[SampleEntry],
    *,
    output_root: Path,
    samtools_path: str,
    threads: int = 1,
    replicates: int = 3,
) -> Tuple[List[SampleEntry], Dict[str, List[str]]]:
    """Split each BAM sequentially into pseudo-replicates.

    Returns the expanded sample list and a mapping from original samples to
    their derived pseudo-replicate names.
    """

    ensure_directory(output_root)
    expanded: List[SampleEntry] = []
    mapping: Dict[str, List[str]] = {}

    for sample in samples:
        logging.info(
            "Generating %d pseudo-replicates for %s", replicates, sample.sample
        )
        sample_out_dir = output_root / sample.sample
        outputs = _sequential_split_bam(
            sample,
            sample_out_dir,
            samtools_path,
            replicates,
            threads,
        )
        mapping[sample.sample] = []
        for idx, path in enumerate(outputs, start=1):
            name = f"{sample.sample}_pr{idx}"
            mapping[sample.sample].append(name)
            expanded.append(
                SampleEntry(
                    sample=name,
                    condition=sample.condition,
                    bam=path,
                    peaks=sample.peaks,
                    peak_type=sample.peak_type,
                    parent_sample=sample.sample,
                )
            )

    return expanded, mapping


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
    """Load peaks into a :class:`pyranges.PyRanges` object."""

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
    """Call MACS2 for a sample and return the resulting peak file path."""

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
    """Ensure every sample has peak calls and return PyRanges per sample."""

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
    """Build consensus peaks across samples with minimum overlap criteria."""

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
    """Load an existing consensus BED file into a ``PyRanges`` object."""

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
        # Guarantee uniqueness in case the BED supplies duplicates
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
    """Implement the MARS method for designs without replicates."""

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


def run_pipeline(
    args: argparse.Namespace,
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

    results_dir = ensure_directory(Path(args.output_dir))
    pseudorep_count = 3

    single_replicate_mode = getattr(args, "single_replicate_mode", "mars")
    conditions = pd.Series({s.sample: s.condition for s in samples})
    replicates_per_condition = conditions.value_counts().min()

    ensure_commands(["samtools"])
    samtools_path = shutil.which("samtools")
    if samtools_path is None:  # pragma: no cover - defensive (ensure_commands guards)
        raise RuntimeError("samtools not found on PATH after validation")

    pseudorep_mapping: Dict[str, List[str]] = {}
    if replicates_per_condition < 2 and single_replicate_mode == "pseudo":
        if conditions.nunique() != 2:
            raise ValueError(
                "Pseudo-replicate mode requires exactly two experimental conditions"
            )

        pseudo_root = results_dir / "pseudoreplicates"
        samples, pseudorep_mapping = generate_pseudo_replicates(
            samples,
            output_root=pseudo_root,
            samtools_path=samtools_path,
            threads=args.threads,
            replicates=pseudorep_count,
        )
        for sample in samples:
            sample.ensure_paths()
        conditions = pd.Series({s.sample: s.condition for s in samples})
        replicates_per_condition = conditions.value_counts().min()

    consensus_arg = getattr(args, "consensus_peaks", None)
    consensus_path = Path(consensus_arg) if consensus_arg else None

    required_cmds = ["multiBamSummary", "samtools"]
    needs_peak_calling = consensus_path is None and any(sample.peaks is None for sample in samples)
    if needs_peak_calling:
        required_cmds.append(get_macs_command())
    ensure_commands(required_cmds)

    for sample in samples:
        sample.is_paired = detect_paired_end_bam(sample.bam, samtools_path, args.threads)

    library_sizes = compute_library_sizes(samples, samtools_path, args.threads)

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

    if prior_registry.enabled:
        consensus_prior_weights = prior_registry.get_consensus_weights(counts_df.index)
    else:
        consensus_prior_weights = pd.Series(0.0, index=counts_df.index, dtype=float)

    diff_res = call_differential_analysis(counts_df, conditions, library_sizes)

    diff_res["prior_weight"] = consensus_prior_weights
    diff_res["prior_overlap"] = consensus_prior_weights > 0

    diff_path = results_dir / "differential_results.tsv"
    diff_res.to_csv(diff_path, sep="\t")

    prior_adjusted_path: Optional[Path] = None
    if prior_registry.enabled and not diff_res.empty:
        overlap_factor = consensus_prior_weights.clip(0, 1)
        shrink_factor = 1.0 - prior_registry.weight * (1.0 - overlap_factor)
        penalty_factor = 1.0 + prior_registry.weight * (1.0 - overlap_factor)
        prior_adjusted = diff_res.copy()
        prior_adjusted["prior_weight"] = overlap_factor
        prior_adjusted["prior_overlap"] = overlap_factor > 0
        prior_adjusted["log2FC_prior"] = diff_res["log2FC"] * shrink_factor
        if "log2FC_shrunk" in diff_res.columns:
            prior_adjusted["log2FC_shrunk_prior"] = diff_res["log2FC_shrunk"] * shrink_factor
        if "lfcSE" in diff_res.columns:
            prior_adjusted["lfcSE_prior"] = diff_res["lfcSE"] * penalty_factor
        if "pvalue" in diff_res.columns:
            prior_adjusted["pvalue_prior"] = np.minimum(1.0, diff_res["pvalue"] * penalty_factor)
            prior_adjusted["padj_prior"] = benjamini_hochberg(prior_adjusted["pvalue_prior"]).fillna(1.0)
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
                "parent_sample": sample.parent_sample,
            }
        )

    args_dict = {key: value for key, value in vars(args).items() if key != "samples"}
    pseudorep_metadata: Optional[Dict[str, object]] = None
    if pseudorep_mapping:
        pseudorep_metadata = {
            "mode": single_replicate_mode,
            "replicates_per_sample": pseudorep_count,
            "mapping": pseudorep_mapping,
        }
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
        "pseudoreplicates": pseudorep_metadata,
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


def build_runmode_samples(args: argparse.Namespace) -> List[SampleEntry]:
    a_peaks = args.a_peaks or []
    b_peaks = args.b_peaks or []

    used_names: Set[str] = set()
    samples = _build_condition_samples(args.condition_a, args.a_bams, a_peaks, used_names)
    samples.extend(_build_condition_samples(args.condition_b, args.b_bams, b_peaks, used_names))
    return samples


def add_common_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--output-dir", default="results", help="Output directory")
    parser.add_argument("--peak-dir", default="peaks", help="Directory for peak calls")
    parser.add_argument(
        "--consensus-peaks",
        help="Use an existing consensus BED instead of building peaks from samples",
    )
    parser.add_argument(
        "--peak-type",
        default="narrow",
        choices=["narrow", "broad"],
        help="Default peak type when calling MACS2",
    )
    parser.add_argument(
        "--peak-extension",
        type=int,
        default=250,
        help="Extension for narrow peaks when building consensus windows (bp)",
    )
    parser.add_argument(
        "--min-overlap",
        type=int,
        default=2,
        help="Minimum samples required for consensus peak",
    )
    parser.add_argument("--macs2-genome", default="hs", help="MACS2 genome size (e.g. hs, mm, 2.7e9)")
    parser.add_argument("--macs2-qvalue", type=float, default=0.01, help="MACS2 q-value cutoff")
    parser.add_argument(
        "--macs2-extra",
        nargs=argparse.REMAINDER,
        default=[],
        help="Additional arguments for MACS2",
    )
    parser.add_argument("--prior-bed", help="BED file describing prior peaks")
    parser.add_argument("--prior-bigwig", help="bigWig of reference signal intensities for priors")
    parser.add_argument("--prior-manifest", help="Manifest (JSON or key=value) describing prior resources")
    parser.add_argument(
        "--prior-weight",
        type=float,
        default=0.3,
        help="Strength of the prior regularisation (0 disables influence)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=16,
        help="Threads for multiBamSummary (--numberOfProcessors)",
    )
    parser.add_argument(
        "--single-replicate-mode",
        choices=["mars", "pseudo"],
        default="mars",
        help=(
            "Strategy for single-replicate (1v1) designs: 'mars' keeps the"
            " existing MARS workflow while 'pseudo' splits each BAM into"
            " sequential pseudo-replicates"
        ),
    )
    parser.add_argument("--gtf", help="Optional GTF file for annotation")
    parser.add_argument("--enrichr", action="store_true", help="Run Enrichr GO Biological Process analysis")
    parser.add_argument("--enrichr-top", type=int, default=200, help="Number of top peaks for enrichment")
    parser.add_argument("--log-level", default="INFO", help="Logging level")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="peakforge",
        description="CUT&Tag / ChIP-seq differential analysis pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    tsv_parser = subparsers.add_parser(
        "tsvmode",
        help="Run the pipeline using a metadata TSV/CSV sheet",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    tsv_parser.add_argument("metadata", help="Sample sheet (TSV/CSV)")
    add_common_arguments(tsv_parser)

    run_parser = subparsers.add_parser(
        "runmode",
        help="Run the pipeline by specifying BAM/peak files directly",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    run_parser.add_argument("--condition-a", required=True, help="Reference condition label")
    run_parser.add_argument("--a-bams", nargs="+", required=True, help="BAM files for condition A")
    run_parser.add_argument(
        "--a-peaks",
        nargs="*",
        default=None,
        help="Optional peak files aligned with --a-bams",
    )
    run_parser.add_argument("--condition-b", required=True, help="Contrast condition label")
    run_parser.add_argument("--b-bams", nargs="+", required=True, help="BAM files for condition B")
    run_parser.add_argument(
        "--b-peaks",
        nargs="*",
        default=None,
        help="Optional peak files aligned with --b-bams",
    )
    add_common_arguments(run_parser)

    peakshape_parser = subparsers.add_parser(
        "peakshape",
        help="Run peak shape profiling for two bigWig tracks",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    peak_shape.add_cli_arguments(peakshape_parser)

    return parser


# ---------------------------------------------------------------------------
# Entrypoint
# ---------------------------------------------------------------------------


def main(argv: Optional[Sequence[str]] = None) -> None:
    ensure_python_version()
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level.upper(), logging.INFO),
        format="[%(asctime)s] %(levelname)s: %(message)s",
    )

    try:
        if args.command == "tsvmode":
            metadata_path = Path(args.metadata)
            samples = load_samples(metadata_path)
            run_pipeline(args, samples=samples, metadata_path=metadata_path)
        elif args.command == "runmode":
            samples = build_runmode_samples(args)
            run_pipeline(args, samples=samples, metadata_path=None)
        elif args.command == "peakshape":
            peak_shape.run_peak_shape(args)
        else:  # pragma: no cover - defensive guard
            parser.print_help()
            sys.exit(1)
    except Exception as exc:  # pragma: no cover - CLI exception reporting
        logging.error("Pipeline failed: %s", exc)
        sys.exit(1)


if __name__ == "__main__":
    main()
