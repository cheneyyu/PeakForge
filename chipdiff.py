#!/usr/bin/env python3
"""chipdiff.py

Two-group narrow-peak chromatin differential analysis and exploratory ranking.

The script expects a sample sheet describing the input BAM files and
(optional) pre-computed peak files.  The sample sheet must be a tab- or
comma-delimited text file with the following required columns:

    sample    Unique sample identifier (no spaces)
    condition    Experimental condition or group label
    bam    Path to the aligned reads in BAM format

Optional columns:

    control_bam    Path to a matched control/input BAM to be passed to
                   MACS2/MACS3 via ``-c`` during automatic peak calling.
                   This is strongly recommended for ChIP-seq and usually
                   omitted for ATAC-seq and CUT&Tag.
    peaks    Path to an existing peak file (narrowPeak or broadPeak)
    peak_type    One of {auto, narrow, broad}.  ``auto`` (default)
                 attempts to infer the peak type from the file name.
                 Use ``narrow`` for punctate enrichments (e.g., TF/CUT&Tag
                 factor peaks) and ``broad`` for diffuse domains (e.g.,
                 many histone-mark ChIP-seq datasets).

Example TSV sample sheet::

    sample  condition   bam                 control_bam         peaks               peak_type
    T1      treated     data/T1.bam         data/input.bam      data/T1_peaks.narrowPeak   narrow
    T2      treated     data/T2.bam         data/input.bam      data/T2_peaks.narrowPeak   narrow
    C1      control     data/C1.bam         data/input.bam      data/C1_peaks.narrowPeak   narrow
    C2      control     data/C2.bam         data/input.bam      data/C2_peaks.narrowPeak   narrow

For a single analysis, using a consistent peak type across samples is
recommended.

The pipeline performs the following steps:

1. Peak calling with MACS2/3 (if required, optionally using matched controls).
2. Construction of consensus peaks across samples.
3. Counting read overlaps per consensus peak using deepTools
   ``multiBamSummary``.
4. Formal differential inference with PyDESeq2 for replicated designs, or
   exploratory effect-size/MARS-derived ranking for exact 1-vs-1 designs.
5. Optional annotation against a GTF file and Enrichr enrichment via
   gseapy.
6. Plot generation (volcano, MA, correlation, heatmap) and metadata
   capture.

Dependencies: numpy, pandas, scipy, statsmodels, matplotlib, seaborn,
pyranges, gseapy, MACS2, deepTools.
"""
from __future__ import annotations

import argparse
import importlib.metadata
import json
import logging
import math
import platform
import shutil
import subprocess
import sys
import warnings
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set

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

try:
    import pysam
except ImportError:  # pragma: no cover - optional dependency
    pysam = None

import peak_shape
from io_utils import ensure_integer_columns, read_bed_frame
from motif_ranking import run_pairwise_motif_ranking


__version__ = "0.2.3"

SAM_PROPER_PAIR = 0x2
SAM_UNMAPPED = 0x4
SAM_FIRST_MATE = 0x40
SAM_SECONDARY = 0x100
SAM_QCFAIL = 0x200
SAM_DUPLICATE = 0x400
SAM_SUPPLEMENTARY = 0x800


def _detect_macs_command() -> str:
    """Resolve the available MACS executable.

    Preference is given to ``macs3`` when it can be executed successfully;
    otherwise a working ``macs2`` is used as fallback. A runtime error is
    raised when neither executable is runnable on ``PATH``.
    """

    for candidate in ("macs3", "macs2"):
        resolved = shutil.which(candidate)
        if resolved is None:
            continue
        try:
            result = subprocess.run(
                [candidate, "--version"],
                check=False,
                capture_output=True,
                text=True,
                timeout=15,
            )
        except subprocess.TimeoutExpired:
            logging.warning("Skipping unusable MACS executable %s at %s (timed out)", candidate, resolved)
            continue
        if result.returncode == 0:
            return candidate
        logging.warning(
            "Skipping unusable MACS executable %s at %s (exit %s): %s",
            candidate,
            resolved,
            result.returncode,
            (result.stderr or result.stdout).strip(),
        )
    raise RuntimeError(
        "Missing a runnable MACS executable. Install MACS via 'pip install macs3' "
        "(recommended) or install a working macs2 binary."
    )


MACS_COMMAND: Optional[str] = None
"""Name of the resolved MACS executable (``macs3`` preferred, ``macs2`` fallback)."""



# ---------------------------------------------------------------------------
# Data classes and utility helpers
# ---------------------------------------------------------------------------


@dataclass
class SampleEntry:
    """Representation of a single sample entry from the metadata sheet."""

    sample: str
    condition: str
    bam: Path
    control_bam: Optional[Path] = None
    peaks: Optional[Path] = None
    peak_type: str = "auto"
    is_paired: Optional[bool] = None

    def ensure_paths(self) -> None:
        if not self.bam.exists():
            raise FileNotFoundError(f"BAM file not found for sample {self.sample}: {self.bam}")
        if self.control_bam is not None and not self.control_bam.exists():
            raise FileNotFoundError(
                f"Control/Input BAM file not found for sample {self.sample}: {self.control_bam}"
            )
        if self.peaks is not None and not self.peaks.exists():
            raise FileNotFoundError(f"Peak file not found for sample {self.sample}: {self.peaks}")


@dataclass(frozen=True)
class CountFilterConfig:
    """Filtering and count-unit definition shared by counting and library size.

    ``fragment`` mode counts the first mate of each proper pair once and asks
    deepTools to extend that mate across the complete paired-end fragment.
    This makes interval counts and the selected library-size denominator use
    the same countable unit. Tn5 insertion-site counting is intentionally not
    implemented.
    """

    count_unit: str = "read"
    min_mapq: int = 0
    exclude_duplicates: bool = False
    proper_pairs_only: bool = False
    sam_flag_exclude: int = 0

    def validate(self) -> None:
        if self.count_unit not in {"read", "fragment"}:
            raise ValueError("count_unit must be either 'read' or 'fragment'")
        if self.min_mapq < 0:
            raise ValueError("min_mapq must be non-negative")
        if self.sam_flag_exclude < 0:
            raise ValueError("sam_flag_exclude must be non-negative")

    @property
    def effective_exclude_flags(self) -> int:
        flags = (
            SAM_UNMAPPED
            | SAM_SECONDARY
            | SAM_QCFAIL
            | SAM_SUPPLEMENTARY
            | int(self.sam_flag_exclude)
        )
        if self.exclude_duplicates:
            flags |= SAM_DUPLICATE
        return flags

    def effective_include_flags(self, *, paired_end: bool) -> int:
        if self.count_unit == "fragment":
            if not paired_end:
                raise ValueError("fragment counting requires a paired-end BAM")
            return SAM_FIRST_MATE | SAM_PROPER_PAIR
        if self.proper_pairs_only:
            if not paired_end:
                raise ValueError("proper-pairs-only filtering requires a paired-end BAM")
            return SAM_PROPER_PAIR
        return 0

    def as_metadata(self) -> Dict[str, object]:
        return {
            "count_unit": self.count_unit,
            "min_mapq": self.min_mapq,
            "exclude_duplicates": self.exclude_duplicates,
            "proper_pairs_only": self.proper_pairs_only or self.count_unit == "fragment",
            "sam_flag_exclude_requested": self.sam_flag_exclude,
            "sam_flag_exclude_effective": self.effective_exclude_flags,
            "tn5_insertion_site_counting": False,
        }


# ---------------------------------------------------------------------------
# File and command helpers
# ---------------------------------------------------------------------------


def ensure_commands(commands: Sequence[str]) -> None:
    missing = [cmd for cmd in commands if shutil.which(cmd) is None]
    if missing:
        joined = ", ".join(sorted(missing))
        raise RuntimeError(
            "Missing required command(s): "
            f"{joined}. Install MACS via 'pip install macs3' (recommended) or a working macs2; deepTools via "
            "'pip install deeptools' "
            "and samtools via 'conda install -c bioconda samtools'."
        )


def get_macs_command() -> str:
    """Return the available MACS executable, preferring ``macs3``.

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


def git_state() -> Dict[str, object]:
    """Return the source commit and dirty-state flag for result provenance."""

    repo = Path(__file__).resolve().parent
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=repo,
        check=False,
        capture_output=True,
        text=True,
    )
    status = subprocess.run(
        ["git", "status", "--porcelain"],
        cwd=repo,
        check=False,
        capture_output=True,
        text=True,
    )
    return {
        "commit": commit.stdout.strip() if commit.returncode == 0 else None,
        "dirty": bool(status.stdout.strip()) if status.returncode == 0 else None,
    }


def _command_version(command: str, *args: str) -> Optional[str]:
    resolved = shutil.which(command)
    if resolved is None:
        return None
    completed = subprocess.run(
        [resolved, *args],
        check=False,
        capture_output=True,
    )
    combined = b"\n".join(part for part in (completed.stdout, completed.stderr) if part).decode(
        "utf-8", errors="replace"
    )
    first_line = next((line.strip() for line in combined.splitlines() if line.strip()), None)
    return first_line


def software_versions() -> Dict[str, object]:
    packages = {}
    for package in (
        "deeptools",
        "gseapy",
        "matplotlib",
        "numpy",
        "pandas",
        "pydeseq2",
        "pyranges",
        "pysam",
        "scipy",
        "seaborn",
        "statsmodels",
    ):
        try:
            packages[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            packages[package] = None
    return {
        "peakforge": __version__,
        "python": platform.python_version(),
        "platform": platform.platform(),
        "python_packages": packages,
        "external_tools": {
            "samtools": _command_version("samtools", "--version"),
            "multiBamSummary": _command_version("multiBamSummary", "--version"),
            "macs3": _command_version("macs3", "--version"),
            "macs2": _command_version("macs2", "--version"),
        },
    }


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


def samtools_count_units(
    bam: Path,
    samtools_path: str,
    *,
    min_mapq: int,
    exclude_flags: int,
    include_flags: int = 0,
    threads: int = 1,
) -> int:
    """Count alignments matching an explicit samtools flag/MAPQ definition."""

    cmd = [samtools_path, "view", "-c", "-q", str(min_mapq), "-F", str(exclude_flags)]
    if include_flags:
        cmd.extend(["-f", str(include_flags)])
    if threads > 1:
        cmd.extend(["-@", str(threads)])
    cmd.append(str(bam))
    result = subprocess.run(cmd, check=False, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"samtools view failed for {bam} with exit code {result.returncode}: "
            f"{result.stderr.strip()}"
        )
    try:
        return int(result.stdout.strip() or 0)
    except ValueError as exc:
        raise RuntimeError(
            f"Unable to parse samtools count for {bam}: {result.stdout!r}"
        ) from exc


def compute_library_size_metrics(
    samples: Sequence[SampleEntry],
    samtools_path: str,
    count_config: CountFilterConfig,
    threads: int = 1,
) -> pd.DataFrame:
    """Return raw, filtered-alignment, and selected count-unit totals."""

    count_config.validate()
    records: List[Dict[str, object]] = []
    for sample in samples:
        logging.info("Estimating library size for sample %s", sample.sample)
        paired_end = bool(sample.is_paired)
        countable_include = count_config.effective_include_flags(paired_end=paired_end)
        alignment_include = (
            SAM_PROPER_PAIR
            if (count_config.proper_pairs_only or count_config.count_unit == "fragment")
            else 0
        )
        raw_mapped = bam_total_mapped_reads(sample.bam, samtools_path, threads)
        filtered_alignments = samtools_count_units(
            sample.bam,
            samtools_path,
            min_mapq=count_config.min_mapq,
            exclude_flags=count_config.effective_exclude_flags,
            include_flags=alignment_include,
            threads=threads,
        )
        filtered_units = samtools_count_units(
            sample.bam,
            samtools_path,
            min_mapq=count_config.min_mapq,
            exclude_flags=count_config.effective_exclude_flags,
            include_flags=countable_include,
            threads=threads,
        )
        if filtered_units <= 0:
            raise ValueError(
                f"No countable {count_config.count_unit} units remain after filtering for {sample.sample}"
            )
        records.append(
            {
                "sample": sample.sample,
                "paired_end": paired_end,
                "raw_mapped_alignments": raw_mapped,
                "filtered_alignments": filtered_alignments,
                "filtered_countable_units": filtered_units,
                "selected_library_size": filtered_units,
                **count_config.as_metadata(),
            }
        )
    return pd.DataFrame.from_records(records).set_index("sample")


def compute_library_sizes(
    samples: Sequence[SampleEntry],
    samtools_path: str,
    threads: int = 1,
    count_config: Optional[CountFilterConfig] = None,
) -> pd.Series:
    """Return selected library sizes under the same filters as interval counts."""

    config = count_config or CountFilterConfig()
    metrics = compute_library_size_metrics(samples, samtools_path, config, threads)
    return metrics["selected_library_size"].astype(float)


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
        control_bam = _normalise_optional_path(
            getattr(row, "control_bam", getattr(row, "input_bam", None))
        )
        peaks = _normalise_optional_path(getattr(row, "peaks", None))
        peak_type_raw = getattr(row, "peak_type", "auto")
        peak_type = str(peak_type_raw).lower() if peak_type_raw is not None else "auto"
        entry = SampleEntry(
            sample=str(sample),
            condition=str(condition),
            bam=Path(str(bam)),
            control_bam=control_bam,
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


def read_peak_file(
    path: Path,
    peak_type: str,
    peak_edge_padding: int,
    *,
    peak_coordinate_mode: str = "edge-padding",
    summit_fixed_width: int = 500,
) -> pr.PyRanges:
    """Load peaks using native/edge-padded or summit-centred coordinates.

    ``edge-padding`` preserves the historical behaviour: the requested number
    of bases is added to both native narrowPeak edges.  ``summit-fixed`` reads
    the zero-based summit offset from narrowPeak column 10 and creates an exact
    total-width interval around that summit.  Broad peaks do not define a
    narrowPeak summit and therefore cannot use ``summit-fixed``.
    """

    if peak_edge_padding < 0:
        raise ValueError("peak edge padding must be non-negative")

    coordinate_mode = peak_coordinate_mode.replace("_", "-").lower()
    if coordinate_mode not in {"edge-padding", "summit-fixed"}:
        raise ValueError(
            "peak coordinate mode must be 'edge-padding' or 'summit-fixed'"
        )

    if coordinate_mode == "summit-fixed":
        if peak_type != "narrow":
            raise ValueError("summit-fixed coordinates require narrowPeak input")
        if summit_fixed_width <= 0:
            raise ValueError("summit fixed width must be a positive integer")
        narrowpeak_columns = (
            "Chromosome",
            "Start",
            "End",
            "Name",
            "Score",
            "Strand",
            "SignalValue",
            "PValue",
            "QValue",
            "SummitOffset",
        )
        frame = read_bed_frame(
            path,
            column_names=narrowpeak_columns,
            min_columns=10,
        )
        frame = ensure_integer_columns(frame, ("Start", "End", "SummitOffset"))
        native_width = frame["End"] - frame["Start"]
        invalid = (frame["SummitOffset"] < 0) | (frame["SummitOffset"] >= native_width)
        if invalid.any():
            invalid_rows = ", ".join(str(index + 1) for index in frame.index[invalid][:5])
            raise ValueError(
                f"narrowPeak summit offset must lie within its interval in {path}; "
                f"invalid data row(s): {invalid_rows}"
            )
        summit = frame["Start"].to_numpy() + frame["SummitOffset"].to_numpy()
        start = np.maximum(summit - summit_fixed_width // 2, 0)
        frame["Start"] = start
        frame["End"] = start + summit_fixed_width
        return pr.PyRanges(frame[["Chromosome", "Start", "End"]])

    frame = read_bed_frame(path)
    frame = ensure_integer_columns(frame, ("Start", "End"))

    if peak_type == "narrow":
        start = frame["Start"].to_numpy() - peak_edge_padding
        end = frame["End"].to_numpy() + peak_edge_padding
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
    macs_format: str = "AUTO",
    macs_nomodel: bool = False,
    macs_shift: Optional[int] = None,
    macs_extsize: Optional[int] = None,
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
    if sample.control_bam is not None:
        cmd.extend(["-c", str(sample.control_bam)])
    resolved_format = macs_format.upper()
    if resolved_format == "AUTO":
        resolved_format = "BAMPE" if sample.is_paired else "BAM"
    if resolved_format not in {"BAM", "BAMPE"}:
        raise ValueError("macs_format must be AUTO, BAM, or BAMPE")
    if resolved_format == "BAMPE":
        if macs_shift not in {None, 0}:
            raise ValueError(
                "MACS paired-end BAMPE mode uses observed fragments and requires "
                "macs_shift to be 0 or omitted"
            )
        if macs_extsize is not None:
            raise ValueError(
                "MACS paired-end BAMPE mode uses observed fragment lengths; "
                "macs_extsize is only supported with read-based BAM mode"
            )
    cmd.extend(["-f", resolved_format])

    if macs_nomodel:
        cmd.append("--nomodel")
    if macs_shift is not None:
        cmd.extend(["--shift", str(macs_shift)])
    if macs_extsize is not None:
        if macs_extsize <= 0:
            raise ValueError("macs_extsize must be positive")
        cmd.extend(["--extsize", str(macs_extsize)])

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
    macs_format: str = "AUTO",
    macs_nomodel: bool = False,
    macs_shift: Optional[int] = None,
    macs_extsize: Optional[int] = None,
    macs2_extra: Optional[List[str]] = None,
) -> Path:
    """Call MACS2 for a sample and return the resulting peak file path."""

    cmd, peak_path = _macs2_command(
        sample,
        output_dir=output_dir,
        macs2_genome=macs2_genome,
        macs2_qval=macs2_qval,
        peak_type=peak_type,
        macs_format=macs_format,
        macs_nomodel=macs_nomodel,
        macs_shift=macs_shift,
        macs_extsize=macs_extsize,
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
    peak_edge_padding: int,
    peak_coordinate_mode: str,
    summit_fixed_width: int,
    default_peak_type: str,
    macs2_params: Dict[str, object],
    peak_output_dir: Path,
    command_log: Optional[List[List[str]]] = None,
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
                macs_format=str(macs2_params.get("format", "AUTO")),
                macs_nomodel=bool(macs2_params.get("nomodel", False)),
                macs_shift=macs2_params.get("shift"),
                macs_extsize=macs2_params.get("extsize"),
                macs2_extra=macs2_params.get("extra", []),
            )
            if command_log is not None:
                command_log.append(list(cmd))
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

            pr_obj = read_peak_file(
                Path(peak_path),
                peak_type,
                peak_edge_padding,
                peak_coordinate_mode=peak_coordinate_mode,
                summit_fixed_width=summit_fixed_width,
            )
            df = pr_obj.df
            df["Sample"] = sample.sample
            pr_obj = pr.PyRanges(df)
            peak_ranges[sample.sample] = pr_obj

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

        pr_obj = read_peak_file(
            Path(job.peak_path),
            job.peak_type,
            peak_edge_padding,
            peak_coordinate_mode=peak_coordinate_mode,
            summit_fixed_width=summit_fixed_width,
        )
        df = pr_obj.df
        df["Sample"] = job.sample.sample
        pr_obj = pr.PyRanges(df)
        peak_ranges[job.sample.sample] = pr_obj
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


def write_consensus_bed(consensus: pr.PyRanges, output_path: Path) -> None:
    """Write consensus intervals to a BED file."""

    ensure_directory(output_path.parent)
    consensus.df[["Chromosome", "Start", "End", "Name"]].to_csv(
        output_path,
        sep="\t",
        header=False,
        index=False,
    )


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


def load_interval_bed(path: Path) -> pr.PyRanges:
    """Load a generic BED file as intervals suitable for overlap filtering."""

    if not path.exists():
        raise FileNotFoundError(f"BED file not found: {path}")
    frame = read_bed_frame(path)
    frame = ensure_integer_columns(frame, ("Start", "End"))
    return pr.PyRanges(frame[["Chromosome", "Start", "End"]])


def filter_consensus_by_blacklist(
    consensus: pr.PyRanges,
    blacklist_path: Path,
) -> tuple[pr.PyRanges, Dict[str, object], pd.DataFrame]:
    """Drop consensus peaks with any overlap against a blacklist BED."""

    blacklist = load_interval_bed(blacklist_path)
    df = consensus.df.copy()
    if df.empty:
        summary = {
            "applied": True,
            "blacklist_bed": str(blacklist_path),
            "total_peaks": 0,
            "kept_peaks": 0,
            "removed_peaks": 0,
            "overlapped_peaks": 0,
            "max_blacklist_overlaps": 0,
        }
        empty_report = pd.DataFrame(
            columns=[
                "Chromosome",
                "Start",
                "End",
                "Name",
                "Support",
                "blacklist_overlaps",
                "keep",
                "filter_reason",
            ]
        )
        return consensus, summary, empty_report

    if len(blacklist) == 0:
        report_df = df[["Chromosome", "Start", "End", "Name", "Support"]].copy()
        report_df["blacklist_overlaps"] = 0
        report_df["keep"] = True
        report_df["filter_reason"] = ""
        summary = {
            "applied": True,
            "blacklist_bed": str(blacklist_path),
            "total_peaks": int(len(report_df)),
            "kept_peaks": int(len(report_df)),
            "removed_peaks": 0,
            "overlapped_peaks": 0,
            "max_blacklist_overlaps": 0,
        }
        return consensus, summary, report_df

    overlap_df = consensus.count_overlaps(blacklist).df.copy()
    overlap_df["blacklist_overlaps"] = (
        pd.to_numeric(overlap_df.get("NumberOverlaps", 0), errors="coerce")
        .fillna(0)
        .astype(int)
    )
    overlap_df["keep"] = overlap_df["blacklist_overlaps"] == 0
    overlap_df["filter_reason"] = np.where(overlap_df["keep"], "", "blacklist")

    kept_df = overlap_df.loc[
        overlap_df["keep"],
        ["Chromosome", "Start", "End", "Name", "Support"],
    ].copy()
    filtered_consensus = pr.PyRanges(kept_df)

    overlapped = int((overlap_df["blacklist_overlaps"] > 0).sum())
    max_overlaps = int(overlap_df["blacklist_overlaps"].max()) if len(overlap_df) else 0
    summary = {
        "applied": True,
        "blacklist_bed": str(blacklist_path),
        "total_peaks": int(len(overlap_df)),
        "kept_peaks": int(len(kept_df)),
        "removed_peaks": int(len(overlap_df) - len(kept_df)),
        "overlapped_peaks": overlapped,
        "max_blacklist_overlaps": max_overlaps,
    }
    report_df = overlap_df[
        [
            "Chromosome",
            "Start",
            "End",
            "Name",
            "Support",
            "blacklist_overlaps",
            "keep",
            "filter_reason",
        ]
    ].copy()
    return filtered_consensus, summary, report_df


def validate_fraction_threshold(name: str, value: Optional[float]) -> Optional[float]:
    """Validate a fraction threshold constrained to the [0, 1] interval."""

    if value is None:
        return None
    if not 0.0 <= value <= 1.0:
        raise ValueError(f"{name} must be between 0 and 1 inclusive")
    return value


def filter_consensus_by_sequence_mask(
    consensus: pr.PyRanges,
    fasta_path: Path,
    *,
    max_n_fraction: Optional[float] = None,
    max_lowercase_fraction: Optional[float] = None,
) -> tuple[pr.PyRanges, Dict[str, object], pd.DataFrame]:
    """Filter consensus peaks by sequence masking metrics from a reference FASTA."""

    if max_n_fraction is None and max_lowercase_fraction is None:
        raise ValueError("At least one sequence-mask threshold must be provided")
    if pysam is None:
        raise ImportError(
            "pysam is required for genome-FASTA masking filters; install it with 'pip install pysam'"
        )
    if not fasta_path.exists():
        raise FileNotFoundError(f"Genome FASTA not found: {fasta_path}")

    df = consensus.df.copy()
    if df.empty:
        summary = {
            "applied": True,
            "genome_fasta": str(fasta_path),
            "max_n_fraction": max_n_fraction,
            "max_lowercase_fraction": max_lowercase_fraction,
            "total_peaks": 0,
            "kept_peaks": 0,
            "removed_peaks": 0,
            "missing_sequence_peaks": 0,
            "mean_n_fraction_kept": None,
            "mean_lowercase_fraction_kept": None,
        }
        empty_report = pd.DataFrame(
            columns=[
                "Chromosome",
                "Start",
                "End",
                "Name",
                "Support",
                "length",
                "n_fraction",
                "lowercase_fraction",
                "missing_sequence",
                "keep",
                "filter_reason",
            ]
        )
        return consensus, summary, empty_report

    records: List[Dict[str, object]] = []
    try:
        fasta_reader = pysam.FastaFile(str(fasta_path))
    except (OSError, ValueError) as exc:
        raise RuntimeError(
            f"Unable to open/index genome FASTA {fasta_path}; ensure it is bgzip/faidx compatible and run 'samtools faidx' if needed"
        ) from exc

    with fasta_reader as fasta:
        references = set(fasta.references)
        for row in df.itertuples(index=False):
            chrom = str(row.Chromosome)
            start = int(row.Start)
            end = int(row.End)
            name = str(row.Name)
            support = getattr(row, "Support", pd.NA)

            seq = ""
            missing_sequence = chrom not in references
            if not missing_sequence:
                try:
                    seq = fasta.fetch(chrom, start, end)
                except (KeyError, ValueError, OSError):
                    missing_sequence = True

            seq_length = len(seq)
            n_fraction = math.nan
            lowercase_fraction = math.nan
            if not missing_sequence and seq_length > 0:
                n_fraction = sum(base.upper() == "N" for base in seq) / seq_length
                lowercase_fraction = sum(base.islower() for base in seq) / seq_length
            elif seq_length == 0:
                missing_sequence = True

            keep = True
            reasons: List[str] = []
            if missing_sequence:
                keep = False
                reasons.append("missing_sequence")
            if max_n_fraction is not None and not math.isnan(n_fraction) and n_fraction > max_n_fraction:
                keep = False
                reasons.append("n_fraction")
            if (
                max_lowercase_fraction is not None
                and not math.isnan(lowercase_fraction)
                and lowercase_fraction > max_lowercase_fraction
            ):
                keep = False
                reasons.append("lowercase_fraction")

            records.append(
                {
                    "Chromosome": chrom,
                    "Start": start,
                    "End": end,
                    "Name": name,
                    "Support": support,
                    "length": max(0, end - start),
                    "n_fraction": n_fraction,
                    "lowercase_fraction": lowercase_fraction,
                    "missing_sequence": missing_sequence,
                    "keep": keep,
                    "filter_reason": ",".join(reasons),
                }
            )

    metrics_df = pd.DataFrame(records)
    missing_sequence_count = int(metrics_df["missing_sequence"].sum())
    if missing_sequence_count == len(metrics_df):
        raise ValueError(
            "No consensus intervals could be fetched from the genome FASTA; check the genome build and chromosome naming."
        )

    kept_metrics = metrics_df[metrics_df["keep"]].copy()
    kept_df = kept_metrics[["Chromosome", "Start", "End", "Name", "Support"]].copy()
    filtered_consensus = pr.PyRanges(kept_df)

    mean_n_fraction = kept_metrics["n_fraction"].dropna().mean()
    mean_lowercase_fraction = kept_metrics["lowercase_fraction"].dropna().mean()
    summary = {
        "applied": True,
        "genome_fasta": str(fasta_path),
        "max_n_fraction": max_n_fraction,
        "max_lowercase_fraction": max_lowercase_fraction,
        "total_peaks": int(len(metrics_df)),
        "kept_peaks": int(len(kept_metrics)),
        "removed_peaks": int(len(metrics_df) - len(kept_metrics)),
        "missing_sequence_peaks": missing_sequence_count,
        "mean_n_fraction_kept": None if pd.isna(mean_n_fraction) else float(mean_n_fraction),
        "mean_lowercase_fraction_kept": None
        if pd.isna(mean_lowercase_fraction)
        else float(mean_lowercase_fraction),
    }
    return filtered_consensus, summary, metrics_df


# ---------------------------------------------------------------------------
# Counting with deepTools
# ---------------------------------------------------------------------------


def canonicalize_count_rows(counts: pd.DataFrame) -> pd.DataFrame:
    """Return interval-count rows in a deterministic genomic order.

    ``multiBamSummary`` may emit the same BED intervals in different row orders
    when worker scheduling changes.  Sorting before statistical analysis keeps
    PyDESeq2 input order and tie-broken exploratory ranks reproducible across
    thread counts and repeated runs.
    """

    required = ["Chromosome", "start", "end"]
    missing = [column for column in required if column not in counts.columns]
    if missing:
        raise ValueError(
            "Counts matrix is missing canonical coordinate columns: "
            + ", ".join(missing)
        )
    return counts.sort_values(required, kind="mergesort").reset_index(drop=True)


def run_multibamsummary(
    consensus_bed: Path,
    samples: List[SampleEntry],
    output_dir: Path,
    *,
    count_config: CountFilterConfig,
    threads: int = 1,
    command_log: Optional[List[List[str]]] = None,
) -> Path:
    """Count interval overlaps using the configured read or fragment unit."""

    count_config.validate()
    include_flags = {
        count_config.effective_include_flags(paired_end=bool(sample.is_paired))
        for sample in samples
    }
    if len(include_flags) != 1:
        raise ValueError(
            "All BAMs in one analysis must resolve to the same count-unit include flags"
        )
    include_flag = include_flags.pop()

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
        "--labels",
        *(sample.sample for sample in samples),
        "--outFileName",
        str(out_npz),
        "--outRawCounts",
        str(out_tsv),
        "--numberOfProcessors",
        str(threads),
    ])
    cmd.extend(["--minMappingQuality", str(count_config.min_mapq)])
    cmd.extend(["--samFlagExclude", str(count_config.effective_exclude_flags)])
    if include_flag:
        cmd.extend(["--samFlagInclude", str(include_flag)])
    if count_config.count_unit == "fragment":
        cmd.append("--extendReads")
    if command_log is not None:
        command_log.append(list(cmd))
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


def pydeseq2_differential(
    counts: pd.DataFrame,
    conditions: pd.Series,
    *,
    n_cpus: int = 1,
) -> pd.DataFrame:
    if DeseqDataSet is None or DeseqStats is None:
        raise ImportError("pydeseq2 is required for the DESeq2 workflow but is not installed")
    if n_cpus < 1:
        raise ValueError("PyDESeq2 n_cpus must be at least 1")

    logging.info("Running PyDESeq2 differential analysis")
    samples = counts.columns.tolist()
    cond = conditions.loc[samples]
    if cond.nunique() != 2:
        raise ValueError("PyDESeq2 differential analysis requires exactly two conditions")

    condition_order = list(dict.fromkeys(cond.tolist()))
    reference = condition_order[0]
    contrast = condition_order[1]

    metadata = pd.DataFrame({"condition": cond.astype("category")}, index=samples)
    dds = DeseqDataSet(
        counts=counts.T.astype(int),
        metadata=metadata,
        design="~condition",
        n_cpus=n_cpus,
    )
    dds.deseq2()

    stats = DeseqStats(
        dds,
        contrast=("condition", contrast, reference),
        n_cpus=n_cpus,
    )
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
    result["analysis_mode"] = "replicate_supported_inference"
    result["interpretation"] = "formal_inference_with_biological_replicates"
    return result


def _rank_subset(
    values: pd.Series,
    eligible: pd.Series,
    *,
    ascending: bool,
) -> pd.Series:
    ranks = pd.Series(pd.NA, index=values.index, dtype="Int64")
    mask = eligible.fillna(False) & values.notna()
    if mask.any():
        # Coordinate/peak identifiers provide a stable tie-break independent
        # of the row order returned by the counting backend.
        ordered = values.loc[mask].sort_index(kind="mergesort")
        ranked = ordered.rank(method="first", ascending=ascending).astype("Int64")
        ranks.loc[ranked.index] = ranked
    return ranks


def single_pair_statistics(
    counts: pd.DataFrame,
    conditions: pd.Series,
    library_sizes: pd.Series,
    *,
    minimum_mean_cpm: float = 1.0,
    pseudocount: float = 0.5,
) -> pd.DataFrame:
    """Compute exploratory effect-size and MARS-derived single-pair rankings.

    Sampling-model p/q values are retained as explicitly named diagnostics.
    They do not estimate biological variability and must not be interpreted as
    formal significance or FDR-controlled discoveries.
    """

    if minimum_mean_cpm < 0:
        raise ValueError("minimum_mean_cpm must be non-negative")
    if pseudocount <= 0:
        raise ValueError("pseudocount must be positive")

    logging.info("Running exploratory single-pair ranking")
    samples = counts.columns.tolist()
    cond_series = conditions.loc[samples]
    unique_conditions = list(dict.fromkeys(cond_series.tolist()))
    if len(unique_conditions) != 2:
        raise ValueError("Single-pair ranking requires exactly two conditions")

    condition_a, condition_b = unique_conditions
    condition_a_cols = cond_series[cond_series == condition_a].index
    condition_b_cols = cond_series[cond_series == condition_b].index
    if len(condition_a_cols) != 1 or len(condition_b_cols) != 1:
        raise ValueError(
            "Exploratory single-pair mode requires exactly one sample per condition"
        )

    count_a = counts.loc[:, condition_a_cols].iloc[:, 0].astype(float)
    count_b = counts.loc[:, condition_b_cols].iloc[:, 0].astype(float)
    if not isinstance(library_sizes, pd.Series):
        library_sizes = pd.Series(library_sizes)
    library_sizes = library_sizes.reindex(samples)
    if library_sizes.isna().any():
        missing = library_sizes[library_sizes.isna()].index.tolist()
        raise ValueError(
            "Library size information missing for samples: " + ", ".join(missing)
        )

    library_size_a = float(library_sizes.loc[list(condition_a_cols)].iloc[0])
    library_size_b = float(library_sizes.loc[list(condition_b_cols)].iloc[0])
    if library_size_a <= 0 or library_size_b <= 0:
        raise ValueError("Library sizes must be positive for both single-pair samples")

    c_a = count_a.to_numpy(dtype=float)
    c_b = count_b.to_numpy(dtype=float)
    if np.any(c_a < 0) or np.any(c_b < 0):
        raise ValueError("Counts must be non-negative")

    cpm_a = c_a / library_size_a * 1_000_000.0
    cpm_b = c_b / library_size_b * 1_000_000.0
    mean_cpm = 0.5 * (cpm_a + cpm_b)
    library_offset = math.log2(library_size_b / library_size_a)
    normalized_log2fc = (
        np.log2((c_b + pseudocount) / (c_a + pseudocount)) - library_offset
    )

    both_zero = (c_a == 0) & (c_b == 0)
    one_zero = (c_a == 0) ^ (c_b == 0)
    positive_counts = (c_a > 0) & (c_b > 0)

    with np.errstate(divide="ignore"):
        log2_b = np.log2(c_b)
        log2_a = np.log2(c_a)
    valid_mask = positive_counts & np.isfinite(log2_b) & np.isfinite(log2_a)

    M = np.full_like(c_a, np.nan, dtype=float)
    A = np.full_like(c_a, np.nan, dtype=float)
    M[valid_mask] = log2_b[valid_mask] - log2_a[valid_mask]
    A[valid_mask] = 0.5 * (log2_b[valid_mask] + log2_a[valid_mask])

    sqrt_total = math.sqrt(library_size_b * library_size_a)
    p = np.full_like(M, np.nan, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        p[valid_mask] = np.exp2(A[valid_mask]) / sqrt_total
    epsilon = 1e-12
    if np.any(valid_mask):
        p[valid_mask] = np.clip(p[valid_mask], epsilon, 1 - epsilon)

    log_factor = math.log(2.0)
    denom = (library_size_b + library_size_a) * p
    with np.errstate(divide="ignore", invalid="ignore"):
        variance = 4.0 * (1.0 - p) / (denom * (log_factor ** 2))
    sd = np.sqrt(variance)
    expected_null_m = np.full_like(M, library_offset, dtype=float)

    z_scores = np.full_like(M, np.nan, dtype=float)
    valid_z = valid_mask & np.isfinite(sd) & (sd > 0)
    z_scores[valid_z] = (M[valid_z] - expected_null_m[valid_z]) / sd[valid_z]
    pvals = np.ones_like(M, dtype=float)
    finite_z = np.isfinite(z_scores)
    pvals[finite_z] = 2.0 * stats.norm.sf(np.abs(z_scores[finite_z]))

    res_df = pd.DataFrame(
        {
            "Peak": counts.index,
            "condition_a": condition_a,
            "condition_b": condition_b,
            "count_condition_a": c_a,
            "count_condition_b": c_b,
            "library_size_a": library_size_a,
            "library_size_b": library_size_b,
            "cpm_condition_a": cpm_a,
            "cpm_condition_b": cpm_b,
            "mean_cpm": mean_cpm,
            "normalized_log2fc": normalized_log2fc,
            "A": A,
            "M": M,
            "expected_null_m": expected_null_m,
            "mars_variance": variance,
            "mars_sd": sd,
            "mars_score": z_scores,
            "sampling_pvalue": pvals,
        }
    ).set_index("Peak")
    res_df["sampling_qvalue"] = benjamini_hochberg(
        res_df["sampling_pvalue"]
    ).fillna(1.0)
    res_df["zero_count_status"] = np.select(
        [both_zero, one_zero],
        ["both_zero", "one_zero"],
        default="both_positive",
    )
    res_df["rank_eligible"] = (~both_zero) & (res_df["mean_cpm"] >= minimum_mean_cpm)
    res_df["rank_up"] = _rank_subset(
        res_df["normalized_log2fc"],
        res_df["rank_eligible"] & (res_df["normalized_log2fc"] > 0),
        ascending=False,
    )
    res_df["rank_down"] = _rank_subset(
        res_df["normalized_log2fc"],
        res_df["rank_eligible"] & (res_df["normalized_log2fc"] < 0),
        ascending=True,
    )
    res_df["rank_absolute_lfc"] = _rank_subset(
        res_df["normalized_log2fc"].abs(),
        res_df["rank_eligible"],
        ascending=False,
    )
    res_df["rank_absolute_mars"] = _rank_subset(
        res_df["mars_score"].abs(),
        res_df["rank_eligible"] & res_df["mars_score"].notna(),
        ascending=False,
    )
    res_df["method"] = "single_pair_effect_ranking_with_mars_diagnostic"
    res_df["analysis_mode"] = "single_pair_exploratory"
    res_df["interpretation"] = "ranking_only_no_biological_variance_estimation"
    return res_df


def mars_differential(
    counts: pd.DataFrame, conditions: pd.Series, library_sizes: pd.Series
) -> pd.DataFrame:
    """Deprecated compatibility wrapper for the exploratory single-pair API."""

    warnings.warn(
        "mars_differential() legacy pvalue/padj/log2FC aliases are deprecated; "
        "use single_pair_statistics() and the explicit exploratory column names.",
        DeprecationWarning,
        stacklevel=2,
    )
    result = single_pair_statistics(
        counts,
        conditions,
        library_sizes,
        minimum_mean_cpm=0.0,
    )
    result["log2FC"] = result["normalized_log2fc"]
    result["log2FC_shrunk"] = result["normalized_log2fc"]
    result["pvalue"] = result["sampling_pvalue"]
    result["padj"] = result["sampling_qvalue"]
    return result


# ---------------------------------------------------------------------------
# Differential orchestration helpers
# ---------------------------------------------------------------------------


def call_differential_analysis(
    counts: pd.DataFrame,
    conditions: pd.Series,
    library_sizes: pd.Series,
    *,
    minimum_mean_cpm: float = 1.0,
    pseudocount: float = 0.5,
    pydeseq2_n_cpus: int = 1,
) -> pd.DataFrame:
    """Select and run the appropriate differential analysis workflow."""

    if counts.empty:
        raise ValueError("Counts matrix is empty; cannot perform differential analysis")
    if conditions.nunique() != 2:
        raise ValueError("Differential analysis requires exactly two experimental conditions")

    group_sizes = conditions.value_counts()
    replicates_per_condition = group_sizes.min()
    if replicates_per_condition >= 2:
        logging.info("Detected replicates per condition; using DESeq2 analysis via PyDESeq2")
        return pydeseq2_differential(counts, conditions, n_cpus=pydeseq2_n_cpus)

    if not (group_sizes == 1).all():
        raise ValueError(
            "PeakForge supports either at least two replicates in both groups or an exact 1-vs-1 "
            "exploratory comparison; mixed replicated/unreplicated group sizes are unsupported"
        )

    logging.info("Detected exact 1-vs-1 design; using exploratory candidate ranking")
    return single_pair_statistics(
        counts,
        conditions,
        library_sizes,
        minimum_mean_cpm=minimum_mean_cpm,
        pseudocount=pseudocount,
    )


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
    single_pair = "normalized_log2fc" in results.columns
    if single_pair:
        base_mean = results["mean_cpm"]
        effect = results["normalized_log2fc"]
    else:
        base_mean = results.get("baseMean")
        if base_mean is None:
            base_mean = counts.mean(axis=1)
        effect = results["log2FC"]
    fig, ax = plt.subplots(figsize=(6, 6))
    A = np.log2(base_mean + 1e-6)
    sns.scatterplot(x=A, y=effect, ax=ax, s=10, color="steelblue")
    ax.axhline(0, color="black", linestyle="--", linewidth=0.8)
    ax.set_xlabel("Mean log2 CPM" if single_pair else "Average log2 count")
    ax.set_ylabel("Normalized log2 fold change" if single_pair else "log2 fold change")
    ax.set_title("Exploratory effect–abundance plot" if single_pair else "MA plot")
    save_plot(fig, output)


def plot_sample_correlation(counts: pd.DataFrame, output: Path) -> None:
    corr = counts.apply(lambda x: np.log1p(x)).corr()
    fig, ax = plt.subplots(figsize=(6, 6))
    sns.heatmap(corr, cmap="viridis", ax=ax, annot=True, fmt=".2f")
    ax.set_title("Sample correlation")
    save_plot(fig, output)


def plot_top_heatmap(counts: pd.DataFrame, results: pd.DataFrame, output: Path, top_n: int = 50) -> None:
    if "rank_absolute_lfc" in results.columns:
        top = results.sort_values("rank_absolute_lfc", na_position="last").head(top_n).index
        title = f"Top {top_n} exploratory candidates by absolute effect"
    else:
        top = results.sort_values("padj").head(top_n).index
        title = f"Top {top_n} differential peaks"
    data = counts.loc[top]
    log_data = np.log2(data + 1)
    norm = log_data.sub(log_data.mean(axis=1), axis=0)
    fig = plt.figure(figsize=(8, max(4, len(top) * 0.2)))
    sns.heatmap(norm, cmap="RdBu_r", center=0)
    plt.title(title)
    plt.ylabel("Peaks")
    plt.xlabel("Samples")
    save_plot(fig, output)


def plot_single_pair_ranking(results: pd.DataFrame, output: Path) -> None:
    """Plot an exploratory effect-size ranking without significance thresholds."""

    required = {"normalized_log2fc", "rank_absolute_lfc", "rank_eligible"}
    missing = required - set(results.columns)
    if missing:
        raise ValueError(
            "Single-pair ranking plot is missing columns: " + ", ".join(sorted(missing))
        )
    ranked = results.loc[results["rank_eligible"].fillna(False)].copy()
    ranked = ranked.dropna(subset=["rank_absolute_lfc", "normalized_log2fc"])
    ranked = ranked.sort_values("rank_absolute_lfc")
    fig, ax = plt.subplots(figsize=(7, 5))
    colors = np.where(ranked["normalized_log2fc"] >= 0, "#B2182B", "#2166AC")
    ax.scatter(
        ranked["rank_absolute_lfc"].astype(float),
        ranked["normalized_log2fc"],
        c=colors,
        s=12,
        alpha=0.65,
        linewidths=0,
    )
    ax.axhline(0, color="black", linestyle="--", linewidth=0.8)
    ax.set_xlabel("Absolute-effect rank")
    ax.set_ylabel("Normalized log2 fold change")
    ax.set_title("Exploratory single-pair candidate ranking")
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
    sns.scatterplot(
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
    if "normalized_log2fc" in results.columns:
        outputs = {
            "sample_correlation": output_dir / "sample_correlation.png",
            "effect_abundance": output_dir / "effect_abundance.png",
            "candidate_ranking": output_dir / "candidate_ranking.png",
            "top_heatmap": output_dir / "top_candidate_heatmap.png",
        }
        plot_sample_correlation(counts, outputs["sample_correlation"])
        plot_ma(results, counts, outputs["effect_abundance"])
        plot_single_pair_ranking(results, outputs["candidate_ranking"])
        plot_top_heatmap(counts, results, outputs["top_heatmap"])
        return outputs

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

    conditions = pd.Series({s.sample: s.condition for s in samples})

    consensus_arg = getattr(args, "consensus_peaks", None)
    consensus_path = Path(consensus_arg) if consensus_arg else None
    blacklist_bed_arg = getattr(args, "blacklist_bed", None)
    blacklist_bed_path = Path(blacklist_bed_arg) if blacklist_bed_arg else None
    genome_fasta_arg = getattr(args, "genome_fasta", None)
    genome_fasta_path = Path(genome_fasta_arg) if genome_fasta_arg else None
    motif_rank_requested = bool(getattr(args, "motif_rank", False))
    motif_genome_fasta_arg = getattr(args, "motif_genome_fasta", None)
    motif_genome_fasta_path = Path(motif_genome_fasta_arg) if motif_genome_fasta_arg else genome_fasta_path
    motif_max_fraction = validate_fraction_threshold(
        "--motif-max-fraction",
        getattr(args, "motif_max_fraction", 0.8),
    )
    max_n_fraction = validate_fraction_threshold(
        "--max-n-fraction",
        getattr(args, "max_n_fraction", None),
    )
    max_lowercase_fraction = validate_fraction_threshold(
        "--max-lowercase-fraction",
        getattr(args, "max_lowercase_fraction", None),
    )
    sequence_filter_requested = max_n_fraction is not None or max_lowercase_fraction is not None
    filters_requested = blacklist_bed_path is not None or sequence_filter_requested
    if sequence_filter_requested and genome_fasta_path is None:
        raise ValueError("--genome-fasta is required when using sequence masking filters")
    if motif_rank_requested and motif_genome_fasta_path is None:
        raise ValueError("--motif-genome-fasta or --genome-fasta is required when using --motif-rank")
    if motif_rank_requested and not motif_genome_fasta_path.exists():
        raise FileNotFoundError(f"Motif genome FASTA not found: {motif_genome_fasta_path}")

    required_cmds = ["multiBamSummary", "samtools"]
    needs_peak_calling = consensus_path is None and any(sample.peaks is None for sample in samples)
    if needs_peak_calling:
        required_cmds.append(get_macs_command())
    ensure_commands(required_cmds)

    samtools_path = shutil.which("samtools")
    if samtools_path is None:  # pragma: no cover - defensive (ensure_commands guards)
        raise RuntimeError("samtools not found on PATH after validation")

    for sample in samples:
        sample.is_paired = detect_paired_end_bam(sample.bam, samtools_path, args.threads)

    results_dir = ensure_directory(Path(args.output_dir))
    count_config = CountFilterConfig(
        count_unit=getattr(args, "count_unit", "read"),
        min_mapq=getattr(args, "min_mapq", 0),
        exclude_duplicates=bool(getattr(args, "exclude_duplicates", False)),
        proper_pairs_only=bool(getattr(args, "proper_pairs_only", False)),
        sam_flag_exclude=getattr(args, "sam_flag_exclude", 0),
    )
    library_size_metrics = compute_library_size_metrics(
        samples,
        samtools_path,
        count_config,
        args.threads,
    )
    library_sizes = library_size_metrics["selected_library_size"].astype(float)
    library_size_metrics_path = results_dir / "library_size_metrics.tsv"
    library_size_metrics.to_csv(library_size_metrics_path, sep="\t")
    executed_commands: List[List[str]] = []

    consensus: pr.PyRanges
    consensus_bed = results_dir / "consensus_peaks.bed"
    consensus_metadata: Dict[str, object] = {
        "source": "generated" if consensus_path is None else "provided",
        "input": str(consensus_path) if consensus_path else None,
        "path": str(consensus_bed),
        "peak_coordinate_definition": {
            "applied_to_sample_peaks": consensus_path is None,
            "mode": getattr(args, "peak_coordinate_mode", "edge-padding"),
            "edge_padding_bp_per_side": getattr(args, "peak_edge_padding", 250),
            "summit_fixed_total_width_bp": getattr(args, "summit_fixed_width", 500),
        },
    }

    if consensus_path is None:
        macs2_params = {
            "genome": args.macs2_genome,
            "qvalue": args.macs2_qvalue,
            "format": getattr(args, "macs_format", "AUTO"),
            "nomodel": bool(getattr(args, "macs_nomodel", False)),
            "shift": getattr(args, "macs_shift", None),
            "extsize": getattr(args, "macs_extsize", None),
            "extra": args.macs2_extra,
        }
        peak_ranges = load_all_peaks(
            samples,
            peak_edge_padding=getattr(args, "peak_edge_padding", 250),
            peak_coordinate_mode=getattr(args, "peak_coordinate_mode", "edge-padding"),
            summit_fixed_width=getattr(args, "summit_fixed_width", 500),
            default_peak_type=args.peak_type,
            macs2_params=macs2_params,
            peak_output_dir=Path(args.peak_dir),
            command_log=executed_commands,
        )

        consensus = build_consensus(peak_ranges, min_overlap=args.min_overlap)
    else:
        logging.info("Using provided consensus peaks: %s", consensus_path)
        consensus = load_consensus_bed(consensus_path)

    prefilter_path: Optional[Path] = results_dir / "consensus_peaks.prefilter.bed" if filters_requested else None
    if prefilter_path is not None:
        if consensus_path is None:
            write_consensus_bed(consensus, prefilter_path)
        else:
            ensure_directory(prefilter_path.parent)
            if consensus_path.resolve() != prefilter_path.resolve():
                shutil.copyfile(consensus_path, prefilter_path)

    blacklist_filter_metadata: Dict[str, object] = {
        "applied": False,
        "blacklist_bed": str(blacklist_bed_path) if blacklist_bed_path else None,
        "report": None,
        "prefilter_path": str(prefilter_path) if prefilter_path and blacklist_bed_path else None,
    }
    if blacklist_bed_path is not None:
        blacklist_report_path = results_dir / "consensus_blacklist_filter.tsv"
        consensus, blacklist_filter_metadata, blacklist_report_df = filter_consensus_by_blacklist(
            consensus,
            blacklist_bed_path,
        )
        blacklist_report_df.to_csv(blacklist_report_path, sep="\t", index=False)
        blacklist_filter_metadata["report"] = str(blacklist_report_path)
        blacklist_filter_metadata["prefilter_path"] = str(prefilter_path) if prefilter_path else None
        logging.info(
            "Blacklist filter retained %d/%d consensus peaks (removed %d; overlapped %d)",
            blacklist_filter_metadata["kept_peaks"],
            blacklist_filter_metadata["total_peaks"],
            blacklist_filter_metadata["removed_peaks"],
            blacklist_filter_metadata["overlapped_peaks"],
        )

    sequence_filter_metadata: Dict[str, object] = {
        "applied": False,
        "genome_fasta": str(genome_fasta_path) if genome_fasta_path else None,
        "max_n_fraction": max_n_fraction,
        "max_lowercase_fraction": max_lowercase_fraction,
        "report": None,
        "prefilter_path": str(prefilter_path) if prefilter_path and sequence_filter_requested else None,
    }
    if sequence_filter_requested:
        metrics_path = results_dir / "consensus_sequence_filter.tsv"
        consensus, sequence_filter_metadata, metrics_df = filter_consensus_by_sequence_mask(
            consensus,
            genome_fasta_path,
            max_n_fraction=max_n_fraction,
            max_lowercase_fraction=max_lowercase_fraction,
        )
        metrics_df.to_csv(metrics_path, sep="\t", index=False)
        sequence_filter_metadata["report"] = str(metrics_path)
        sequence_filter_metadata["prefilter_path"] = str(prefilter_path)
        logging.info(
            "Sequence mask filter retained %d/%d consensus peaks (removed %d; missing sequence for %d)",
            sequence_filter_metadata["kept_peaks"],
            sequence_filter_metadata["total_peaks"],
            sequence_filter_metadata["removed_peaks"],
            sequence_filter_metadata["missing_sequence_peaks"],
        )

    consensus_metadata["blacklist_filter"] = blacklist_filter_metadata
    consensus_metadata["sequence_filter"] = sequence_filter_metadata

    if len(consensus) == 0:
        raise ValueError("Consensus peak set is empty; cannot proceed with counting")

    if filters_requested or consensus_path is None:
        write_consensus_bed(consensus, consensus_bed)
    else:
        ensure_directory(consensus_bed.parent)
        if consensus_path.resolve() != consensus_bed.resolve():
            shutil.copyfile(consensus_path, consensus_bed)

    counts_tsv = run_multibamsummary(
        consensus_bed,
        samples,
        results_dir / "counts",
        count_config=count_config,
        threads=args.threads,
        command_log=executed_commands,
    )
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
    raw_counts = canonicalize_count_rows(raw_counts)
    # Replace the backend-specific raw ordering/header with the canonical count
    # matrix used downstream so the persisted artifact is itself reproducible.
    raw_counts.to_csv(counts_tsv, sep="\t", index=False)

    consensus_df = consensus.df.copy().rename(columns={"Start": "start", "End": "end"})
    merged = raw_counts.merge(consensus_df[["Chromosome", "start", "end", "Name"]],
                              on=["Chromosome", "start", "end"], how="left")
    merged["Peak"] = merged["Name"].fillna(
        merged["Chromosome"].astype(str)
        + ":"
        + merged["start"].astype(int).astype(str)
        + "-"
        + merged["end"].astype(int).astype(str)
    )
    merged.sort_values(
        ["Chromosome", "start", "end", "Peak"],
        kind="mergesort",
        inplace=True,
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

    diff_res = call_differential_analysis(
        counts_df,
        conditions,
        library_sizes,
        minimum_mean_cpm=getattr(args, "single_pair_min_mean_cpm", 1.0),
        pseudocount=getattr(args, "single_pair_pseudocount", 0.5),
        pydeseq2_n_cpus=getattr(args, "pydeseq2_cpus", 1),
    )

    diff_path = results_dir / "differential_results.tsv"
    diff_res.to_csv(diff_path, sep="\t")

    motif_ranking_paths = None
    if motif_rank_requested:
        condition_order = list(dict.fromkeys(sample.condition for sample in samples))
        if len(condition_order) != 2:
            raise ValueError("Motif ranking requires exactly two conditions")
        negative_condition, positive_condition = condition_order
        motif_file_arg = getattr(args, "motif_file", None)
        motif_file_path = Path(motif_file_arg) if motif_file_arg else None
        motif_ranking_paths = run_pairwise_motif_ranking(
            consensus_bed=consensus_bed,
            diff_res=diff_res,
            output_dir=results_dir / "motif_ranking",
            genome_fasta=motif_genome_fasta_path,
            motif_file=motif_file_path,
            positive_condition=positive_condition,
            negative_condition=negative_condition,
            score_metric=args.motif_score_metric,
            gsea_weight=args.motif_gsea_weight,
            min_peaks=args.motif_min_peaks,
            max_fraction=0.8 if motif_max_fraction is None else motif_max_fraction,
            top_peaks=args.motif_top_peaks,
            threads=args.threads,
        )

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
            if "rank_absolute_lfc" in diff_res.columns:
                top_peaks = (
                    diff_res.loc[diff_res["rank_eligible"].fillna(False)]
                    .sort_values("rank_absolute_lfc", na_position="last")
                    .head(args.enrichr_top)
                    .index
                )
                logging.warning(
                    "Single-pair Enrichr input is an exploratory effect-size ranking, not a set of significant peaks"
                )
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
                "control_bam": str(sample.control_bam) if sample.control_bam is not None else None,
                "peaks": str(sample.peaks) if sample.peaks is not None else None,
                "peak_type": sample.peak_type,
                "library_size": float(library_sizes.loc[sample.sample]),
                "library_size_metrics": {
                    key: (
                        value.item()
                        if isinstance(value, np.generic)
                        else value
                    )
                    for key, value in library_size_metrics.loc[sample.sample].to_dict().items()
                },
            }
        )

    args_dict = {key: value for key, value in vars(args).items() if key != "samples"}
    metadata = {
        "timestamp": datetime.utcnow().isoformat(),
        "peakforge_version": __version__,
        "git": git_state(),
        "software_versions": software_versions(),
        "invocation": [sys.executable, *sys.argv],
        "args": args_dict,
        "metadata_sheet": str(metadata_path) if metadata_path else None,
        "samples": sample_metadata,
        "counts_matrix": str(counts_tsv),
        "differential_results": str(diff_path),
        "plots": {key: str(path) for key, path in plot_paths.items()},
        "motif_ranking": {key: str(path) for key, path in motif_ranking_paths.items()}
        if motif_ranking_paths
        else None,
        "annotation": str(annotation_path) if annotation_path else None,
        "enrichr": str(enrichr_path) if enrichr_path else None,
        "library_sizes": {name: float(value) for name, value in library_sizes.to_dict().items()},
        "library_size_metrics": str(library_size_metrics_path),
        "count_filter": count_config.as_metadata(),
        "executed_commands": executed_commands,
        "analysis_mode": str(diff_res["analysis_mode"].iloc[0])
        if "analysis_mode" in diff_res.columns and not diff_res.empty
        else None,
        "consensus": consensus_metadata,
    }
    save_metadata(metadata, results_dir / "metadata.json")


# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------

def _resolve_optional_paths(
    values: Optional[Sequence[str]],
    *,
    expected_len: int,
    label: str,
    allow_broadcast: bool = False,
) -> List[Optional[Path]]:
    if values is None or len(values) == 0:
        return [None] * expected_len

    if allow_broadcast and len(values) == 1 and expected_len > 1:
        values = list(values) * expected_len

    if len(values) != expected_len:
        if allow_broadcast:
            raise ValueError(
                f"{label} must contain either one path or exactly {expected_len} paths"
            )
        raise ValueError(f"{label} must contain exactly {expected_len} paths")

    resolved: List[Optional[Path]] = []
    for value in values:
        if value in {"", "-", "None", "none", "NA", "na"}:
            resolved.append(None)
        else:
            resolved.append(Path(value))
    return resolved


def _build_condition_samples(
    condition: str,
    bam_files: Sequence[str],
    peak_files: Optional[Sequence[str]],
    control_bams: Optional[Sequence[str]],
    used_names: Set[str],
) -> List[SampleEntry]:
    if not bam_files:
        raise ValueError(f"No BAM files supplied for condition {condition}")

    resolved_peaks = _resolve_optional_paths(
        peak_files,
        expected_len=len(bam_files),
        label=f"Peak files for condition {condition}",
        allow_broadcast=False,
    )
    resolved_controls = _resolve_optional_paths(
        control_bams,
        expected_len=len(bam_files),
        label=f"Control/Input BAM files for condition {condition}",
        allow_broadcast=True,
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

        entries.append(
            SampleEntry(
                sample=candidate,
                condition=condition,
                bam=bam_path,
                control_bam=resolved_controls[idx],
                peaks=resolved_peaks[idx],
                peak_type="auto",
            )
        )

    return entries


def build_runmode_samples(args: argparse.Namespace) -> List[SampleEntry]:
    a_peaks = args.a_peaks or []
    b_peaks = args.b_peaks or []
    a_controls = args.a_control_bams or []
    b_controls = args.b_control_bams or []

    used_names: Set[str] = set()
    samples = _build_condition_samples(
        args.condition_a,
        args.a_bams,
        a_peaks,
        a_controls,
        used_names,
    )
    samples.extend(
        _build_condition_samples(
            args.condition_b,
            args.b_bams,
            b_peaks,
            b_controls,
            used_names,
        )
    )
    return samples


def write_sample_sheet(samples: Sequence[SampleEntry], output_path: Path) -> Path:
    """Export sample metadata as TSV/CSV based on output extension."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    sep = "\t" if output_path.suffix.lower() != ".csv" else ","
    frame = pd.DataFrame(
        [
            {
                "sample": sample.sample,
                "condition": sample.condition,
                "bam": str(sample.bam),
                "control_bam": str(sample.control_bam) if sample.control_bam else "",
                "peaks": str(sample.peaks) if sample.peaks else "",
                "peak_type": sample.peak_type,
            }
            for sample in samples
        ]
    )
    frame.to_csv(output_path, sep=sep, index=False)
    return output_path


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
        "--peak-edge-padding",
        "--peak-extension",
        dest="peak_edge_padding",
        type=int,
        default=250,
        help=(
            "Base pairs added to each original narrow-peak interval edge before consensus "
            "construction; does not use the midpoint or summit"
        ),
    )
    parser.add_argument(
        "--peak-coordinate-mode",
        choices=["edge-padding", "summit-fixed"],
        default="edge-padding",
        help=(
            "Coordinate definition used before consensus construction: edge-padding keeps "
            "native narrowPeak boundaries plus --peak-edge-padding; summit-fixed uses "
            "narrowPeak column 10 as the centre of an exact-width interval"
        ),
    )
    parser.add_argument(
        "--summit-fixed-width",
        type=int,
        default=500,
        help=(
            "Exact total interval width in bp when --peak-coordinate-mode summit-fixed; "
            "ignored in edge-padding mode"
        ),
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
        "--macs-format",
        choices=["AUTO", "BAM", "BAMPE"],
        default="AUTO",
        help="MACS input format; AUTO resolves BAMPE for paired-end BAMs and BAM otherwise",
    )
    parser.add_argument(
        "--macs-nomodel",
        action="store_true",
        help="Pass --nomodel to MACS",
    )
    parser.add_argument(
        "--macs-shift",
        type=int,
        default=None,
        help="Pass an assay-specific --shift value to MACS",
    )
    parser.add_argument(
        "--macs-extsize",
        type=int,
        default=None,
        help="Pass an assay-specific --extsize value to MACS",
    )
    parser.add_argument(
        "--macs2-extra",
        "--macs-extra",
        dest="macs2_extra",
        nargs=argparse.REMAINDER,
        default=[],
        help="Additional arguments for MACS2",
    )
    parser.add_argument(
        "--count-unit",
        choices=["read", "fragment"],
        default="read",
        help=(
            "Count aligned reads or paired-end DNA fragments; fragment mode counts one "
            "proper pair once and extends the first mate across the fragment"
        ),
    )
    parser.add_argument(
        "--min-mapq",
        type=int,
        default=0,
        help="Minimum mapping quality used for both interval counts and library size",
    )
    parser.add_argument(
        "--exclude-duplicates",
        action="store_true",
        help="Exclude alignments marked as PCR/optical duplicates from counts and library size",
    )
    parser.add_argument(
        "--proper-pairs-only",
        action="store_true",
        help="Restrict read mode to proper pairs; fragment mode always requires proper pairs",
    )
    parser.add_argument(
        "--sam-flag-exclude",
        type=int,
        default=0,
        help=(
            "Additional SAM flag bitmask to exclude; unmapped, secondary, QC-fail, and "
            "supplementary alignments are always excluded"
        ),
    )
    parser.add_argument(
        "--single-pair-min-mean-cpm",
        type=float,
        default=1.0,
        help="Minimum mean CPM for eligibility in exploratory single-pair ranking",
    )
    parser.add_argument(
        "--single-pair-pseudocount",
        type=float,
        default=0.5,
        help="Count-scale pseudocount for exploratory normalized single-pair log2 fold change",
    )
    parser.add_argument(
        "--blacklist-bed",
        help="BED file of blacklist regions to exclude from the consensus set before counting",
    )
    parser.add_argument(
        "--genome-fasta",
        help="Reference FASTA used for optional consensus masking filters",
    )
    parser.add_argument(
        "--motif-rank",
        action="store_true",
        help="Run pair-wise motif ranking on all consensus peaks using a preranked, GSEA-like motif score",
    )
    parser.add_argument(
        "--motif-genome-fasta",
        help="Reference FASTA used for motif scanning; defaults to --genome-fasta when supplied",
    )
    parser.add_argument(
        "--motif-file",
        help="Motif database file for pair-wise motif ranking (defaults to HOMER vertebrate known motifs)",
    )
    parser.add_argument(
        "--motif-score-metric",
        default="auto",
        choices=["auto", "signed_product", "signed_log10p", "signed_lfc", "signed_mars"],
        help=(
            "Peak ranking metric used before motif enrichment; auto uses signed normalized "
            "log2FC for single-pair output and signed effect-by-pvalue for replicated output"
        ),
    )
    parser.add_argument(
        "--motif-gsea-weight",
        type=float,
        default=1.0,
        help="Weight exponent for the GSEA-like motif running enrichment score",
    )
    parser.add_argument(
        "--motif-min-peaks",
        type=int,
        default=10,
        help="Minimum number of motif-bearing peaks required to score a motif",
    )
    parser.add_argument(
        "--motif-max-fraction",
        type=float,
        default=0.8,
        help="Skip motifs present in more than this fraction of consensus peaks",
    )
    parser.add_argument(
        "--motif-top-peaks",
        type=int,
        default=25,
        help="Number of top supporting peaks to report per motif",
    )
    parser.add_argument(
        "--max-n-fraction",
        type=float,
        default=None,
        help="Drop consensus peaks whose reference sequence has N fraction above this threshold",
    )
    parser.add_argument(
        "--max-lowercase-fraction",
        type=float,
        default=None,
        help="Drop consensus peaks whose reference sequence has soft-masked lowercase fraction above this threshold",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=16,
        help="Threads for multiBamSummary (--numberOfProcessors)",
    )
    parser.add_argument(
        "--pydeseq2-cpus",
        type=int,
        default=1,
        help="Worker processes used by PyDESeq2 in replicate-supported mode",
    )
    parser.add_argument("--gtf", help="Optional GTF file for annotation")
    parser.add_argument("--enrichr", action="store_true", help="Run Enrichr GO Biological Process analysis")
    parser.add_argument("--enrichr-top", type=int, default=200, help="Number of top peaks for enrichment")
    parser.add_argument("--log-level", default="INFO", help="Logging level")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="peakforge",
        description=(
            "Two-group narrow-peak workflow with replicate-supported inference and "
            "exploratory exact-1-vs-1 ranking"
        ),
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
    run_parser.add_argument(
        "--a-control-bams",
        "--a-controls",
        nargs="*",
        default=None,
        help="Optional control/input BAMs for condition A; supply one path to share across all A samples or one per BAM",
    )
    run_parser.add_argument("--condition-b", required=True, help="Contrast condition label")
    run_parser.add_argument("--b-bams", nargs="+", required=True, help="BAM files for condition B")
    run_parser.add_argument(
        "--b-peaks",
        nargs="*",
        default=None,
        help="Optional peak files aligned with --b-bams",
    )
    run_parser.add_argument(
        "--b-control-bams",
        "--b-controls",
        nargs="*",
        default=None,
        help="Optional control/input BAMs for condition B; supply one path to share across all B samples or one per BAM",
    )
    add_common_arguments(run_parser)

    peakshape_parser = subparsers.add_parser(
        "peakshape",
        help="Run peak shape profiling for two bigWig tracks",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    peak_shape.add_cli_arguments(peakshape_parser)

    sheet_parser = subparsers.add_parser(
        "makesheet",
        help="Generate a metadata TSV/CSV from runmode-style arguments",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    sheet_parser.add_argument("--condition-a", required=True, help="First condition label")
    sheet_parser.add_argument("--a-bams", nargs="+", required=True, help="BAM files for condition A")
    sheet_parser.add_argument(
        "--a-peaks",
        nargs="*",
        default=None,
        help="Optional peak files aligned with --a-bams",
    )
    sheet_parser.add_argument(
        "--a-control-bams",
        "--a-controls",
        nargs="*",
        default=None,
        help="Optional control/input BAMs for condition A; supply one path to share across all A samples or one per BAM",
    )
    sheet_parser.add_argument("--condition-b", required=True, help="Second condition label")
    sheet_parser.add_argument("--b-bams", nargs="+", required=True, help="BAM files for condition B")
    sheet_parser.add_argument(
        "--b-peaks",
        nargs="*",
        default=None,
        help="Optional peak files aligned with --b-bams",
    )
    sheet_parser.add_argument(
        "--b-control-bams",
        "--b-controls",
        nargs="*",
        default=None,
        help="Optional control/input BAMs for condition B; supply one path to share across all B samples or one per BAM",
    )
    sheet_parser.add_argument(
        "--output",
        required=True,
        help="Output metadata file path (.tsv recommended, .csv supported)",
    )
    sheet_parser.add_argument("--log-level", default="INFO", help="Logging level")

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
        elif args.command == "makesheet":
            samples = build_runmode_samples(args)
            output_path = write_sample_sheet(samples, Path(args.output))
            logging.info("Metadata sheet written to %s", output_path)
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
