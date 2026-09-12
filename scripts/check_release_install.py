#!/usr/bin/env python3
"""Smoke-test the installed wheel, not modules in a source checkout."""

from __future__ import annotations

import argparse
from importlib.metadata import version
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import chipdiff
import numpy as np
import pandas as pd
import pysam


def write_bam(path: Path, starts: list[int]) -> None:
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 1000}],
    }
    with pysam.AlignmentFile(path, "wb", header=header) as output:
        for index, start in enumerate(starts):
            for flag, offset, mate_offset, length in (
                (99, 0, 50, 100),
                (147, 50, 0, -100),
            ):
                read = pysam.AlignedSegment(output.header)
                read.query_name = f"fragment_{index}"
                read.query_sequence = "A" * 50
                read.flag = flag
                read.reference_id = read.next_reference_id = 0
                read.reference_start = start + offset
                read.next_reference_start = start + mate_offset
                read.template_length = length
                read.mapping_quality = 60
                read.cigar = ((0, 50),)
                read.query_qualities = pysam.qualitystring_to_array("I" * 50)
                output.write(read)
    pysam.index(str(path))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", default="0.2.3")
    args = parser.parse_args()
    prefix = Path(sys.prefix).resolve()
    module = Path(chipdiff.__file__).resolve()
    assert module.is_relative_to(prefix), (
        f"Imported checkout instead of installed wheel: {module}"
    )
    assert version("peakforge") == args.version
    # A tiny deterministic count matrix exercises the replicated dependency stack.
    rng = np.random.default_rng(20260912)
    means = np.repeat(rng.uniform(20, 200, size=(128, 1)), 6, axis=1)
    means[:16, 3:] *= 2
    names = ["A1", "A2", "A3", "B1", "B2", "B3"]
    counts = pd.DataFrame(rng.negative_binomial(10, 10 / (10 + means)), columns=names)
    replicated = chipdiff.pydeseq2_differential(
        counts, pd.Series(["A"] * 3 + ["B"] * 3, index=names), n_cpus=1
    )
    assert len(replicated) == 128
    assert {"pvalue", "padj"}.issubset(replicated)
    assert set(replicated.analysis_mode) == {"replicate_supported_inference"}
    assert np.isfinite(replicated.log2FC).all()
    environment = dict(
        os.environ, PATH=str(prefix / "bin") + os.pathsep + os.environ.get("PATH", "")
    )
    with tempfile.TemporaryDirectory(prefix="peakforge-installed-smoke-") as temporary:
        work = Path(temporary)
        write_bam(work / "a.bam", [100, 200])
        write_bam(work / "b.bam", [100, 300])
        (work / "consensus.bed").write_text(
            "chr1\t90\t180\tp1\nchr1\t190\t280\tp2\nchr1\t290\t380\tp3\n"
        )
        subprocess.run(
            [
                str(prefix / "bin/peakforge"),
                "runmode",
                "--condition-a",
                "A",
                "--a-bams",
                str(work / "a.bam"),
                "--condition-b",
                "B",
                "--b-bams",
                str(work / "b.bam"),
                "--consensus-peaks",
                str(work / "consensus.bed"),
                "--count-unit",
                "fragment",
                "--threads",
                "1",
                "--output-dir",
                str(work / "results"),
            ],
            check=True,
            cwd=work,
            env=environment,
        )
        results = pd.read_csv(work / "results/differential_results.tsv", sep="\t")
        assert len(results) == 3
        assert set(results.analysis_mode) == {"single_pair_exploratory"}
        assert "padj" not in results
        metadata = json.loads((work / "results/metadata.json").read_text())
        assert metadata["peakforge_version"] == args.version
    print(f"Installed-wheel smoke test passed: PeakForge {args.version} ({module})")


if __name__ == "__main__":
    main()
