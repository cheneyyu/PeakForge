from __future__ import annotations

import shutil
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pysam
import pandas as pd

from chipdiff import (
    CountFilterConfig,
    SampleEntry,
    canonicalize_count_rows,
    compute_library_size_metrics,
    run_multibamsummary,
)


def _segment(
    header: pysam.AlignmentHeader,
    name: str,
    flag: int,
    start: int,
    *,
    mate_start: int = -1,
    template_length: int = 0,
    mapq: int = 60,
) -> pysam.AlignedSegment:
    read = pysam.AlignedSegment(header)
    read.query_name = name
    read.query_sequence = "A" * 50
    read.flag = flag
    read.reference_id = 0
    read.reference_start = start
    read.mapping_quality = mapq
    read.cigar = ((0, 50),)
    read.next_reference_id = 0 if mate_start >= 0 else -1
    read.next_reference_start = mate_start
    read.template_length = template_length
    read.query_qualities = pysam.qualitystring_to_array("I" * 50)
    return read


def _write_synthetic_bam(path: Path) -> None:
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    )
    reads = [
        _segment(header, "good", 99, 100, mate_start=150, template_length=100),
        _segment(header, "good", 147, 150, mate_start=100, template_length=-100),
        _segment(header, "duplicate", 99 | 1024, 200, mate_start=250, template_length=100),
        _segment(header, "duplicate", 147 | 1024, 250, mate_start=200, template_length=-100),
        _segment(header, "lowmap", 99, 300, mate_start=350, template_length=100, mapq=10),
        _segment(header, "lowmap", 147, 350, mate_start=300, template_length=-100, mapq=10),
        _segment(header, "secondary", 256, 400),
        _segment(header, "supplementary", 2048, 450),
        _segment(header, "improper", 65, 500, mate_start=550, template_length=100),
        _segment(header, "improper", 129, 550, mate_start=500, template_length=-100),
    ]
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for read in reads:
            out.write(read)
    pysam.index(str(path))


@unittest.skipUnless(shutil.which("samtools"), "samtools is required")
class CountingTests(unittest.TestCase):
    def test_count_rows_are_canonicalized_independent_of_backend_order(self) -> None:
        frame = pd.DataFrame(
            {
                "Chromosome": ["chr2", "chr1", "chr1"],
                "start": [20, 30, 10],
                "end": [30, 40, 20],
                "sample": [3.0, 2.0, 1.0],
            }
        )
        expected = canonicalize_count_rows(frame)
        observed = canonicalize_count_rows(frame.iloc[::-1])
        pd.testing.assert_frame_equal(expected, observed)
        self.assertEqual(
            list(zip(expected["Chromosome"], expected["start"])),
            [("chr1", 10), ("chr1", 30), ("chr2", 20)],
        )

    def test_fragment_library_size_uses_one_filtered_proper_pair(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            bam = Path(tmpdir) / "tiny.bam"
            _write_synthetic_bam(bam)
            sample = SampleEntry("tiny", "A", bam, is_paired=True)
            config = CountFilterConfig(
                count_unit="fragment",
                min_mapq=30,
                exclude_duplicates=True,
            )
            metrics = compute_library_size_metrics(
                [sample], shutil.which("samtools"), config, threads=1
            )

        row = metrics.loc["tiny"]
        self.assertEqual(int(row["raw_mapped_alignments"]), 10)
        self.assertEqual(int(row["filtered_alignments"]), 2)
        self.assertEqual(int(row["filtered_countable_units"]), 1)
        self.assertEqual(int(row["selected_library_size"]), 1)

    def test_read_library_size_counts_both_mates_and_improper_reads(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            bam = Path(tmpdir) / "tiny.bam"
            _write_synthetic_bam(bam)
            sample = SampleEntry("tiny", "A", bam, is_paired=True)
            config = CountFilterConfig(
                count_unit="read",
                min_mapq=30,
                exclude_duplicates=True,
                proper_pairs_only=False,
            )
            metrics = compute_library_size_metrics(
                [sample], shutil.which("samtools"), config, threads=1
            )

        self.assertEqual(int(metrics.loc["tiny", "selected_library_size"]), 4)

    def test_multibamsummary_fragment_command_matches_library_filters(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            bed = tmp / "regions.bed"
            bed.write_text("chr1\t90\t210\tpeak1\n")
            bam = tmp / "tiny.bam"
            bam.touch()
            sample = SampleEntry("tiny", "A", bam, is_paired=True)
            config = CountFilterConfig(
                count_unit="fragment",
                min_mapq=30,
                exclude_duplicates=True,
            )
            commands: list[list[str]] = []

            def fake_run(cmd, **kwargs):
                raw = Path(cmd[cmd.index("--outRawCounts") + 1])
                raw.parent.mkdir(parents=True, exist_ok=True)
                raw.write_text("#chr\tstart\tend\ttiny\nchr1\t90\t210\t1\n")

            with patch("chipdiff.run_command", side_effect=fake_run):
                run_multibamsummary(
                    bed,
                    [sample],
                    tmp / "counts",
                    count_config=config,
                    command_log=commands,
                )

        cmd = commands[0]
        self.assertEqual(cmd[cmd.index("--samFlagInclude") + 1], "66")
        self.assertEqual(cmd[cmd.index("--samFlagExclude") + 1], "3844")
        self.assertEqual(cmd[cmd.index("--minMappingQuality") + 1], "30")
        self.assertIn("--extendReads", cmd)


if __name__ == "__main__":
    unittest.main()
