from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd
import pyranges as pr

from chipdiff import build_consensus, read_peak_file


class PeakHandlingTests(unittest.TestCase):
    def test_edge_padding_extends_original_boundaries_and_clips_zero(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            peaks = Path(tmpdir) / "sample.narrowPeak"
            peaks.write_text("chr1\t50\t150\tp1\t100\t.\t1\t1\t1\t25\n")
            padded = read_peak_file(peaks, "narrow", 100).df.iloc[0]

        self.assertEqual(int(padded.Start), 0)
        self.assertEqual(int(padded.End), 250)

    def test_broad_peak_boundaries_are_not_padded(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            peaks = Path(tmpdir) / "sample.broadPeak"
            peaks.write_text("chr1\t50\t150\tp1\t100\t.\t1\t1\t1\n")
            loaded = read_peak_file(peaks, "broad", 100).df.iloc[0]

        self.assertEqual(int(loaded.Start), 50)
        self.assertEqual(int(loaded.End), 150)

    def test_summit_fixed_uses_narrowpeak_offset_and_exact_total_width(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            peaks = Path(tmpdir) / "sample.narrowPeak"
            peaks.write_text("chr1\t50\t150\tp1\t100\t.\t1\t1\t1\t25\n")
            loaded = read_peak_file(
                peaks,
                "narrow",
                250,
                peak_coordinate_mode="summit-fixed",
                summit_fixed_width=101,
            ).df.iloc[0]

        self.assertEqual((int(loaded.Start), int(loaded.End)), (25, 126))
        self.assertEqual(int(loaded.End) - int(loaded.Start), 101)

    def test_summit_fixed_clips_start_and_preserves_requested_width(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            peaks = Path(tmpdir) / "sample.narrowPeak"
            peaks.write_text("chr1\t0\t50\tp1\t100\t.\t1\t1\t1\t5\n")
            loaded = read_peak_file(
                peaks,
                "narrow",
                0,
                peak_coordinate_mode="summit-fixed",
                summit_fixed_width=100,
            ).df.iloc[0]

        self.assertEqual((int(loaded.Start), int(loaded.End)), (0, 100))

    def test_summit_fixed_rejects_missing_or_invalid_summits(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            missing = Path(tmpdir) / "missing.narrowPeak"
            missing.write_text("chr1\t50\t150\n")
            with self.assertRaisesRegex(ValueError, "at least 10 columns"):
                read_peak_file(
                    missing,
                    "narrow",
                    0,
                    peak_coordinate_mode="summit-fixed",
                )

            invalid = Path(tmpdir) / "invalid.narrowPeak"
            invalid.write_text("chr1\t50\t150\tp1\t100\t.\t1\t1\t1\t-1\n")
            with self.assertRaisesRegex(ValueError, "summit offset"):
                read_peak_file(
                    invalid,
                    "narrow",
                    0,
                    peak_coordinate_mode="summit-fixed",
                )

    def test_summit_fixed_rejects_broad_peaks(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            peaks = Path(tmpdir) / "sample.broadPeak"
            peaks.write_text("chr1\t50\t150\tp1\t100\t.\t1\t1\t1\n")
            with self.assertRaisesRegex(ValueError, "require narrowPeak"):
                read_peak_file(
                    peaks,
                    "broad",
                    0,
                    peak_coordinate_mode="summit-fixed",
                )

    def test_consensus_support_counts_distinct_samples(self) -> None:
        a = pr.PyRanges(
            pd.DataFrame(
                {"Chromosome": ["chr1"], "Start": [100], "End": [200], "Sample": ["a"]}
            )
        )
        b = pr.PyRanges(
            pd.DataFrame(
                {"Chromosome": ["chr1"], "Start": [150], "End": [250], "Sample": ["b"]}
            )
        )
        consensus = build_consensus({"a": a, "b": b}, min_overlap=2).df
        self.assertEqual(len(consensus), 1)
        self.assertEqual(int(consensus.iloc[0].Support), 2)
        self.assertEqual((int(consensus.iloc[0].Start), int(consensus.iloc[0].End)), (100, 250))


if __name__ == "__main__":
    unittest.main()
