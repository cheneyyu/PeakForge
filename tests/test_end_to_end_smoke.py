from __future__ import annotations

import os
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd
import pysam


def _write_bam(path: Path, starts: list[int]) -> None:
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    )
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for idx, start in enumerate(starts):
            name = f"fragment_{idx}"
            for flag, read_start, mate_start, tlen in (
                (99, start, start + 50, 100),
                (147, start + 50, start, -100),
            ):
                read = pysam.AlignedSegment(header)
                read.query_name = name
                read.query_sequence = "A" * 50
                read.flag = flag
                read.reference_id = 0
                read.reference_start = read_start
                read.mapping_quality = 60
                read.cigar = ((0, 50),)
                read.next_reference_id = 0
                read.next_reference_start = mate_start
                read.template_length = tlen
                read.query_qualities = pysam.qualitystring_to_array("I" * 50)
                out.write(read)
    pysam.index(str(path))


@unittest.skipUnless(
    shutil.which("multiBamSummary") and shutil.which("samtools"),
    "deepTools and samtools are required",
)
class EndToEndSmokeTests(unittest.TestCase):
    def test_single_pair_fragment_pipeline_is_thread_invariant(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            bam_a = tmp / "a.bam"
            bam_b = tmp / "b.bam"
            _write_bam(bam_a, [100, 200])
            _write_bam(bam_b, [100, 300])
            consensus = tmp / "consensus.bed"
            consensus.write_text(
                "chr1\t90\t180\tp1\n"
                "chr1\t190\t280\tp2\n"
                "chr1\t290\t380\tp3\n"
            )
            outputs = []
            script = Path(__file__).resolve().parents[1] / "chipdiff.py"
            for threads in (1, 2):
                output = tmp / f"out_t{threads}"
                cmd = [
                    sys.executable,
                    str(script),
                    "runmode",
                    "--condition-a",
                    "A",
                    "--a-bams",
                    str(bam_a),
                    "--condition-b",
                    "B",
                    "--b-bams",
                    str(bam_b),
                    "--consensus-peaks",
                    str(consensus),
                    "--count-unit",
                    "fragment",
                    "--min-mapq",
                    "30",
                    "--exclude-duplicates",
                    "--threads",
                    str(threads),
                    "--output-dir",
                    str(output),
                ]
                completed = subprocess.run(
                    cmd,
                    check=False,
                    capture_output=True,
                    text=True,
                    env=os.environ.copy(),
                )
                self.assertEqual(completed.returncode, 0, completed.stderr)
                outputs.append(output)

            count_frames = [
                pd.read_csv(path / "counts" / "counts.tsv", sep="\t") for path in outputs
            ]
            result_frames = [
                pd.read_csv(path / "differential_results.tsv", sep="\t") for path in outputs
            ]
            pd.testing.assert_frame_equal(count_frames[0], count_frames[1])
            pd.testing.assert_frame_equal(result_frames[0], result_frames[1])
            self.assertNotIn("padj", result_frames[0].columns)
            self.assertTrue(
                (result_frames[0]["analysis_mode"] == "single_pair_exploratory").all()
            )
            metadata = json.loads((outputs[0] / "metadata.json").read_text())
            self.assertEqual(metadata["peakforge_version"], "0.2.3")
            self.assertIn("commit", metadata["git"])
            self.assertEqual(metadata["count_filter"]["count_unit"], "fragment")
            self.assertTrue(metadata["executed_commands"])
            self.assertIn("samtools", metadata["software_versions"]["external_tools"])


if __name__ == "__main__":
    unittest.main()
