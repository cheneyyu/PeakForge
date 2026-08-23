from __future__ import annotations

import tempfile
import unittest
from argparse import Namespace
from pathlib import Path
from unittest.mock import patch

from chipdiff import _macs2_command, build_runmode_samples, load_samples


class ChipdiffControlTests(unittest.TestCase):
    def test_load_samples_reads_control_bam_column(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            bam = tmp / "sample.bam"
            control = tmp / "input.bam"
            bam.touch()
            control.touch()
            sheet = tmp / "samples.tsv"
            sheet.write_text(
                "sample\tcondition\tbam\tcontrol_bam\n"
                f"s1\tA\t{bam}\t{control}\n"
            )

            samples = load_samples(sheet)

        self.assertEqual(len(samples), 1)
        self.assertEqual(samples[0].control_bam, control)

    def test_build_runmode_samples_broadcasts_single_control_per_condition(self) -> None:
        args = Namespace(
            condition_a="K562",
            a_bams=["/tmp/k562_rep1.bam", "/tmp/k562_rep2.bam"],
            a_peaks=None,
            a_control_bams=["/tmp/k562_input.bam"],
            condition_b="HepG2",
            b_bams=["/tmp/hepg2_rep1.bam", "/tmp/hepg2_rep2.bam"],
            b_peaks=None,
            b_control_bams=["/tmp/hepg2_input.bam"],
        )

        samples = build_runmode_samples(args)

        self.assertEqual(len(samples), 4)
        self.assertEqual([str(sample.control_bam) for sample in samples[:2]], ["/tmp/k562_input.bam"] * 2)
        self.assertEqual([str(sample.control_bam) for sample in samples[2:]], ["/tmp/hepg2_input.bam"] * 2)

    def test_macs_command_includes_control_bam(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            bam = tmp / "sample.bam"
            control = tmp / "input.bam"
            bam.touch()
            control.touch()

            sample = build_runmode_samples(
                Namespace(
                    condition_a="K562",
                    a_bams=[str(bam)],
                    a_peaks=None,
                    a_control_bams=[str(control)],
                    condition_b="HepG2",
                    b_bams=[str(tmp / "other.bam")],
                    b_peaks=None,
                    b_control_bams=None,
                )
            )[0]
            sample.is_paired = False
            with patch("chipdiff.get_macs_command", return_value="macs3"):
                cmd, _ = _macs2_command(
                    sample,
                    output_dir=tmp,
                    macs2_genome="hs",
                    macs2_qval=0.01,
                    peak_type="narrow",
                )

        self.assertIn("-c", cmd)
        self.assertIn(str(control), cmd)

    def test_macs_command_exposes_read_based_assay_configuration(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            bam = tmp / "atac.bam"
            bam.touch()
            sample = build_runmode_samples(
                Namespace(
                    condition_a="A",
                    a_bams=[str(bam)],
                    a_peaks=None,
                    a_control_bams=None,
                    condition_b="B",
                    b_bams=[str(tmp / "other.bam")],
                    b_peaks=None,
                    b_control_bams=None,
                )
            )[0]
            sample.is_paired = False
            with patch("chipdiff.get_macs_command", return_value="macs3"):
                cmd, _ = _macs2_command(
                    sample,
                    output_dir=tmp,
                    macs2_genome="hs",
                    macs2_qval=0.01,
                    peak_type="narrow",
                    macs_format="BAM",
                    macs_nomodel=True,
                    macs_shift=-100,
                    macs_extsize=200,
                )

        self.assertEqual(cmd[cmd.index("-f") + 1], "BAM")
        self.assertIn("--nomodel", cmd)
        self.assertEqual(cmd[cmd.index("--shift") + 1], "-100")
        self.assertEqual(cmd[cmd.index("--extsize") + 1], "200")

    def test_macs_bampe_rejects_read_shift_and_extension(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            bam = tmp / "atac.bam"
            bam.touch()
            sample = build_runmode_samples(
                Namespace(
                    condition_a="A",
                    a_bams=[str(bam)],
                    a_peaks=None,
                    a_control_bams=None,
                    condition_b="B",
                    b_bams=[str(tmp / "other.bam")],
                    b_peaks=None,
                    b_control_bams=None,
                )
            )[0]
            sample.is_paired = True
            with patch("chipdiff.get_macs_command", return_value="macs3"):
                with self.assertRaisesRegex(ValueError, "requires macs_shift"):
                    _macs2_command(
                        sample,
                        output_dir=tmp,
                        macs2_genome="hs",
                        macs2_qval=0.01,
                        peak_type="narrow",
                        macs_format="BAMPE",
                        macs_shift=-100,
                    )
                with self.assertRaisesRegex(ValueError, "observed fragment lengths"):
                    _macs2_command(
                        sample,
                        output_dir=tmp,
                        macs2_genome="hs",
                        macs2_qval=0.01,
                        peak_type="narrow",
                        macs_format="BAMPE",
                        macs_extsize=200,
                    )


if __name__ == "__main__":
    unittest.main()
