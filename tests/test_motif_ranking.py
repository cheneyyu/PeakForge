from __future__ import annotations

import gzip
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from motif_ranking import (
    MOTIF_COLUMN_SUFFIX,
    build_ranked_peak_table,
    load_motif_presence,
    score_motif_families,
    score_motif_presence,
)


class MotifRankingTests(unittest.TestCase):
    def setUp(self) -> None:
        self.diff_res = pd.DataFrame(
            {
                "log2FC": [2.5, 1.8, 1.1, -1.2, -1.7, -2.3],
                "log2FC_shrunk": [2.4, 1.7, 1.0, -1.1, -1.6, -2.2],
                "pvalue": [1e-7, 1e-6, 1e-4, 1e-4, 1e-6, 1e-7],
                "padj": [1e-5, 1e-4, 1e-3, 1e-3, 1e-4, 1e-5],
            },
            index=["peak1", "peak2", "peak3", "peak4", "peak5", "peak6"],
        )
        self.ranked = build_ranked_peak_table(self.diff_res)
        self.presence = pd.DataFrame(
            {
                "TumorTF/exp1": [True, True, True, False, False, False],
                "NormalTF/exp1": [False, False, False, True, True, True],
                "SharedFamily/exp1": [True, True, False, False, False, False],
                "SharedFamily/exp2": [False, False, False, False, True, True],
            },
            index=["peak1", "peak2", "peak3", "peak4", "peak5", "peak6"],
        )

    def test_ranked_peaks_keep_positive_tail_on_top(self) -> None:
        ordered = self.ranked["Peak"].tolist()
        self.assertEqual(ordered[:3], ["peak1", "peak2", "peak3"])
        self.assertEqual(ordered[-3:], ["peak4", "peak5", "peak6"])

    def test_exact_motif_scoring_recovers_expected_direction(self) -> None:
        scored = score_motif_presence(
            self.ranked,
            self.presence,
            positive_condition="tumor",
            negative_condition="normal",
            min_peaks=2,
            max_fraction=0.95,
            top_peaks=2,
        ).set_index("motif_name")

        self.assertEqual(scored.loc["TumorTF/exp1", "direction"], "tumor")
        self.assertGreater(scored.loc["TumorTF/exp1", "delta_auc"], 0)
        self.assertEqual(scored.loc["NormalTF/exp1", "direction"], "normal")
        self.assertLess(scored.loc["NormalTF/exp1", "delta_auc"], 0)

    def test_family_scoring_unions_member_peak_sets(self) -> None:
        family = score_motif_families(
            self.ranked,
            self.presence,
            positive_condition="tumor",
            negative_condition="normal",
            min_peaks=2,
            max_fraction=0.95,
            top_peaks=2,
        ).set_index("motif_name")

        self.assertIn("SharedFamily", family.index)
        self.assertEqual(family.loc["SharedFamily", "member_motif_count"], 2)
        self.assertEqual(family.loc["SharedFamily", "n_motif_peaks"], 4)

    def test_load_motif_presence_reads_homer_counts(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            scan_path = Path(tmpdir) / "motif_scan.tsv.gz"
            header = [
                "PeakID (cmd=annotatePeaks.pl fake.bed fake.fa)",
                f"TumorTF/exp1{MOTIF_COLUMN_SUFFIX}",
                f"NormalTF/exp1{MOTIF_COLUMN_SUFFIX}",
            ]
            rows = [
                ["peak1", "2", "0"],
                ["peak2", "0", "3"],
            ]
            with gzip.open(scan_path, "wt") as handle:
                handle.write("\t".join(header) + "\n")
                for row in rows:
                    handle.write("\t".join(row) + "\n")

            loaded = load_motif_presence(scan_path)

        self.assertEqual(loaded.index.tolist(), ["peak1", "peak2"])
        self.assertEqual(loaded.columns.tolist(), ["TumorTF/exp1", "NormalTF/exp1"])
        self.assertTrue(bool(loaded.loc["peak1", "TumorTF/exp1"]))
        self.assertFalse(bool(loaded.loc["peak1", "NormalTF/exp1"]))
        self.assertFalse(bool(loaded.loc["peak2", "TumorTF/exp1"]))
        self.assertTrue(bool(loaded.loc["peak2", "NormalTF/exp1"]))

    def test_single_pair_auto_ranking_uses_effect_not_sampling_pvalue(self) -> None:
        single_pair = pd.DataFrame(
            {
                "normalized_log2fc": [2.0, 1.0, -3.0],
                "mars_score": [1.0, 100.0, -2.0],
                "sampling_pvalue": [0.9, 1e-20, 0.8],
                "sampling_qvalue": [0.9, 3e-20, 0.9],
                "mean_cpm": [10.0, 10.0, 10.0],
                "rank_eligible": [True, True, True],
                "analysis_mode": ["single_pair_exploratory"] * 3,
            },
            index=["effect_first", "tiny_pvalue", "negative_tail"],
        )

        ranked = build_ranked_peak_table(single_pair, score_metric="auto")

        self.assertEqual(ranked.iloc[0]["Peak"], "effect_first")
        self.assertEqual(ranked.iloc[0]["motif_rank_metric"], "signed_normalized_log2fc")
        with self.assertRaisesRegex(ValueError, "cannot be used"):
            build_ranked_peak_table(single_pair, score_metric="signed_product")


if __name__ == "__main__":
    unittest.main()
