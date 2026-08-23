from __future__ import annotations

import math
import unittest

import pandas as pd

from chipdiff import call_differential_analysis, single_pair_statistics


class SinglePairStatisticsTests(unittest.TestCase):
    def setUp(self) -> None:
        self.counts = pd.DataFrame(
            {
                "sample_a": [0, 0, 10, 100, 50],
                "sample_b": [0, 8, 20, 400, 50],
            },
            index=["both_zero", "one_zero", "positive", "depth_matched", "equal_raw"],
        )
        self.conditions = pd.Series({"sample_a": "A", "sample_b": "B"})
        self.library_sizes = pd.Series({"sample_a": 1_000.0, "sample_b": 4_000.0})

    def test_canonical_formula_and_schema(self) -> None:
        result = single_pair_statistics(
            self.counts,
            self.conditions,
            self.library_sizes,
            minimum_mean_cpm=1.0,
            pseudocount=0.5,
        )

        expected = math.log2((400.5 / 4_000.0) / (100.5 / 1_000.0))
        self.assertAlmostEqual(
            result.loc["depth_matched", "normalized_log2fc"], expected, places=12
        )
        self.assertNotIn("pvalue", result.columns)
        self.assertNotIn("padj", result.columns)
        self.assertIn("sampling_pvalue", result.columns)
        self.assertIn("sampling_qvalue", result.columns)
        self.assertTrue((result["analysis_mode"] == "single_pair_exploratory").all())
        self.assertTrue(
            (
                result["interpretation"]
                == "ranking_only_no_biological_variance_estimation"
            ).all()
        )

    def test_zero_count_policy(self) -> None:
        result = single_pair_statistics(
            self.counts,
            self.conditions,
            self.library_sizes,
            minimum_mean_cpm=1.0,
        )

        self.assertEqual(result.loc["both_zero", "zero_count_status"], "both_zero")
        self.assertFalse(bool(result.loc["both_zero", "rank_eligible"]))
        self.assertTrue(pd.isna(result.loc["both_zero", "mars_score"]))
        self.assertEqual(result.loc["both_zero", "sampling_pvalue"], 1.0)
        self.assertEqual(result.loc["one_zero", "zero_count_status"], "one_zero")
        self.assertTrue(bool(result.loc["one_zero", "rank_eligible"]))
        self.assertTrue(pd.isna(result.loc["one_zero", "mars_score"]))

    def test_primary_rank_uses_normalized_effect(self) -> None:
        result = single_pair_statistics(
            self.counts,
            self.conditions,
            self.library_sizes,
            minimum_mean_cpm=1.0,
        )
        up = result.dropna(subset=["rank_up"]).sort_values("rank_up")
        down = result.dropna(subset=["rank_down"]).sort_values("rank_down")
        self.assertEqual(up.index[0], "one_zero")
        self.assertEqual(down.index[0], "equal_raw")

    def test_tied_ranks_use_peak_identifier_as_stable_tie_break(self) -> None:
        counts = pd.DataFrame(
            {"sample_a": [10, 10, 10], "sample_b": [20, 20, 20]},
            index=["peak_b", "peak_a", "peak_c"],
        )
        conditions = pd.Series({"sample_a": "A", "sample_b": "B"})
        library_sizes = pd.Series({"sample_a": 1_000.0, "sample_b": 1_000.0})

        forward = single_pair_statistics(counts, conditions, library_sizes)
        reverse = single_pair_statistics(
            counts.iloc[::-1], conditions, library_sizes
        )
        rank_columns = ["rank_up", "rank_down", "rank_absolute_lfc"]
        pd.testing.assert_frame_equal(
            forward[rank_columns].sort_index(),
            reverse[rank_columns].sort_index(),
        )
        self.assertEqual(int(forward.loc["peak_a", "rank_up"]), 1)

    def test_mixed_group_sizes_are_rejected(self) -> None:
        counts = pd.DataFrame(
            {"a1": [10], "a2": [11], "b1": [20]},
            index=["peak"],
        )
        conditions = pd.Series({"a1": "A", "a2": "A", "b1": "B"})
        library_sizes = pd.Series({"a1": 1000, "a2": 1000, "b1": 1000})
        with self.assertRaisesRegex(ValueError, "exact 1-vs-1"):
            call_differential_analysis(counts, conditions, library_sizes)


if __name__ == "__main__":
    unittest.main()
