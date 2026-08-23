from __future__ import annotations

import math
import unittest

import pandas as pd

from chipdiff import mars_differential


class MarsDifferentialTests(unittest.TestCase):
    def test_mars_reports_library_normalized_log2fc(self) -> None:
        counts = pd.DataFrame(
            {
                "pc": [50, 100],
                "pn": [50, 400],
            },
            index=["peak_equal_raw_counts", "peak_depth_matched_counts"],
        )
        conditions = pd.Series({"pc": "PC", "pn": "PN"})
        library_sizes = pd.Series({"pc": 1_000.0, "pn": 4_000.0})

        result = mars_differential(counts, conditions, library_sizes)

        expected_equal_raw = math.log2((50.5 / 4_000.0) / (50.5 / 1_000.0))
        expected_depth_matched = math.log2((400.5 / 4_000.0) / (100.5 / 1_000.0))

        self.assertAlmostEqual(result.loc["peak_equal_raw_counts", "log2FC"], expected_equal_raw, places=10)
        self.assertAlmostEqual(result.loc["peak_equal_raw_counts", "log2FC_shrunk"], expected_equal_raw, places=10)
        self.assertAlmostEqual(result.loc["peak_depth_matched_counts", "log2FC"], expected_depth_matched, places=10)
        self.assertAlmostEqual(result.loc["peak_depth_matched_counts", "log2FC_shrunk"], expected_depth_matched, places=10)

        # If counts scale with library size, the normalized effect should be ~0.
        self.assertAlmostEqual(result.loc["peak_depth_matched_counts", "log2FC"], 0.0, delta=0.01)


if __name__ == "__main__":
    unittest.main()
