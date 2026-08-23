from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd


def load_module():
    path = (
        Path(__file__).resolve().parents[1]
        / "example"
        / "summarize_lopo_validation.py"
    )
    spec = importlib.util.spec_from_file_location("summarize_lopo_validation", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_lopo_uses_single_pair_ranking_without_sampling_qvalue() -> None:
    module = load_module()
    reference = pd.DataFrame(
        {
            "Peak": ["p1", "p2", "p3", "p4"],
            "log2FC": [2.0, -2.0, 0.1, -0.1],
            "padj": [0.01, 0.02, 0.8, 0.9],
        }
    ).set_index("Peak", drop=False)
    heldout = pd.DataFrame(
        {
            "Peak": ["p1", "p2", "p3", "p4"],
            "normalized_log2fc": [3.0, -2.5, 0.2, -0.1],
            "rank_eligible": [True, True, True, True],
        }
    ).set_index("Peak", drop=False)

    metrics = module.evaluate_fold(
        1,
        reference,
        heldout,
        alpha=0.05,
        lfc=0.0,
        top_fraction=0.5,
    )

    assert metrics["auroc"] == 1.0
    assert metrics["auprc"] == 1.0
    assert metrics["precision_at_top_fraction"] == 1.0
    assert metrics["recall_at_top_fraction"] == 1.0
    assert metrics["analysis_interpretation"].endswith("not_ground_truth")
