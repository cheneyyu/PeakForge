#!/usr/bin/env python3
"""Summarize fixed-consensus held-out ranking comparisons.

Each fold compares an exploratory exact-1-vs-1 ranking with formal calls from
a replicate-supported 2-vs-2 reference. The reference is empirical, not
independent ground truth.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize three fixed-consensus held-out PeakForge folds"
    )
    parser.add_argument(
        "--results-root", required=True, help="Root directory containing fold outputs"
    )
    parser.add_argument(
        "--output-dir", required=True, help="Directory for summary tables"
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.05,
        help="Adjusted-p cutoff for replicate-supported reference calls",
    )
    parser.add_argument(
        "--lfc",
        type=float,
        default=0.0,
        help="Minimum absolute reference log2 fold change",
    )
    parser.add_argument(
        "--top-fraction",
        type=float,
        default=0.05,
        help="Fraction of eligible held-out peaks used for top-rank precision/recall",
    )
    args = parser.parse_args()
    if not 0.0 < args.alpha <= 1.0:
        parser.error("--alpha must be in (0, 1]")
    if args.lfc < 0.0:
        parser.error("--lfc must be non-negative")
    if not 0.0 < args.top_fraction <= 1.0:
        parser.error("--top-fraction must be in (0, 1]")
    return args


def clean_header(value: object) -> str:
    return str(value).replace("#", "").replace("'", "").replace('"', "").strip()


def load_results(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(path)
    frame = pd.read_csv(path, sep="\t")
    frame.columns = [clean_header(column) for column in frame.columns]
    if "Peak" not in frame.columns:
        first = frame.columns[0]
        if first.startswith("Unnamed"):
            frame = frame.rename(columns={first: "Peak"})
    if "Peak" not in frame.columns:
        raise ValueError(f"Results file {path} is missing Peak identifiers")
    frame["Peak"] = frame["Peak"].astype(str)
    if frame["Peak"].duplicated().any():
        raise ValueError(f"Results file {path} contains duplicate Peak identifiers")
    return frame.set_index("Peak", drop=False)


def boolean_series(values: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(values):
        return values.fillna(False)
    return values.astype(str).str.lower().isin({"true", "1", "yes"})


def auroc(labels: np.ndarray, scores: np.ndarray) -> float:
    positives = int(labels.sum())
    negatives = int(len(labels) - positives)
    if positives == 0 or negatives == 0:
        return float("nan")
    ranks = rankdata(scores, method="average")
    rank_sum = float(ranks[labels].sum())
    return (rank_sum - positives * (positives + 1) / 2.0) / (
        positives * negatives
    )


def average_precision(labels: np.ndarray, scores: np.ndarray) -> float:
    positives = int(labels.sum())
    if positives == 0:
        return float("nan")
    order = np.argsort(-scores, kind="mergesort")
    ordered_scores = scores[order]
    ordered_labels = labels[order]
    true_positives = 0
    false_positives = 0
    previous_recall = 0.0
    result = 0.0
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and ordered_scores[end] == ordered_scores[start]:
            end += 1
        group = ordered_labels[start:end]
        true_positives += int(group.sum())
        false_positives += int(len(group) - group.sum())
        recall = true_positives / positives
        precision = true_positives / (true_positives + false_positives)
        result += (recall - previous_recall) * precision
        previous_recall = recall
        start = end
    return result


def precision_at(labels: np.ndarray, scores: np.ndarray, size: int) -> float:
    if len(labels) == 0:
        return float("nan")
    selected = np.argsort(-scores, kind="mergesort")[: min(size, len(labels))]
    return float(labels[selected].mean())


def json_ready(value: object) -> object:
    if isinstance(value, dict):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, list):
        return [json_ready(item) for item in value]
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def evaluate_fold(
    fold: int,
    reference: pd.DataFrame,
    heldout: pd.DataFrame,
    *,
    alpha: float,
    lfc: float,
    top_fraction: float,
) -> dict[str, object]:
    reference_required = {"Peak", "log2FC", "padj"}
    heldout_required = {"Peak", "normalized_log2fc", "rank_eligible"}
    if missing := reference_required - set(reference.columns):
        raise ValueError(f"Reference fold {fold} is missing: {sorted(missing)}")
    if missing := heldout_required - set(heldout.columns):
        raise ValueError(f"Held-out fold {fold} is missing: {sorted(missing)}")

    shared = reference.index.intersection(heldout.index, sort=False)
    if shared.empty:
        raise ValueError(f"Fold {fold} has no shared fixed-consensus intervals")
    frame = pd.DataFrame(
        {
            "reference_effect": pd.to_numeric(
                reference.loc[shared, "log2FC"], errors="coerce"
            ),
            "reference_padj": pd.to_numeric(
                reference.loc[shared, "padj"], errors="coerce"
            ),
            "heldout_effect": pd.to_numeric(
                heldout.loc[shared, "normalized_log2fc"], errors="coerce"
            ),
            "rank_eligible": boolean_series(
                heldout.loc[shared, "rank_eligible"]
            ).to_numpy(),
        },
        index=shared,
    )
    frame["reference_positive"] = frame["reference_padj"].lt(alpha) & frame[
        "reference_effect"
    ].abs().ge(lfc)
    evaluable = frame.loc[
        frame["rank_eligible"]
        & frame["reference_padj"].notna()
        & frame["reference_effect"].notna()
        & frame["heldout_effect"].notna()
    ].copy()
    if evaluable.empty:
        raise ValueError(f"Fold {fold} has no rank-eligible shared intervals")

    labels = evaluable["reference_positive"].to_numpy(dtype=bool)
    scores = evaluable["heldout_effect"].abs().to_numpy(dtype=float)
    top_size = max(1, math.ceil(len(evaluable) * top_fraction))
    top_order = np.argsort(-scores, kind="mergesort")[:top_size]
    top_true = int(labels[top_order].sum())
    total_positive = int(labels.sum())

    return {
        "fold": fold,
        "shared_intervals": int(len(frame)),
        "evaluation_intervals": int(len(evaluable)),
        "reference_formal_calls_shared": int(frame["reference_positive"].sum()),
        "reference_formal_calls_evaluable": total_positive,
        "positive_prevalence": total_positive / len(evaluable),
        "heldout_rank_eligible_intervals": int(frame["rank_eligible"].sum()),
        "pearson_effect": float(
            evaluable["reference_effect"].corr(
                evaluable["heldout_effect"], method="pearson"
            )
        ),
        "spearman_effect": float(
            evaluable["reference_effect"].corr(
                evaluable["heldout_effect"], method="spearman"
            )
        ),
        "sign_concordance": float(
            (
                np.sign(evaluable["reference_effect"])
                == np.sign(evaluable["heldout_effect"])
            ).mean()
        ),
        "auroc": auroc(labels, scores),
        "auprc": average_precision(labels, scores),
        "precision_at_100": precision_at(labels, scores, 100),
        "precision_at_500": precision_at(labels, scores, 500),
        "top_fraction": top_fraction,
        "precision_at_top_fraction": top_true / top_size,
        "recall_at_top_fraction": (
            top_true / total_positive if total_positive else float("nan")
        ),
        "analysis_interpretation": (
            "exploratory_ranking_vs_replicate_supported_reference_not_ground_truth"
        ),
    }


def main() -> None:
    args = parse_args()
    results_root = Path(args.results_root)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for fold in (1, 2, 3):
        reference_path = (
            results_root
            / f"fold{fold}_reference_2v2"
            / "differential_results.tsv"
        )
        heldout_path = (
            results_root / f"fold{fold}_heldout_1v1" / "differential_results.tsv"
        )
        rows.append(
            evaluate_fold(
                fold,
                load_results(reference_path),
                load_results(heldout_path),
                alpha=args.alpha,
                lfc=args.lfc,
                top_fraction=args.top_fraction,
            )
        )

    frame = pd.DataFrame(rows)
    frame.to_csv(output_dir / "lopo_summary.tsv", sep="\t", index=False)
    numeric = frame.select_dtypes(include=[np.number])
    summary = {
        "reference_status": "replicate_supported_reference_not_ground_truth",
        "single_pair_status": "exploratory_ranking_no_biological_variance_estimation",
        "alpha": args.alpha,
        "lfc": args.lfc,
        "top_fraction": args.top_fraction,
        "folds": rows,
        "mean_metrics": numeric.mean(numeric_only=True).to_dict(),
    }
    (output_dir / "lopo_summary.json").write_text(
        json.dumps(json_ready(summary), indent=2, allow_nan=False) + "\n"
    )


if __name__ == "__main__":
    main()
