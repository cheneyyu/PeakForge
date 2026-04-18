#!/usr/bin/env python3
"""Summarize 3v3 leave-one-pair-out validation outputs."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize held-out 3v3 PeakForge validation folds")
    parser.add_argument("--results-root", required=True, help="Root directory containing fold outputs")
    parser.add_argument("--output-dir", required=True, help="Directory where summary files will be written")
    parser.add_argument("--alpha", type=float, default=0.05, help="Adjusted p-value cutoff for gold-standard hits")
    parser.add_argument("--lfc", type=float, default=0.0, help="Absolute log2FC cutoff for gold-standard hits")
    return parser.parse_args()


def clean_header(value: str | None) -> str:
    return (value or "").replace("#", "").replace("'", "").replace('"', "").strip()


def parse_float(value: str | None, default: float = 0.0) -> float:
    if value in (None, ""):
        return default
    try:
        return float(value)
    except ValueError:
        return default


def load_results(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows: dict[str, dict[str, str]] = {}
        for row in reader:
            cleaned = {clean_header(key): value for key, value in row.items()}
            peak = cleaned.get("Peak")
            if not peak:
                raise ValueError(f"Results file {path} is missing Peak identifiers")
            rows[peak] = cleaned
    return rows


def auc(labels: list[bool], scores: list[float]) -> float:
    paired = [(float(score), 1 if label else 0) for label, score in zip(labels, scores)]
    positives = sum(label for _, label in paired)
    negatives = len(paired) - positives
    if positives == 0 or negatives == 0:
        return float("nan")

    paired.sort(key=lambda item: item[0])
    rank_sum = 0.0
    rank = 1
    idx = 0
    while idx < len(paired):
        end = idx + 1
        while end < len(paired) and paired[end][0] == paired[idx][0]:
            end += 1
        avg_rank = (rank + (rank + (end - idx) - 1)) / 2.0
        positives_in_tie = sum(label for _, label in paired[idx:end])
        rank_sum += avg_rank * positives_in_tie
        rank += end - idx
        idx = end

    return (rank_sum - positives * (positives + 1) / 2.0) / (positives * negatives)


def pearsonr(x: list[float], y: list[float]) -> float:
    if len(x) != len(y) or not x:
        return float("nan")
    mean_x = statistics.fmean(x)
    mean_y = statistics.fmean(y)
    numerator = sum((a - mean_x) * (b - mean_y) for a, b in zip(x, y))
    denom_x = math.sqrt(sum((a - mean_x) ** 2 for a in x))
    denom_y = math.sqrt(sum((b - mean_y) ** 2 for b in y))
    if denom_x == 0.0 or denom_y == 0.0:
        return float("nan")
    return numerator / (denom_x * denom_y)


def rankdata(values: list[float]) -> list[float]:
    indexed = sorted(enumerate(values), key=lambda item: item[1])
    ranks = [0.0] * len(values)
    start = 0
    while start < len(indexed):
        end = start + 1
        while end < len(indexed) and indexed[end][1] == indexed[start][1]:
            end += 1
        avg_rank = (start + 1 + end) / 2.0
        for idx in range(start, end):
            ranks[indexed[idx][0]] = avg_rank
        start = end
    return ranks


def spearmanr(x: list[float], y: list[float]) -> float:
    if len(x) != len(y) or not x:
        return float("nan")
    return pearsonr(rankdata(x), rankdata(y))


def evaluate_fold(
    fold: int,
    gold_rows: dict[str, dict[str, str]],
    heldout_rows: dict[str, dict[str, str]],
    *,
    alpha: float,
    lfc: float,
) -> dict[str, float | int]:
    shared_peaks = sorted(set(gold_rows) & set(heldout_rows))
    if not shared_peaks:
        raise ValueError(f"Fold {fold} has no shared peaks between gold and held-out runs")

    gold_lfc = [parse_float(gold_rows[peak].get("log2FC")) for peak in shared_peaks]
    heldout_lfc = [parse_float(heldout_rows[peak].get("log2FC")) for peak in shared_peaks]

    labels = [
        parse_float(gold_rows[peak].get("padj"), 1.0) <= alpha
        and abs(parse_float(gold_rows[peak].get("log2FC"))) >= lfc
        for peak in shared_peaks
    ]
    gold_sig = sum(labels)
    heldout_sig = sum(parse_float(heldout_rows[peak].get("padj"), 1.0) <= alpha for peak in shared_peaks)
    sign_concordance = (
        sum((math.copysign(1.0, g) == math.copysign(1.0, h)) for g, h in zip(gold_lfc, heldout_lfc)) / len(shared_peaks)
        if shared_peaks
        else float("nan")
    )

    return {
        "fold": fold,
        "n_shared_peaks": len(shared_peaks),
        "gold_significant": gold_sig,
        "heldout_significant": heldout_sig,
        "pearson_log2fc": pearsonr(gold_lfc, heldout_lfc),
        "spearman_log2fc": spearmanr(gold_lfc, heldout_lfc),
        "sign_concordance": sign_concordance,
        "heldout_abs_lfc_auroc": auc(labels, [abs(value) for value in heldout_lfc]),
    }


def write_tsv(path: Path, rows: list[dict[str, float | int]]) -> None:
    if not rows:
        raise ValueError("No rows available to write")
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def mean_ignore_nan(values: list[float]) -> float:
    clean = [value for value in values if not math.isnan(value)]
    if not clean:
        return float("nan")
    return statistics.fmean(clean)


def main() -> None:
    args = parse_args()
    results_root = Path(args.results_root)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, float | int]] = []
    for fold in (1, 2, 3):
        gold_path = results_root / f"fold{fold}_gold_2v2" / "differential_results.tsv"
        heldout_path = results_root / f"fold{fold}_heldout_1v1" / "differential_results.tsv"
        if not gold_path.exists():
            raise FileNotFoundError(f"Missing fold {fold} gold results: {gold_path}")
        if not heldout_path.exists():
            raise FileNotFoundError(f"Missing fold {fold} held-out results: {heldout_path}")
        rows.append(
            evaluate_fold(
                fold,
                load_results(gold_path),
                load_results(heldout_path),
                alpha=args.alpha,
                lfc=args.lfc,
            )
        )

    write_tsv(output_dir / "lopo_summary.tsv", rows)

    summary = {
        "alpha": args.alpha,
        "lfc": args.lfc,
        "folds": rows,
        "mean_shared_peaks": mean_ignore_nan([float(row["n_shared_peaks"]) for row in rows]),
        "mean_gold_significant": mean_ignore_nan([float(row["gold_significant"]) for row in rows]),
        "mean_heldout_significant": mean_ignore_nan([float(row["heldout_significant"]) for row in rows]),
        "mean_pearson_log2fc": mean_ignore_nan([float(row["pearson_log2fc"]) for row in rows]),
        "mean_spearman_log2fc": mean_ignore_nan([float(row["spearman_log2fc"]) for row in rows]),
        "mean_sign_concordance": mean_ignore_nan([float(row["sign_concordance"]) for row in rows]),
        "mean_heldout_abs_lfc_auroc": mean_ignore_nan([float(row["heldout_abs_lfc_auroc"]) for row in rows]),
    }
    (output_dir / "lopo_summary.json").write_text(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
