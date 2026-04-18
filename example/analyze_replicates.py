#!/usr/bin/env python3
"""Assess reproducibility between 2v2 and 1v1 PeakForge runs."""
from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd
import pyranges as pr
from scipy import stats


@dataclass
class ResultBundle:
    label: str
    path: Path
    full: pd.DataFrame
    significant: pd.DataFrame


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize replicate concordance for PeakForge outputs")
    parser.add_argument("--main", required=True, help="Differential results from the 2v2 run")
    parser.add_argument("--one-vs-one", nargs="+", required=True,
                        help="Differential results from 1v1 runs (at least two files)")
    parser.add_argument("--output-dir", required=True, help="Directory where reports will be written")
    parser.add_argument("--alpha", type=float, default=0.05, help="Adjusted p-value cutoff")
    parser.add_argument("--lfc", type=float, default=1.0, help="Absolute log2 fold-change cutoff")
    parser.add_argument("--top-n", type=int, default=200,
                        help="Top-N peaks (by padj) from the 2v2 run to evaluate overlap")
    return parser.parse_args()


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def load_results(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Results file not found: {path}")
    df = pd.read_csv(path, sep="\t", index_col=0)
    required = {"log2FC", "padj"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Results file {path} missing required columns: {', '.join(sorted(missing))}")
    df.index = df.index.astype(str)
    return df


def filter_significant(df: pd.DataFrame, alpha: float, lfc: float) -> pd.DataFrame:
    mask = (df["padj"] <= alpha) & (df["log2FC"].abs() >= lfc)
    return df.loc[mask].copy()


def parse_peak(name: str) -> Tuple[str, int, int]:
    chrom, coords = name.split(":", 1)
    start_str, end_str = coords.split("-", 1)
    return chrom, int(start_str.replace(",", "")), int(end_str.replace(",", ""))


def peaks_to_ranges(df: pd.DataFrame) -> pr.PyRanges:
    if df.empty:
        return pr.PyRanges()
    records = []
    for peak in df.index:
        try:
            chrom, start, end = parse_peak(peak)
        except Exception as exc:
            raise ValueError(f"Failed to parse peak name '{peak}': {exc}") from exc
        records.append({
            "Chromosome": chrom,
            "Start": int(start),
            "End": int(end),
            "Peak": peak,
        })
    return pr.PyRanges(pd.DataFrame.from_records(records))


def overlap_pairs(main: pd.DataFrame, sub: pd.DataFrame) -> pd.DataFrame:
    main_ranges = peaks_to_ranges(main)
    sub_ranges = peaks_to_ranges(sub)
    if len(main_ranges) == 0 or len(sub_ranges) == 0:
        return pd.DataFrame(columns=["Peak_main", "Peak_sub"])
    joined = main_ranges.join(sub_ranges, suffix="_sub")
    if joined.df.empty:
        return pd.DataFrame(columns=["Peak_main", "Peak_sub"])
    cols = [col for col in joined.df.columns if col.startswith("Peak")]
    mapping = joined.df[[cols[0], cols[1]]].drop_duplicates()
    mapping.columns = ["Peak_main", "Peak_sub"]
    return mapping


def bp_length(ranges: pr.PyRanges) -> int:
    if len(ranges) == 0:
        return 0
    merged = ranges.merge()
    if merged.df.empty:
        return 0
    lengths = (merged.df["End"] - merged.df["Start"]).astype(int)
    return int(lengths.sum())


def evaluate_pair(main_bundle: ResultBundle, sub_bundle: ResultBundle, top_n: int) -> Dict[str, float]:
    pairs = overlap_pairs(main_bundle.significant, sub_bundle.significant)
    shared_main = set(pairs["Peak_main"]) if not pairs.empty else set()
    shared_sub = set(pairs["Peak_sub"]) if not pairs.empty else set()

    main_count = len(main_bundle.significant)
    sub_count = len(sub_bundle.significant)
    overlap_count = len(shared_main)
    precision = overlap_count / sub_count if sub_count else 0.0
    recall = overlap_count / main_count if main_count else 0.0

    concordance = float("nan")
    spearman_r = float("nan")
    spearman_p = float("nan")
    mean_abs_diff = float("nan")
    if overlap_count:
        merged = pairs.drop_duplicates()
        merged["log2FC_main"] = main_bundle.significant.loc[merged["Peak_main"], "log2FC"].to_numpy()
        merged["log2FC_sub"] = sub_bundle.significant.loc[merged["Peak_sub"], "log2FC"].to_numpy()
        concordance = float((np.sign(merged["log2FC_main"]) == np.sign(merged["log2FC_sub"])).mean())
        if overlap_count > 1:
            spearman_r, spearman_p = stats.spearmanr(merged["log2FC_main"], merged["log2FC_sub"])
        mean_abs_diff = float(np.mean(np.abs(merged["log2FC_main"] - merged["log2FC_sub"])))

    main_ranges = peaks_to_ranges(main_bundle.significant)
    sub_ranges = peaks_to_ranges(sub_bundle.significant)
    overlap_bp = bp_length(main_ranges.intersect(sub_ranges)) if len(main_ranges) and len(sub_ranges) else 0
    union_ranges = pr.PyRanges(pd.concat([main_ranges.df, sub_ranges.df], ignore_index=True)) if overlap_bp else None
    union_bp = bp_length(union_ranges) if union_ranges is not None else (
        bp_length(main_ranges) + bp_length(sub_ranges)
    )
    jaccard = overlap_bp / union_bp if union_bp else 0.0

    top_n = min(top_n, len(main_bundle.significant))
    top_overlap_fraction = 0.0
    if top_n:
        top_main = main_bundle.significant.sort_values("padj").head(top_n)
        top_pairs = overlap_pairs(top_main, sub_bundle.significant)
        unique_top = set(top_pairs["Peak_main"]) if not top_pairs.empty else set()
        top_overlap_fraction = len(unique_top) / top_n if top_n else 0.0

    return {
        "comparison": sub_bundle.label,
        "n_significant_main": main_count,
        "n_significant_1v1": sub_count,
        "n_overlap": overlap_count,
        "precision": precision,
        "recall": recall,
        "f1": (2 * precision * recall / (precision + recall)) if (precision + recall) else 0.0,
        "sign_concordance": concordance,
        "spearman_r": spearman_r,
        "spearman_p": spearman_p,
        "mean_abs_log2fc_diff": mean_abs_diff,
        "overlap_bp": overlap_bp,
        "union_bp": union_bp,
        "jaccard_bp": jaccard,
        "top_n": top_n,
        "top_overlap_fraction": top_overlap_fraction,
    }


def load_bundle(label: str, path: Path, alpha: float, lfc: float) -> ResultBundle:
    full = load_results(path)
    significant = filter_significant(full, alpha, lfc)
    return ResultBundle(label=label, path=path, full=full, significant=significant)


def main() -> None:
    args = parse_args()
    output_dir = ensure_dir(Path(args.output_dir))

    main_bundle = load_bundle("2v2", Path(args.main), args.alpha, args.lfc)
    comparisons: List[ResultBundle] = []
    for path in args.one_vs_one:
        label = Path(path).parent.name
        comparisons.append(load_bundle(label, Path(path), args.alpha, args.lfc))

    summary_rows = [evaluate_pair(main_bundle, bundle, args.top_n) for bundle in comparisons]
    summary_df = pd.DataFrame(summary_rows)
    summary_path = output_dir / "replicate_overlap_summary.tsv"
    summary_df.to_csv(summary_path, sep="\t", index=False)

    review = {
        "alpha": args.alpha,
        "lfc": args.lfc,
        "top_n": args.top_n,
        "main_results": str(main_bundle.path),
        "comparisons": [row for row in summary_rows],
    }
    (output_dir / "replicate_overlap_summary.json").write_text(json.dumps(review, indent=2))

    global_matrix = build_global_matrix(main_bundle, comparisons)
    global_matrix.to_csv(output_dir / "global_log2fc_matrix.tsv", sep="\t")


def build_global_matrix(main_bundle: ResultBundle, comparisons: Iterable[ResultBundle]) -> pd.DataFrame:
    frames = {main_bundle.label: main_bundle.full["log2FC"]}
    for bundle in comparisons:
        frames[bundle.label] = bundle.full["log2FC"]
    combined = pd.concat(frames, axis=1)
    combined.index.name = "Peak"
    return combined


if __name__ == "__main__":
    main()
