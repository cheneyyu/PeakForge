from __future__ import annotations

import gzip
import logging
import os
import shutil
import subprocess
from pathlib import Path
from typing import Dict, Iterable, Optional, Sequence

import numpy as np
import pandas as pd
from scipy import stats


MOTIF_COLUMN_SUFFIX = " Distance From Peak(sequence,strand,conservation)"


def _repo_root() -> Path:
    return Path(__file__).resolve().parent


def resolve_homer_tool(tool_name: str) -> str:
    """Return a HOMER executable path with pragmatic local fallbacks."""

    direct = shutil.which(tool_name)
    if direct:
        return direct

    candidates = []
    homer_home = os.environ.get("HOMER_HOME")
    if homer_home:
        home = Path(homer_home)
        candidates.extend(
            [
                home / "bin" / tool_name,
                home / tool_name,
                home / "share" / "homer" / "bin" / tool_name,
            ]
        )

    repo = _repo_root()
    candidates.append(repo / ".micromamba" / "envs" / "homer" / "share" / "homer" / "bin" / tool_name)

    for candidate in candidates:
        if candidate.exists():
            return str(candidate)

    raise RuntimeError(
        f"Unable to locate required HOMER executable '{tool_name}'. "
        "Install HOMER and place it on PATH, set HOMER_HOME, or use the bundled micromamba environment."
    )


def resolve_default_known_motif_file() -> Path:
    """Resolve HOMER's vertebrate known motif database."""

    candidates = []
    homer_home = os.environ.get("HOMER_HOME")
    if homer_home:
        home = Path(homer_home)
        candidates.extend(
            [
                home / "data" / "knownTFs" / "vertebrates" / "known.motifs",
                home / "knownTFs" / "vertebrates" / "known.motifs",
                home / "share" / "homer" / "data" / "knownTFs" / "vertebrates" / "known.motifs",
            ]
        )

    annotate_tool = Path(resolve_homer_tool("annotatePeaks.pl")).resolve()
    if len(annotate_tool.parents) >= 2:
        candidates.append(annotate_tool.parents[1] / "data" / "knownTFs" / "vertebrates" / "known.motifs")

    candidates.append(
        _repo_root()
        / ".micromamba"
        / "envs"
        / "homer"
        / "share"
        / "homer"
        / "data"
        / "knownTFs"
        / "vertebrates"
        / "known.motifs"
    )

    for candidate in candidates:
        if candidate.exists():
            return candidate

    raise RuntimeError(
        "Unable to resolve HOMER vertebrate known motif database. "
        "Supply --motif-file explicitly or install HOMER with the knownTFs database."
    )


def build_homer_runtime_env(annotate_tool: Path) -> Dict[str, str]:
    """Prepare an environment where annotatePeaks.pl can find HOMER helper binaries."""

    env = os.environ.copy()
    search_paths: list[str] = []

    def add_if_exists(path: Path) -> None:
        if path.exists():
            text = str(path)
            if text not in search_paths:
                search_paths.append(text)

    homer_home = os.environ.get("HOMER_HOME")
    if homer_home:
        home = Path(homer_home)
        add_if_exists(home)
        add_if_exists(home / "bin")
        add_if_exists(home / "share" / "homer" / "bin")

    add_if_exists(annotate_tool.parent)
    for parent in annotate_tool.parents:
        add_if_exists(parent / "bin")
        add_if_exists(parent / "share" / "homer" / "bin")

    existing_path = env.get("PATH", "")
    env["PATH"] = os.pathsep.join(search_paths + ([existing_path] if existing_path else []))
    return env


def benjamini_hochberg(pvalues: pd.Series) -> pd.Series:
    """Benjamini-Hochberg FDR correction preserving the original index."""

    if pvalues.empty:
        return pvalues.copy()

    clipped = pvalues.astype(float).fillna(1.0).clip(lower=0.0, upper=1.0)
    order = np.argsort(clipped.to_numpy())
    ranked = clipped.to_numpy()[order]
    n = len(ranked)
    adjusted = np.empty(n, dtype=float)
    prev = 1.0
    for idx in range(n - 1, -1, -1):
        rank = idx + 1
        value = ranked[idx] * n / rank
        prev = min(prev, value)
        adjusted[idx] = prev
    result = np.empty(n, dtype=float)
    result[order] = np.clip(adjusted, 0.0, 1.0)
    return pd.Series(result, index=pvalues.index)


def motif_family_name(motif_name: str) -> str:
    return motif_name.split("/")[0].strip()


def motif_name_from_column(column_name: str) -> str:
    if column_name.endswith(MOTIF_COLUMN_SUFFIX):
        return column_name[: -len(MOTIF_COLUMN_SUFFIX)]
    return column_name


def build_ranked_peak_table(
    diff_res: pd.DataFrame,
    *,
    score_metric: str = "signed_product",
) -> pd.DataFrame:
    """Create a signed, preranked peak table suitable for motif scoring."""

    required = {"log2FC", "pvalue"}
    missing = required - set(diff_res.columns)
    if missing:
        raise ValueError(
            "Differential results are missing required columns for motif ranking: "
            + ", ".join(sorted(missing))
        )

    peak_labels = (
        diff_res["Peak"].astype(str)
        if "Peak" in diff_res.columns
        else pd.Series(diff_res.index.astype(str), index=diff_res.index)
    )
    ranked = diff_res.copy().reset_index(drop=True)
    ranked["Peak"] = peak_labels.to_numpy()
    effect_col = "log2FC_shrunk" if "log2FC_shrunk" in ranked.columns else "log2FC"
    effect = ranked[effect_col].astype(float).fillna(0.0)

    pvalues = ranked["pvalue"].astype(float).replace([np.inf, -np.inf], np.nan).fillna(1.0)
    positive = pvalues[pvalues > 0]
    floor = float(positive.min() / 10.0) if not positive.empty else 1e-300
    floor = max(floor, 1e-300)
    clipped_p = pvalues.clip(lower=floor, upper=1.0)
    neglog10_p = -np.log10(clipped_p)

    if score_metric == "signed_product":
        rank_score = effect * neglog10_p
    elif score_metric == "signed_log10p":
        rank_score = np.sign(effect) * neglog10_p
    elif score_metric == "signed_lfc":
        rank_score = effect
    else:
        raise ValueError(
            f"Unsupported motif score metric '{score_metric}'. "
            "Expected one of: signed_product, signed_log10p, signed_lfc."
        )

    ranked["motif_rank_score"] = rank_score
    ranked["motif_rank_weight"] = np.abs(rank_score)
    ranked = ranked.sort_values(
        ["motif_rank_score", effect_col, "padj", "pvalue", "Peak"],
        ascending=[False, False, True, True, True],
        na_position="last",
    ).reset_index(drop=True)
    ranked["motif_rank"] = np.arange(1, len(ranked) + 1)
    ranked["motif_rank_percentile"] = 1.0 - (ranked["motif_rank"] - 1) / max(1, len(ranked) - 1)
    return ranked


def _running_enrichment_score(
    scores_sorted: np.ndarray,
    hit_mask_sorted: np.ndarray,
    *,
    weight: float = 1.0,
) -> float:
    hit_total = int(hit_mask_sorted.sum())
    miss_total = int((~hit_mask_sorted).sum())
    if hit_total == 0 or miss_total == 0:
        return float("nan")

    hit_weights = np.abs(scores_sorted[hit_mask_sorted]) ** weight
    if not np.isfinite(hit_weights).all() or float(hit_weights.sum()) <= 0:
        hit_weights = np.ones(hit_total, dtype=float)

    running = np.full(scores_sorted.shape[0], -1.0 / miss_total, dtype=float)
    running[hit_mask_sorted] = hit_weights / float(hit_weights.sum())
    running = np.cumsum(running)

    pos = float(np.nanmax(running))
    neg = float(np.nanmin(running))
    return pos if abs(pos) >= abs(neg) else neg


def _top_supporting_peaks(
    peak_ids_sorted: np.ndarray,
    hit_mask_sorted: np.ndarray,
    *,
    positive_tail: bool,
    top_peaks: int,
) -> str:
    hit_ids = peak_ids_sorted[hit_mask_sorted]
    if hit_ids.size == 0 or top_peaks <= 0:
        return ""
    if not positive_tail:
        hit_ids = hit_ids[::-1]
    return ",".join(hit_ids[:top_peaks].tolist())


def _score_membership(
    name: str,
    membership: np.ndarray,
    ranked_peaks: pd.DataFrame,
    *,
    positive_condition: str,
    negative_condition: str,
    gsea_weight: float,
    top_peaks: int,
    member_motifs: Optional[Sequence[str]] = None,
) -> Dict[str, object]:
    scores = ranked_peaks["motif_rank_score"].to_numpy(dtype=float)
    peak_ids = ranked_peaks["Peak"].to_numpy(dtype=str)
    n_total = int(scores.size)
    n_hits = int(membership.sum())
    n_misses = n_total - n_hits
    hit_scores = scores[membership]
    miss_scores = scores[~membership]

    if n_hits == 0 or n_misses == 0:
        raise ValueError("Motif scoring requires both motif-bearing and background peaks")

    try:
        u_stat, pvalue = stats.mannwhitneyu(
            hit_scores,
            miss_scores,
            alternative="two-sided",
            method="asymptotic",
        )
    except ValueError:
        u_stat, pvalue = float("nan"), 1.0

    auc = float(u_stat) / float(n_hits * n_misses) if np.isfinite(u_stat) else float("nan")
    delta_auc = 2.0 * (auc - 0.5) if np.isfinite(auc) else float("nan")
    direction = (
        positive_condition
        if delta_auc > 0
        else negative_condition
        if delta_auc < 0
        else "balanced"
    )
    es = _running_enrichment_score(scores, membership, weight=gsea_weight)
    support = _top_supporting_peaks(
        peak_ids,
        membership,
        positive_tail=not np.isfinite(delta_auc) or delta_auc >= 0,
        top_peaks=top_peaks,
    )
    hit_ranks = ranked_peaks.loc[membership, "motif_rank"].to_numpy(dtype=float)

    row: Dict[str, object] = {
        "motif_name": name,
        "n_motif_peaks": n_hits,
        "motif_peak_fraction": n_hits / float(n_total),
        "mean_peak_score": float(np.mean(hit_scores)),
        "median_peak_score": float(np.median(hit_scores)),
        "mean_rank": float(np.mean(hit_ranks)),
        "median_rank": float(np.median(hit_ranks)),
        "rank_auc": auc,
        "delta_auc": delta_auc,
        "pvalue": float(pvalue) if np.isfinite(pvalue) else 1.0,
        "gsea_es": es,
        "direction": direction,
        "positive_condition": positive_condition,
        "negative_condition": negative_condition,
        "top_supporting_peaks": support,
    }
    if member_motifs is not None:
        row["member_motif_count"] = len(member_motifs)
        row["member_motifs"] = ";".join(member_motifs)
    return row


def score_motif_presence(
    ranked_peaks: pd.DataFrame,
    motif_presence: pd.DataFrame,
    *,
    positive_condition: str,
    negative_condition: str,
    gsea_weight: float = 1.0,
    min_peaks: int = 10,
    max_fraction: float = 0.8,
    top_peaks: int = 25,
) -> pd.DataFrame:
    """Score exact motifs using the ranked peak list."""

    if "Peak" not in ranked_peaks.columns:
        raise ValueError("ranked_peaks must include a Peak column")

    motif_presence = motif_presence.reindex(ranked_peaks["Peak"], fill_value=False)
    rows = []
    total = len(ranked_peaks)
    max_hits = max(1, int(total * max_fraction))
    for motif_name in motif_presence.columns:
        membership = motif_presence[motif_name].to_numpy(dtype=bool)
        n_hits = int(membership.sum())
        if n_hits < min_peaks or n_hits > max_hits:
            continue
        row = _score_membership(
            motif_name,
            membership,
            ranked_peaks,
            positive_condition=positive_condition,
            negative_condition=negative_condition,
            gsea_weight=gsea_weight,
            top_peaks=top_peaks,
        )
        row["motif_family"] = motif_family_name(motif_name)
        rows.append(row)

    result = pd.DataFrame(rows)
    if result.empty:
        return result

    result["fdr"] = benjamini_hochberg(result["pvalue"])
    result.sort_values(
        ["fdr", "pvalue", "delta_auc", "gsea_es", "motif_name"],
        ascending=[True, True, False, False, True],
        inplace=True,
    )
    result.reset_index(drop=True, inplace=True)
    return result


def score_motif_families(
    ranked_peaks: pd.DataFrame,
    motif_presence: pd.DataFrame,
    *,
    positive_condition: str,
    negative_condition: str,
    gsea_weight: float = 1.0,
    min_peaks: int = 10,
    max_fraction: float = 0.8,
    top_peaks: int = 25,
) -> pd.DataFrame:
    """Aggregate exact motifs to motif families via peak-set union and rescore."""

    family_map: Dict[str, list[str]] = {}
    for motif_name in motif_presence.columns:
        family_map.setdefault(motif_family_name(motif_name), []).append(motif_name)

    motif_presence = motif_presence.reindex(ranked_peaks["Peak"], fill_value=False)
    rows = []
    total = len(ranked_peaks)
    max_hits = max(1, int(total * max_fraction))
    for family_name, members in family_map.items():
        membership = motif_presence[members].any(axis=1).to_numpy(dtype=bool)
        n_hits = int(membership.sum())
        if n_hits < min_peaks or n_hits > max_hits:
            continue
        row = _score_membership(
            family_name,
            membership,
            ranked_peaks,
            positive_condition=positive_condition,
            negative_condition=negative_condition,
            gsea_weight=gsea_weight,
            top_peaks=top_peaks,
            member_motifs=members,
        )
        row["motif_family"] = family_name
        rows.append(row)

    result = pd.DataFrame(rows)
    if result.empty:
        return result

    result["fdr"] = benjamini_hochberg(result["pvalue"])
    result.sort_values(
        ["fdr", "pvalue", "delta_auc", "gsea_es", "motif_name"],
        ascending=[True, True, False, False, True],
        inplace=True,
    )
    result.reset_index(drop=True, inplace=True)
    return result


def run_motif_scan(
    consensus_bed: Path,
    genome_fasta: Path,
    motif_file: Path,
    output_path: Path,
    *,
    threads: int = 1,
) -> Path:
    """Scan consensus peaks with HOMER and save motif counts to a gzipped TSV."""

    annotate_tool = Path(resolve_homer_tool("annotatePeaks.pl")).resolve()
    env = build_homer_runtime_env(annotate_tool)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(annotate_tool),
        str(consensus_bed),
        str(genome_fasta),
        "-size",
        "given",
        "-m",
        str(motif_file),
        "-nmotifs",
        "-noann",
        "-cpu",
        str(max(1, threads)),
    ]
    logging.info("Running motif scan: %s", " ".join(cmd))
    temp_output = output_path.with_suffix("") if output_path.suffix == ".gz" else output_path
    with open(temp_output, "w", encoding="utf-8") as handle:
        result = subprocess.run(
            cmd,
            check=False,
            stdout=handle,
            stderr=subprocess.PIPE,
            text=True,
            env=env,
        )
    if result.returncode != 0:
        raise RuntimeError(
            f"HOMER motif scan failed with exit code {result.returncode}: "
            f"{result.stderr.strip()}"
        )
    with open(temp_output, "r", encoding="utf-8") as handle:
        header = handle.readline()
        first_data_row = handle.readline()
    if not header.strip():
        raise RuntimeError(f"HOMER motif scan produced an empty output file: {temp_output}")
    if not first_data_row:
        raise RuntimeError(
            "HOMER motif scan produced only a header row. "
            "This usually means annotatePeaks.pl could not find helper binaries on PATH. "
            f"stderr: {result.stderr.strip()}"
        )
    if output_path.suffix == ".gz":
        with open(temp_output, "rb") as src, gzip.open(output_path, "wb") as dst:
            shutil.copyfileobj(src, dst)
        temp_output.unlink()
    if result.stderr.strip():
        logging.info("annotatePeaks.pl: %s", result.stderr.strip())
    return output_path


def load_motif_presence(scan_path: Path) -> pd.DataFrame:
    """Load motif counts from HOMER output and convert them to presence/absence."""

    open_fn = gzip.open if scan_path.suffix == ".gz" else open
    compression = "gzip" if scan_path.suffix == ".gz" else None

    with open_fn(scan_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")

    peak_col = next((col for col in header if col.startswith("PeakID")), None)
    if peak_col is None:
        raise ValueError(f"Motif scan output is missing a PeakID column: {scan_path}")
    motif_cols = [col for col in header if col.endswith(MOTIF_COLUMN_SUFFIX)]
    if not motif_cols:
        raise ValueError(f"Motif scan output does not contain motif columns: {scan_path}")

    usecols = [peak_col] + motif_cols
    scan_df = pd.read_csv(
        scan_path,
        sep="\t",
        compression=compression,
        usecols=usecols,
        low_memory=False,
    )
    scan_df.rename(columns={peak_col: "Peak"}, inplace=True)
    scan_df["Peak"] = scan_df["Peak"].astype(str)
    rename_map = {col: motif_name_from_column(col) for col in motif_cols}
    scan_df.rename(columns=rename_map, inplace=True)
    motif_df = scan_df.set_index("Peak")
    for col in motif_df.columns:
        motif_df[col] = pd.to_numeric(motif_df[col], errors="coerce").fillna(0).astype(np.int16)
    return motif_df.gt(0)


def run_pairwise_motif_ranking(
    *,
    consensus_bed: Path,
    diff_res: pd.DataFrame,
    output_dir: Path,
    genome_fasta: Path,
    motif_file: Optional[Path] = None,
    positive_condition: str,
    negative_condition: str,
    score_metric: str = "signed_product",
    gsea_weight: float = 1.0,
    min_peaks: int = 10,
    max_fraction: float = 0.8,
    top_peaks: int = 25,
    threads: int = 1,
) -> Dict[str, Path]:
    """Run motif scan + pair-wise ranked enrichment on a single contrast."""

    motif_file = Path(motif_file) if motif_file is not None else resolve_default_known_motif_file()
    output_dir.mkdir(parents=True, exist_ok=True)

    ranked_peaks = build_ranked_peak_table(diff_res, score_metric=score_metric)
    peak_ranking_path = output_dir / "peak_ranking.tsv"
    ranked_peaks.to_csv(peak_ranking_path, sep="\t", index=False)

    scan_path = output_dir / "motif_scan_counts.tsv.gz"
    run_motif_scan(
        consensus_bed,
        genome_fasta,
        motif_file,
        scan_path,
        threads=threads,
    )
    motif_presence = load_motif_presence(scan_path)
    exact_scores = score_motif_presence(
        ranked_peaks,
        motif_presence,
        positive_condition=positive_condition,
        negative_condition=negative_condition,
        gsea_weight=gsea_weight,
        min_peaks=min_peaks,
        max_fraction=max_fraction,
        top_peaks=top_peaks,
    )
    family_scores = score_motif_families(
        ranked_peaks,
        motif_presence,
        positive_condition=positive_condition,
        negative_condition=negative_condition,
        gsea_weight=gsea_weight,
        min_peaks=min_peaks,
        max_fraction=max_fraction,
        top_peaks=top_peaks,
    )

    exact_path = output_dir / "motif_pairwise_scores.tsv"
    family_path = output_dir / "motif_family_pairwise_scores.tsv"
    exact_scores.to_csv(exact_path, sep="\t", index=False)
    family_scores.to_csv(family_path, sep="\t", index=False)

    summary_path = output_dir / "motif_pairwise_summary.tsv"
    summary_rows = []
    for label, frame in (("exact", exact_scores), ("family", family_scores)):
        if frame.empty:
            continue
        top = frame.iloc[0]
        summary_rows.append(
            {
                "level": label,
                "motif_name": top["motif_name"],
                "direction": top["direction"],
                "fdr": top["fdr"],
                "delta_auc": top["delta_auc"],
                "gsea_es": top["gsea_es"],
                "n_motif_peaks": top["n_motif_peaks"],
            }
        )
    pd.DataFrame(summary_rows).to_csv(summary_path, sep="\t", index=False)

    return {
        "peak_ranking": peak_ranking_path,
        "scan_counts": scan_path,
        "exact_scores": exact_path,
        "family_scores": family_path,
        "summary": summary_path,
    }
