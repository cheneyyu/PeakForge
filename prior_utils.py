from __future__ import annotations

import json
import logging
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import numpy as np
import pandas as pd
import pyranges as pr
from io_utils import ensure_integer_columns, read_bed_frame

try:  # pragma: no cover - optional dependency
    import pyBigWig  # type: ignore
except ImportError:  # pragma: no cover - optional dependency
    pyBigWig = None  # type: ignore

LOGGER = logging.getLogger(__name__)


def _ensure_directory(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def load_prior_manifest(path: Path) -> Dict[str, object]:
    """Parse a manifest describing prior resources."""

    if not path.exists():
        raise FileNotFoundError(f"Prior manifest not found: {path}")

    text = path.read_text().strip()
    if not text:
        return {}

    lowered = path.suffix.lower()
    if lowered in {".json", ".js", ".json5"} or text.lstrip().startswith("{"):
        try:
            data = json.loads(text)
        except json.JSONDecodeError as exc:
            raise ValueError(f"Failed to parse JSON prior manifest {path}: {exc}") from exc
        if not isinstance(data, dict):
            raise ValueError(f"Prior manifest {path} must contain a JSON object")
        return data

    result: Dict[str, object] = {}
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if "=" in stripped:
            key, value = stripped.split("=", 1)
        elif "\t" in stripped:
            key, value = stripped.split("\t", 1)
        else:
            parts = stripped.split(None, 1)
            if len(parts) != 2:
                raise ValueError(
                    f"Cannot parse manifest line '{line}' in {path}; expected 'key value'"
                )
            key, value = parts
        result[key.strip()] = value.strip()
    return result


@dataclass
class _ShapeStat:
    mean: Optional[float]
    std: Optional[float]
    median: Optional[float]


class PriorRegistry:
    """Utility class that coordinates the use of public priors."""

    def __init__(
        self,
        *,
        prior_bed: Optional[Path | str] = None,
        prior_bigwig: Optional[Path | str] = None,
        prior_stats: Optional[Path | str] = None,
        weight: float = 0.3,
    ) -> None:
        self.prior_bed_path = Path(prior_bed) if prior_bed else None
        self.prior_bigwig_path = Path(prior_bigwig) if prior_bigwig else None
        self.prior_stats_path = Path(prior_stats) if prior_stats else None
        self.weight = float(weight)

        self.prior_ranges: Optional[pr.PyRanges] = None
        self.distributions: Dict[str, object] = {}
        self.shape_stats: Dict[str, _ShapeStat] = {}
        self.sample_tables: List[pd.DataFrame] = []
        self.consensus_table: Optional[pd.DataFrame] = None
        self.tables_paths: Dict[str, Path] = {}
        self.distribution_path: Optional[Path] = None
        self.prior_plot_path: Optional[Path] = None
        self._summary_cache: Optional[Dict[str, object]] = None

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------

    @property
    def enabled(self) -> bool:
        return any([self.prior_bed_path, self.prior_bigwig_path, self.prior_stats_path])

    @property
    def has_shape_priors(self) -> bool:
        return bool(self.shape_stats)

    # ------------------------------------------------------------------
    # Loading routines
    # ------------------------------------------------------------------

    def load_prior_distributions(self) -> None:
        if not self.enabled:
            return

        if self.prior_bed_path:
            self.prior_ranges = self._load_prior_bed(self.prior_bed_path)
            widths = (self.prior_ranges.df["End"] - self.prior_ranges.df["Start"]).astype(float)
            self.distributions["width_distribution"] = widths
            self.distributions["width_mean"] = float(widths.mean()) if not widths.empty else math.nan
            self.distributions["width_std"] = float(widths.std(ddof=0)) if widths.size > 1 else math.nan

        if self.prior_stats_path:
            stats = self._load_prior_stats(self.prior_stats_path)
            self.shape_stats = stats
            self.distributions["shape_stats"] = {
                key: {
                    "mean": value.mean,
                    "std": value.std,
                    "median": value.median,
                }
                for key, value in stats.items()
            }

        if self.prior_bigwig_path and self.prior_ranges is not None:
            intensities = self._summarise_bigwig(self.prior_bigwig_path, self.prior_ranges)
            if intensities is not None and not intensities.empty:
                self.distributions["intensity_distribution"] = intensities
                self.distributions["intensity_mean"] = float(intensities.mean())
                self.distributions["intensity_std"] = float(intensities.std(ddof=0)) if intensities.size > 1 else math.nan

    # ------------------------------------------------------------------
    # Prior helpers
    # ------------------------------------------------------------------

    def _load_prior_bed(self, path: Path) -> pr.PyRanges:
        frame = read_bed_frame(path)
        frame = ensure_integer_columns(frame, ("Start", "End"))
        return pr.PyRanges(frame)

    def _load_prior_stats(self, path: Path) -> Dict[str, _ShapeStat]:
        lowered = path.suffix.lower()
        data: Dict[str, _ShapeStat] = {}
        if lowered in {".json", ".js", ".json5"}:
            try:
                payload = json.loads(path.read_text())
            except json.JSONDecodeError as exc:
                raise ValueError(f"Failed to parse JSON prior stats {path}: {exc}") from exc
            if not isinstance(payload, dict):
                raise ValueError(f"Prior stats file {path} must contain a JSON object")
            for key, value in payload.items():
                if isinstance(value, dict):
                    data[key] = _ShapeStat(
                        mean=self._coerce_optional_float(value.get("mean")),
                        std=self._coerce_optional_float(value.get("std")),
                        median=self._coerce_optional_float(value.get("median")),
                    )
                elif isinstance(value, (list, tuple)):
                    series = pd.Series(value, dtype=float)
                    data[key] = _ShapeStat(
                        mean=float(series.mean()),
                        std=float(series.std(ddof=0)) if series.size > 1 else math.nan,
                        median=float(series.median()),
                    )
                else:
                    data[key] = _ShapeStat(mean=self._coerce_optional_float(value), std=None, median=None)
            return data

        # Delimited text fallback
        sep = "\t" if lowered.endswith("tsv") else ","
        df = pd.read_csv(path, sep=sep)
        if "metric" in df.columns:
            for _, row in df.iterrows():
                metric = str(row["metric"]).strip()
                mean = self._coerce_optional_float(row.get("mean"))
                std = self._coerce_optional_float(row.get("std"))
                median = self._coerce_optional_float(row.get("median"))
                data[metric] = _ShapeStat(mean=mean, std=std, median=median)
            return data

        # Wide-format numeric columns
        for column in df.columns:
            series = pd.to_numeric(df[column], errors="coerce")
            series = series.dropna()
            if series.empty:
                continue
            data[column] = _ShapeStat(
                mean=float(series.mean()),
                std=float(series.std(ddof=0)) if series.size > 1 else math.nan,
                median=float(series.median()),
            )
        return data

    def _summarise_bigwig(self, bigwig_path: Path, ranges: pr.PyRanges) -> Optional[pd.Series]:
        if pyBigWig is None:
            LOGGER.warning(
                "pyBigWig not installed; cannot use prior bigWig %s for intensity summaries", bigwig_path
            )
            return None
        try:
            handle = pyBigWig.open(str(bigwig_path))
        except Exception as exc:  # pragma: no cover - IO heavy
            LOGGER.warning("Failed to open prior bigWig %s: %s", bigwig_path, exc)
            return None
        try:
            chroms = handle.chroms()
            if not chroms:
                LOGGER.warning("Prior bigWig %s does not expose chromosome sizes", bigwig_path)
                return None
            df = ranges.df
            values: List[float] = []
            limit = min(len(df), 5000)
            for idx in range(limit):
                row = df.iloc[idx]
                chrom = str(row["Chromosome"])
                if chrom not in chroms:
                    continue
                start = int(row["Start"])
                end = int(row["End"])
                try:
                    stat = handle.stats(chrom, start, end, nBins=1, exact=True)[0]
                except Exception:
                    stat = None
                if stat is None or not isinstance(stat, (float, int)):
                    continue
                values.append(float(stat))
            if not values:
                return None
            return pd.Series(values, dtype=float)
        finally:
            try:
                handle.close()
            except Exception:  # pragma: no cover - best effort
                pass

    # ------------------------------------------------------------------
    # Matching
    # ------------------------------------------------------------------

    def match_prior(self, peaks_df: pd.DataFrame) -> pd.DataFrame:
        if peaks_df.empty:
            return pd.DataFrame()

        base = peaks_df[["Chromosome", "Start", "End"]].copy()
        base = base.reset_index(drop=True)
        base["Start"] = pd.to_numeric(base["Start"], errors="coerce")
        base["End"] = pd.to_numeric(base["End"], errors="coerce")
        widths = base["End"] - base["Start"]

        width_mean = self.distributions.get("width_mean")
        width_std = self.distributions.get("width_std")
        if width_std in {0, None}:
            width_std = math.nan

        if width_mean is None:
            width_mean = float(widths.mean()) if widths.notna().any() else math.nan

        if math.isnan(width_std) or width_std == 0:
            width_z = pd.Series(np.nan, index=base.index)
        else:
            width_z = (widths - float(width_mean)) / float(width_std)

        overlap = pd.Series(False, index=base.index)
        overlap_count = pd.Series(0, index=base.index, dtype=int)

        if self.prior_ranges is not None and len(self.prior_ranges) > 0:
            augmented = base.copy()
            augmented["__idx"] = np.arange(len(base))
            current = pr.PyRanges(augmented)
            counted = current.count_overlaps(self.prior_ranges)
            counts_df = counted.df.sort_values("__idx").set_index("__idx")
            overlap_count = counts_df["NumberOverlaps"].reindex(base.index).fillna(0).astype(int)
            overlap = overlap_count > 0

        abs_z = np.abs(width_z.fillna(0.0).to_numpy())
        novelty_penalty_values = np.clip(abs_z, 0, 6) * self.weight
        novelty_penalty_values = np.where(overlap.to_numpy(), 0.0, novelty_penalty_values)
        novelty_penalty = pd.Series(novelty_penalty_values, index=base.index)
        effective_weight_values = np.where(
            overlap.to_numpy(),
            self.weight,
            np.clip(self.weight - novelty_penalty_values, 0.0, self.weight),
        )
        effective_weight = pd.Series(effective_weight_values, index=base.index)

        return pd.DataFrame(
            {
                "Width": widths,
                "WidthZ": width_z,
                "PriorOverlap": overlap,
                "OverlapCount": overlap_count,
                "PriorWeight": effective_weight,
                "NoveltyPenalty": novelty_penalty,
            }
        )

    # ------------------------------------------------------------------
    # Recording helpers
    # ------------------------------------------------------------------

    def record_sample(self, sample_name: str, peaks_df: pd.DataFrame) -> None:
        if not self.enabled:
            return
        match = self.match_prior(peaks_df)
        if match.empty:
            return
        base = peaks_df[["Chromosome", "Start", "End"]].reset_index(drop=True).copy()
        base["Sample"] = sample_name
        table = pd.concat([base, match], axis=1)
        self.sample_tables.append(table)
        self._summary_cache = None

    def record_consensus(self, consensus_df: pd.DataFrame) -> None:
        if not self.enabled or consensus_df.empty:
            return
        base = consensus_df[["Chromosome", "Start", "End", "Name"]].reset_index(drop=True).copy()
        base["PeakId"] = (
            base["Chromosome"].astype(str)
            + ":"
            + base["Start"].astype(int).astype(str)
            + "-"
            + base["End"].astype(int).astype(str)
        )
        match = self.match_prior(base)
        if match.empty:
            return
        combined = pd.concat([base, match], axis=1)
        combined["Peak"] = combined["Name"].fillna(combined["PeakId"])
        self.consensus_table = combined
        self._summary_cache = None

    # ------------------------------------------------------------------
    # Adjustments and scoring
    # ------------------------------------------------------------------

    def adjust_scores(self, scores: pd.Series, *, metric: Optional[str] = None) -> pd.Series:
        if not self.enabled:
            return scores
        if metric and metric in self.shape_stats:
            stat = self.shape_stats[metric]
            target = stat.mean if stat.mean is not None else float(scores.mean())
        else:
            target = self.distributions.get("intensity_mean")
            if target is None:
                target = float(scores.median()) if scores.notna().any() else 0.0
        adjusted = scores * (1.0 - self.weight) + float(target) * self.weight
        return adjusted

    def zscore(self, metric: str, values: pd.Series) -> pd.Series:
        if metric not in self.shape_stats:
            return pd.Series(np.nan, index=values.index, dtype=float)
        stat = self.shape_stats[metric]
        if stat.mean is None or stat.std in {None, 0}:
            return pd.Series(np.nan, index=values.index, dtype=float)
        return (values - float(stat.mean)) / float(stat.std)

    def get_consensus_weights(self, index: Iterable[str]) -> pd.Series:
        index = pd.Index(index)
        if not self.enabled or self.consensus_table is None:
            return pd.Series(0.0, index=index)
        mapping: Dict[str, float] = {}
        for _, row in self.consensus_table.iterrows():
            weight = float(row.get("PriorWeight", 0.0) or 0.0)
            name = row.get("Name")
            peak_id = row.get("PeakId")
            if isinstance(name, str) and name:
                mapping[name] = weight
            if isinstance(peak_id, str) and peak_id:
                mapping[peak_id] = weight
        return pd.Series([mapping.get(str(item), 0.0) for item in index], index=index, dtype=float)

    # ------------------------------------------------------------------
    # Persistence and reporting
    # ------------------------------------------------------------------

    def write_peak_tables(self, output_dir: Path) -> Dict[str, Path]:
        if not self.sample_tables:
            return {}
        combined = pd.concat(self.sample_tables, ignore_index=True)
        overlap = combined[combined["PriorOverlap"]]
        novel = combined[~combined["PriorOverlap"]]

        _ensure_directory(output_dir)
        overlap_path = output_dir / "peaks_prior_overlap.tsv"
        novel_path = output_dir / "peaks_prior_novel.tsv"
        combined_path = output_dir / "peaks_prior_all.tsv"

        combined.to_csv(combined_path, sep="\t", index=False)
        overlap.to_csv(overlap_path, sep="\t", index=False)
        novel.to_csv(novel_path, sep="\t", index=False)

        self.tables_paths = {
            "combined": combined_path,
            "overlap": overlap_path,
            "novel": novel_path,
        }
        return self.tables_paths

    def save_distributions(self, path: Path) -> Optional[Path]:
        if not self.enabled:
            return None
        payload: Dict[str, object] = {
            "weight": self.weight,
            "prior_bed": str(self.prior_bed_path) if self.prior_bed_path else None,
            "prior_bigwig": str(self.prior_bigwig_path) if self.prior_bigwig_path else None,
            "prior_stats": str(self.prior_stats_path) if self.prior_stats_path else None,
        }
        width_series = self.distributions.get("width_distribution")
        if isinstance(width_series, pd.Series) and not width_series.empty:
            payload["width_mean"] = float(width_series.mean())
            payload["width_std"] = float(width_series.std(ddof=0)) if width_series.size > 1 else math.nan
            payload["width_quantiles"] = {
                q: float(np.quantile(width_series, q)) for q in (0.1, 0.5, 0.9)
            }
        if "intensity_mean" in self.distributions:
            payload["intensity_mean"] = self.distributions["intensity_mean"]
            payload["intensity_std"] = self.distributions.get("intensity_std")
        if self.shape_stats:
            payload["shape_stats"] = {
                key: {"mean": stat.mean, "std": stat.std, "median": stat.median}
                for key, stat in self.shape_stats.items()
            }

        _ensure_directory(path.parent)
        with path.open("w") as fh:
            json.dump(payload, fh, indent=2)
        self.distribution_path = path
        return path

    def plot_prior_vs_observed(self, observed: Dict[str, pd.Series], output_path: Path) -> Optional[Path]:
        if not self.enabled:
            return None
        metrics = []
        for name, series in observed.items():
            prior_series = self.distributions.get(f"{name}_distribution")
            shape_stat = self.shape_stats.get(name)
            if prior_series is None and shape_stat is None:
                continue
            metrics.append((name, series, prior_series, shape_stat))
        if not metrics:
            return None

        import matplotlib.pyplot as plt  # Local import to avoid hard dependency in headless contexts
        import seaborn as sns

        rows = len(metrics)
        fig, axes = plt.subplots(rows, 1, figsize=(6, 3 * rows), squeeze=False)
        for idx, (name, observed_series, prior_series, shape_stat) in enumerate(metrics):
            ax = axes[idx, 0]
            obs = observed_series.replace([np.inf, -np.inf], np.nan).dropna()
            if obs.empty:
                ax.text(0.5, 0.5, "No observed data", ha="center", va="center")
                ax.set_title(name)
                continue
            sns.kdeplot(obs, ax=ax, fill=True, color="#4c72b0", label="Observed")
            if isinstance(prior_series, pd.Series) and not prior_series.empty:
                sns.kdeplot(prior_series.dropna(), ax=ax, color="#dd8452", label="Prior", linestyle="--")
            elif shape_stat and shape_stat.mean is not None and shape_stat.std not in {None, 0}:
                xs = np.linspace(obs.min(), obs.max(), 200)
                mean = float(shape_stat.mean)
                std = float(shape_stat.std)
                ys = (1 / (std * math.sqrt(2 * math.pi))) * np.exp(-0.5 * ((xs - mean) / std) ** 2)
                ys = ys / ys.max() * obs.size / max(obs.std(ddof=0), 1e-6)
                ax.plot(xs, ys, color="#dd8452", linestyle="--", label="Prior (μ±σ)")
                ax.axvline(mean, color="#dd8452", linestyle=":")
                if std and std > 0:
                    ax.axvline(mean - std, color="#dd8452", linestyle=":", alpha=0.6)
                    ax.axvline(mean + std, color="#dd8452", linestyle=":", alpha=0.6)
            ax.set_title(f"{name} prior vs observed")
            ax.set_xlabel(name)
            ax.legend(loc="best")
        fig.tight_layout()
        _ensure_directory(output_path.parent)
        fig.savefig(output_path, dpi=300)
        plt.close(fig)
        self.prior_plot_path = output_path
        return output_path

    # ------------------------------------------------------------------
    # Summary / metadata
    # ------------------------------------------------------------------

    def summarize(self) -> Dict[str, object]:
        if self._summary_cache is not None:
            return self._summary_cache
        summary: Dict[str, object] = {
            "enabled": self.enabled,
            "weight": self.weight,
        }
        if self.sample_tables:
            combined = pd.concat(self.sample_tables, ignore_index=True)
            total = len(combined)
            overlap = int(combined["PriorOverlap"].sum())
            summary["sample_peaks_total"] = total
            summary["sample_peaks_overlap"] = overlap
            summary["sample_peaks_fraction_overlap"] = overlap / total if total else 0.0
        if self.consensus_table is not None:
            cons = self.consensus_table
            total = len(cons)
            overlap = int(cons["PriorOverlap"].sum())
            summary["consensus_peaks_total"] = total
            summary["consensus_peaks_overlap"] = overlap
            summary["consensus_peaks_fraction_overlap"] = overlap / total if total else 0.0
            summary["consensus_average_weight"] = float(cons["PriorWeight"].mean()) if total else 0.0
        self._summary_cache = summary
        return summary

    def metadata_entry(self) -> Optional[Dict[str, object]]:
        if not self.enabled:
            return None
        summary = self.summarize()
        entry = {
            "weight": self.weight,
            "prior_bed": str(self.prior_bed_path) if self.prior_bed_path else None,
            "prior_bigwig": str(self.prior_bigwig_path) if self.prior_bigwig_path else None,
            "prior_stats": str(self.prior_stats_path) if self.prior_stats_path else None,
            "tables": {key: str(path) for key, path in self.tables_paths.items()},
            "distribution_json": str(self.distribution_path) if self.distribution_path else None,
            "plot": str(self.prior_plot_path) if self.prior_plot_path else None,
        }
        entry.update(summary)
        return entry

    # ------------------------------------------------------------------
    # Utilities
    # ------------------------------------------------------------------

    @staticmethod
    def _coerce_optional_float(value: object) -> Optional[float]:
        if value is None or (isinstance(value, float) and math.isnan(value)):
            return None
        try:
            return float(value)
        except (TypeError, ValueError):
            return None
