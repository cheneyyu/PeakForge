"""Shared IO utilities for PeakForge."""
from __future__ import annotations

from pathlib import Path
from typing import Mapping, MutableMapping, Sequence

import pandas as pd

BED_COLUMNS: tuple[str, str, str] = ("Chromosome", "Start", "End")


def _resolve_path(path: Path | str) -> Path:
    if isinstance(path, Path):
        return path
    return Path(path)


def read_bed_frame(
    path: Path | str,
    *,
    column_names: Sequence[str] = BED_COLUMNS,
    min_columns: int | None = None,
    comment: str = "#",
    dtype: Mapping[int, object] | None = None,
) -> pd.DataFrame:
    """Load a BED-like table into a :class:`pandas.DataFrame`.

    Parameters
    ----------
    path:
        Location of the file to read.
    column_names:
        Names to assign to the first ``len(column_names)`` columns.
    min_columns:
        Minimum number of columns that must be present. Defaults to the
        number of ``column_names``.
    comment:
        Comment indicator passed to :func:`pandas.read_csv`.
    dtype:
        Optional dtype overrides for individual columns.
    """

    target = _resolve_path(path)
    required = len(column_names) if min_columns is None else min_columns
    overrides: MutableMapping[int, object] = {0: str}
    if dtype:
        overrides.update(dtype)

    try:
        frame = pd.read_csv(
            target,
            sep="\t",
            comment=comment,
            header=None,
            dtype=overrides,
        )
    except Exception as exc:  # pragma: no cover - surface informative error
        raise RuntimeError(f"Failed to read BED-like file {target}: {exc}") from exc

    if frame.shape[1] < required:
        raise ValueError(
            f"BED-like file {target} must have at least {required} columns;"
            f" found {frame.shape[1]}"
        )

    base = frame.iloc[:, : len(column_names)].copy()
    base.columns = list(column_names)
    return base


def ensure_integer_columns(frame: pd.DataFrame, columns: Sequence[str]) -> pd.DataFrame:
    """Return a copy of ``frame`` with specified columns coerced to integers."""

    result = frame.copy()
    for column in columns:
        result[column] = pd.to_numeric(result[column], errors="raise").astype(int)
    return result
