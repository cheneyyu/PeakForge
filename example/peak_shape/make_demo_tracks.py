#!/usr/bin/env python3
"""Generate synthetic bigWig tracks for the peak_shape example.

The script creates two small bigWig files (`demo_sample_A.bw` and
`demo_sample_B.bw`) alongside a BED file of candidate regions.  The generated
tracks live on a toy chromosome ``chrDemo`` (length 20 kb) so the example can be
executed entirely offline.

Running this script repeatedly will overwrite the existing demo artefacts.
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Tuple

import numpy as np
import pyBigWig


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default=Path(__file__).resolve().parent / "data",
        type=Path,
        help="Directory where the demo bigWigs and BED file will be written",
    )
    return parser.parse_args()


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def synthetic_profile(center: int, width: float, length: int) -> np.ndarray:
    """Return a Gaussian bump centred at ``center`` with the given ``width``."""

    positions = np.arange(length)
    return np.exp(-0.5 * ((positions - center) / width) ** 2)


def write_bigwig(path: Path, chrom: str, values: np.ndarray) -> None:
    ensure_dir(path.parent)
    bw = pyBigWig.open(str(path), "w")
    try:
        bw.addHeader([(chrom, len(values))])
        starts = np.arange(len(values))
        ends = starts + 1
        bw.addEntries([chrom] * len(values), starts.tolist(), ends=ends.tolist(), values=values.tolist())
    finally:
        bw.close()


def create_demo_tracks(output: Path) -> Tuple[Path, Path]:
    chrom = "chrDemo"
    length = 20_000

    base_signal = synthetic_profile(7_500, 900, length)
    secondary = synthetic_profile(12_500, 600, length)

    sample_a = 2.0 * base_signal + 0.4 * secondary
    sample_b = 1.2 * base_signal + 1.8 * synthetic_profile(13_200, 750, length)

    # Introduce subtle asymmetry to make centroid/skewness informative.
    slope = np.linspace(0.9, 1.1, length)
    sample_a *= slope
    sample_b *= slope[::-1]

    # Add small baseline noise so flank windows remain non-zero after normalisation.
    rng = np.random.default_rng(seed=13)
    sample_a += rng.normal(0, 0.02, length)
    sample_b += rng.normal(0, 0.02, length)

    # Floor negative values at zero.
    sample_a = np.clip(sample_a, 0, None)
    sample_b = np.clip(sample_b, 0, None)

    bw_a = output / "demo_sample_A.bw"
    bw_b = output / "demo_sample_B.bw"

    write_bigwig(bw_a, chrom, sample_a)
    write_bigwig(bw_b, chrom, sample_b)

    bed_path = output / "demo_regions.bed"
    with bed_path.open("w") as bed:
        bed.write("chrDemo\t7000\t8200\tcore_peak\n")
        bed.write("chrDemo\t11500\t13000\tsecondary_peak\n")
        bed.write("chrDemo\t15000\t17000\tflank_associated\n")

    return bw_a, bw_b


def main() -> None:
    args = parse_args()
    ensure_dir(args.output)
    bw_a, bw_b = create_demo_tracks(args.output)
    bed_path = args.output / "demo_regions.bed"

    print("Demo bigWig files written to:")
    print(f"  {bw_a}")
    print(f"  {bw_b}")
    print(f"Regions BED: {bed_path}")


if __name__ == "__main__":
    main()
