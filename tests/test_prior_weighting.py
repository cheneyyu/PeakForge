import numpy as np
import pandas as pd

from chipdiff import compute_prior_weights, weighted_bh, load_shape_covariates
from prior_utils import PriorRegistry


def test_weighted_bh_shape_and_range():
    p = np.array([0.001, 0.02, 0.2, 0.8], dtype=float)
    w = np.array([1.2, 1.0, 0.8, 1.0], dtype=float)
    p_w, q_w, sig = weighted_bh(p, w, alpha=0.05)
    assert len(p_w) == len(p)
    assert len(q_w) == len(p)
    assert len(sig) == len(p)
    assert np.all((p_w >= 0) & (p_w <= 1))
    assert np.all((q_w >= 0) & (q_w <= 1))


def test_no_prior_like_behavior_when_covariates_missing():
    df = pd.DataFrame({"pvalue": [0.01, 0.2, 0.8]})
    w = compute_prior_weights(df)
    assert np.allclose(w, np.ones(len(df)))


def test_prior_bed_covariates_only():
    peaks = pd.DataFrame(
        {
            "Chromosome": ["chr1", "chr1"],
            "Start": [100, 400],
            "End": [200, 500],
        }
    )
    reg = PriorRegistry(weight=0.3)
    reg.distributions["width_mean"] = 100.0
    reg.distributions["width_std"] = 10.0
    out = reg.match_prior(peaks)
    assert {"PriorWeight", "NoveltyPenalty", "WidthZ", "PriorOverlap"}.issubset(out.columns)


def test_prior_bed_plus_bigwig_covariates_placeholder_without_file():
    consensus = pd.DataFrame(
        {
            "Chromosome": ["chr1"],
            "Start": [100],
            "End": [200],
            "Name": ["peak_1"],
        }
    )
    reg = PriorRegistry(prior_bigwig="/non/existing/file.bw", weight=0.3)
    out = reg.compute_consensus_intensity_z(consensus)
    assert set(["Name", "PeakId", "PriorIntensity", "IntZ"]).issubset(out.columns)
    assert len(out) == 1


def test_prior_bed_plus_shape_tsv_merge_tolerates_unmatched_rows(tmp_path):
    tsv = tmp_path / "shape.tsv"
    pd.DataFrame(
        {
            "peak_id": ["chr1:10-20"],
            "ShapeZ": [1.5],
        }
    ).to_csv(tsv, sep="\t", index=False)

    peak_index = pd.Index(["chr1:999-1999"])  # unmatched on purpose
    consensus_df = pd.DataFrame(
        {
            "Name": ["peak_unmatched"],
            "Chromosome": ["chr1"],
            "Start": [999],
            "End": [1999],
        },
        index=peak_index,
    )
    merged = load_shape_covariates(tsv, peak_index, consensus_df)
    assert len(merged) == 1
    assert "ShapeZ" in merged.columns
    assert merged["ShapeZ"].isna().all()
