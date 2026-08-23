# PeakForge

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/cheneyyu/PeakForge/blob/main/colab/PeakForge_Colab_Quickstart.ipynb)

PeakForge is a Python-native toolkit for two-group differential analysis of narrow-peak chromatin profiles. It implements two evidence levels:

- Replicate-supported mode uses PyDESeq2 and is the formal inferential mode. It reports `pvalue`, `padj`, and differential-peak calls.
- Exact `1 vs 1` mode cannot estimate biological variability. It reports exploratory candidate rankings led by library-size-normalized log2 fold change, with a MARS-derived score and explicitly labelled sampling-model diagnostics.

The tested scope of this release is transcription-factor ChIP-seq and replicated ATAC-seq with narrow peaks. CUT&Tag and broad-peak files may be accepted technically, but their performance has not been systematically validated. PeakForge currently supports simple two-group contrasts only; it does not fit paired terms, batches, covariates, interactions, or arbitrary design matrices.

## Installation

PeakForge requires Python 3.10 or newer, `samtools`, deepTools, and a working MACS2/MACS3 installation when peaks are not supplied.

```bash
git clone https://github.com/cheneyyu/PeakForge.git
cd PeakForge
curl -LsSf https://astral.sh/uv/install.sh | sh
uv sync --extra macs3 --extra dev
uv run peakforge --help
samtools --version
macs3 --version
```

Python dependencies are resolved in `pyproject.toml` and pinned in `uv.lock`. Each run also records the detected Python, package, and external-tool versions in `metadata.json`.

## Inputs

PeakForge accepts either a metadata table (`tsvmode`) or direct two-condition arguments (`runmode`). Every sample requires a BAM. A narrowPeak/broadPeak file can be supplied; otherwise PeakForge calls MACS using the selected options.

Matched input/control BAMs are passed only to MACS peak calling with `-c`. They are not subtracted from the differential count matrix and are not used as an offset. Differential counting always uses the IP/assay BAMs.

## Replicate-supported TF ChIP-seq example

```bash
uv run peakforge runmode \
  --condition-a K562 \
  --a-bams K562_rep1.bam K562_rep2.bam K562_rep3.bam \
  --a-controls K562_input.bam \
  --condition-b HepG2 \
  --b-bams HepG2_rep1.bam HepG2_rep2.bam HepG2_rep3.bam \
  --b-controls HepG2_input.bam \
  --macs-format BAM \
  --macs2-qvalue 0.01 \
  --peak-edge-padding 250 \
  --count-unit read \
  --min-mapq 30 \
  --exclude-duplicates \
  --output-dir results/myc_3x3
```

These ENCODE MYC benchmark BAMs are single-end. For paired-end ChIP-seq, use `--macs-format BAMPE --count-unit fragment` when DNA fragments are the intended count unit.

## Replicate-supported paired-end ATAC-seq example

```bash
uv run peakforge runmode \
  --condition-a condition_A \
  --a-bams A_rep1.bam A_rep2.bam A_rep3.bam \
  --condition-b condition_B \
  --b-bams B_rep1.bam B_rep2.bam B_rep3.bam \
  --macs-format BAMPE \
  --macs2-qvalue 0.01 \
  --peak-edge-padding 0 \
  --count-unit fragment \
  --min-mapq 30 \
  --exclude-duplicates \
  --output-dir results/atac_3x3
```

MACS settings are assay- and preprocessing-dependent. `--macs-format`, `--macs-nomodel`, `--macs-shift`, `--macs-extsize`, `--peak-type`, and `--macs-extra` expose the relevant choices. PeakForge records the fully expanded MACS command. There is no universal ATAC/ChIP preset.

`BAMPE` peak calling uses observed paired fragments; PeakForge therefore rejects
nonzero `--macs-shift` or any `--macs-extsize` with `BAMPE`. Read-based,
cut-site-focused MACS peak calling can instead use `--macs-format BAM
--macs-nomodel --macs-shift -50 --macs-extsize 100`. These MACS options affect
peak calling only and do not implement Tn5 insertion-site quantification.

## Counting and library size

The same filter definition is used for interval counts and the library-size denominator:

- `--count-unit read` counts each retained alignment.
- `--count-unit fragment` requires paired-end data, retains the first mate of each proper pair once, and extends it across the paired fragment for interval overlap counting.
- Unmapped, secondary, supplementary, and QC-failed alignments are always excluded.
- `--min-mapq`, `--exclude-duplicates`, `--proper-pairs-only`, and `--sam-flag-exclude` apply consistently to counting and library size.
- `library_size_metrics.tsv` records raw mapped alignments, filtered alignments, filtered countable units, and the selected denominator.
- Tn5 insertion-site counting is not implemented.

Overlapping mates therefore count twice in read mode and once in fragment mode. A fragment may contribute to more than one genomic interval if it overlaps more than one interval.

PeakForge canonicalizes persisted count matrices by chromosome, start, and end
before statistical analysis. Exploratory ranking uses the stable peak identifier
to break exact score ties. These rules make row order and tied ranks independent
of backend thread scheduling.

## Peak consensus coordinates

`--peak-edge-padding N` adds `N` bp to each original narrow-peak edge before clustering peaks into consensus intervals. It does not create a fixed-width interval and does not use the peak midpoint or narrowPeak summit. Broad-peak boundaries are not padded. The historical option name `--peak-extension` remains as an alias.

Alternatively, `--peak-coordinate-mode summit-fixed --summit-fixed-width W` reads the zero-based summit offset from narrowPeak column 10 and replaces every sample peak with an exact total-width `W`-bp interval centered on that summit before consensus clustering. Near chromosome position zero, the interval is shifted right to retain the requested width. This mode requires valid narrowPeak summit fields and is rejected for broadPeak inputs. Coordinate transformation is not applied when `--consensus-peaks` supplies an already constructed interval set.

```bash
uv run peakforge tsvmode samples.tsv \
  --peak-coordinate-mode summit-fixed \
  --summit-fixed-width 500 \
  --min-overlap 2 \
  --output-dir results/summit_500
```

The 250-bp-per-edge setting is an operational ChIP-seq example, not a universal biological optimum. Native boundaries (`0` bp), configurable edge padding, and summit-centered fixed-width intervals should be evaluated for the assay and biological question at hand.

## Exact single-pair exploratory mode

An exact one-sample-per-condition analysis produces columns including:

```text
count_condition_a
count_condition_b
library_size_a
library_size_b
cpm_condition_a
cpm_condition_b
mean_cpm
normalized_log2fc
mars_score
rank_up
rank_down
rank_absolute_lfc
rank_absolute_mars
sampling_pvalue
sampling_qvalue
rank_eligible
zero_count_status
analysis_mode
interpretation
```

`normalized_log2fc` is the primary ranking quantity after the configured mean-CPM filter. `mars_score` is secondary. `sampling_pvalue` and `sampling_qvalue` describe only the conditional sampling model; they do not estimate biological variability, do not provide formal significance, and must not be used to select discoveries.

For backward-compatible Python API use, `mars_differential()` temporarily exposes legacy aliases with a deprecation warning. New CLI outputs do not use generic `pvalue` or `padj` names in single-pair mode.

## Reusing a fixed consensus

```bash
uv run peakforge runmode \
  --condition-a K562 \
  --a-bams K562_rep1.bam K562_rep2.bam K562_rep3.bam \
  --condition-b HepG2 \
  --b-bams HepG2_rep1.bam HepG2_rep2.bam HepG2_rep3.bam \
  --consensus-peaks external_consensus.bed \
  --count-unit read \
  --min-mapq 30 \
  --exclude-duplicates \
  --output-dir results/fixed_consensus
```

This isolates interval construction from downstream counting/statistical testing. Agreement with DiffBind or csaw on real data is interpreted as concordance, not accuracy; known-truth accuracy is evaluated by simulation.

## Output provenance

Every run writes `metadata.json` containing:

- PeakForge version, Git commit, and dirty-state flag;
- Python/package and external-tool versions;
- full invocation and expanded counting/MACS commands;
- sample/control/peak paths;
- count filters and per-sample library-size metrics;
- consensus source and filtering metadata;
- requested peak-coordinate mode, edge padding, summit width, and whether the transformation was applied;
- explicit analysis mode and interpretation.

For reproducible production runs, use a clean version-controlled commit; `metadata.json` records the commit and whether the working tree was dirty.

## Limitations

- Replicate-supported inference requires at least two biological replicates in both groups.
- Exploratory mode requires exactly one sample in each group; mixed `1 vs N` designs are rejected.
- Pairing, batches, covariates, interactions, and general design matrices are not supported.
- Tn5 insertion-site quantification is not supported.
- CUT&Tag and broad-peak performance are not validated in the current release.
- Single-pair sampling-model p/q values are non-inferential diagnostics.

## Development tests

```bash
uv run pytest -q
```

The test suite includes hand-calculated formula cases, zero-count behavior, read/fragment filtering on a synthetic BAM, edge padding, summit-centered exact-width coordinates and invalid-summit handling, consensus support, output schema, thread invariance, and an end-to-end tiny-data smoke test.

## License

PeakForge is distributed under the license in [LICENSE](LICENSE).
