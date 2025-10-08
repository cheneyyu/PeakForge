# PeakForge

PeakForge is a Python-native, DiffBind-style toolkit for end-to-end ATAC-seq, CUT&Tag, and ChIP-seq differential analysis. It ingests BAM files or pre-called MACS2 peaks, builds consensus intervals, quantifies coverage, and performs replicate-aware or no-replicate differential testing. Single-sample contrasts remain supported through an implementation of the MARS test so motif and regulon shifts remain discoverable even when only rich individual libraries are available.

---

## Key capabilities

### Flexible inputs
- Accepts paired-end or single-end BAM files.
- Consumes existing MACS2 peak files (`summits.bed`, `narrowPeak`, or `broadPeak`).
- Automatically launches MACS2 when only BAMs are supplied, with support for narrow, summit, or broad peak modes.

### Consensus peak management
- Enforces a configurable minimum overlap between samples.
- Expands summit and narrow peaks symmetrically (default ±250 bp) while leaving broad calls intact.
- Optionally reuses an existing consensus BED to guarantee identical genomic intervals between runs.

### Counting and quantification
- Uses deepTools `multiBamSummary BED-file` to build a counts matrix that is exported as TSV alongside the `.npz` archive.
- Calculates library sizes via `samtools idxstats` for single-sample MARS testing.

### Differential analysis
- Automatically chooses between the PyDESeq2 workflow (replicates present) and the MARS test (no replicates).
- Supports contrasts with multiple replicates per condition as well as 1 vs 1 comparisons.
- Emits `differential_results.tsv` plus optional annotations and enrichment tables.

### Optional prior integration
PeakForge can blend public resources into the analysis to stabilise peak statistics or down-weight outliers that deviate from
well-characterised regulatory regions.  The `PriorRegistry` orchestrates this behaviour and is shared by both the main pipeline
and the `peakforge peakshape` subcommand.

#### Supported prior artefacts
- **BED intervals** (`--prior-bed` or manifest `prior_bed`): define a catalogue of curated peaks.  The loader normalises the
  coordinates, converts them into a `PyRanges` object, and computes width distributions alongside per-peak overlap counts.
- **bigWig tracks** (`--prior-bigwig` or manifest `prior_bigwig`): summarise signal intensity across the prior intervals using
  `pyBigWig`.  The median and standard deviation of these summaries are recorded when available.
- **Shape statistics** (`--prior-shape` / `--prior-stats` or manifest `prior_stats`): provide distributional expectations for
  per-peak metrics such as summit sharpness or shoulder ratios.  JSON, TSV, CSV, and wide-format tables are supported.
- **Manifest** (`--prior-manifest`): centralises the above paths and the default mixing `prior_weight`.  Relative paths resolve
  relative to the manifest location, making it easy to ship priors alongside a project.

#### How priors influence scoring
1. Peak widths from the prior BED are used to compute a reference mean and standard deviation.  Each observed peak receives a
   z-score (`WidthZ`).
2. Peaks that overlap at least one prior interval are treated as familiar and inherit the configured `prior_weight` directly.
   Non-overlapping peaks apply a *novelty penalty* that scales with the width z-score, reducing (but not eliminating) the weight
   assigned to novel events.
3. When shape statistics or bigWig intensities are provided, the registry mixes observed scores with the prior expectations,
   yielding adjusted metrics via `adjust_scores`.  This ensures that poorly covered peaks can still be ranked sensibly.
4. The per-peak weights are exposed through `get_consensus_weights`, allowing downstream logic to incorporate them when
   prioritising differential hits or plotting ranked lists.

#### Outputs and reporting
- For every sample PeakForge records overlap tables (`peaks_prior_all.tsv`, `peaks_prior_overlap.tsv`, `peaks_prior_novel.tsv`).
- Consensus peaks inherit the same statistics and are saved with effective weights plus overlap counts.
- Summary JSON files capture aggregate overlap fractions, mean weights, and provenance of prior artefacts for reproducibility.
- Optional distribution JSON (`prior_distributions.json`) and comparative density plots (`prior_vs_observed.png`) can be
  produced for audit trails or supplementary figures.

### Outputs and visualisation
- Volcano plots, MA plots, sample correlation heatmaps, and top-peak heatmaps summarise differential results.
- JSON metadata captures run configuration, library sizes, and provenance of priors for reproducibility.
- The `peakforge peakshape` subcommand emits peak-shape delta metrics and comparison plots for signal profiling.

### Pipeline summary
1. Validate inputs and (optionally) call MACS2 to ensure each sample has peaks.
2. Merge peak calls into a consensus catalogue that respects the requested overlap threshold.
3. Quantify read counts per consensus interval with deepTools and estimate library sizes via `samtools idxstats`.
4. Automatically select PyDESeq2 (replicates present) or the MARS test (no replicates) for differential analysis.
5. Generate plots, summary tables, and optional annotations/enrichment reports under the output directory.

## Installation

### Requirements
- Python ≥ 3.9.
- Python libraries: `numpy`, `pandas`, `scipy`, `statsmodels`, `matplotlib`, `seaborn`, `pyranges`, `gseapy`, `pydeseq2`.
- External tools: [MACS2][macs2], [deepTools][deeptools], and [samtools][samtools].

### Install commands
```bash
pip install numpy pandas scipy statsmodels matplotlib seaborn pyranges gseapy pydeseq2
conda install -c bioconda macs2 deeptools samtools
```

### Verify installation
After installing the dependencies, confirm that the CLI is reachable and that optional tooling is on your `PATH`:

```bash
./peakforge --help
which macs2 samtools multiBamSummary
```

If any of the external tools are missing you can re-run the `conda install` command or add them to an existing environment.

---

## Quick start

1. Prepare a metadata file (`samples.tsv`) describing your libraries.
2. Run the pipeline: `./peakforge tsvmode samples.tsv --output-dir results` (or `python chipdiff.py tsvmode ...`).
3. Review tables, plots, and metadata under `results/`.

### Sample sheet format

The sheet can be tab- or comma-delimited and must include `sample`, `condition`, and `bam`. Optional columns let you supply pre-called peaks.

| sample | condition | bam             | peaks                 | peak_type |
|--------|-----------|-----------------|-----------------------|-----------|
| S1     | treated   | data/S1.bam     | data/S1_summits.bed   | summit    |
| S2     | treated   | data/S2.bam     | -                     | -         |
| C1     | control   | data/C1.bam     | data/C1_peaks.bed     | narrow    |

### Example command
```bash
./peakforge tsvmode samples.tsv \
  --output-dir results \
  --peak-dir peaks \
  --min-overlap 2 \
  --gtf annotations.gtf \
  --enrichr
```
This run performs peak calling (if necessary), constructs consensus intervals, quantifies counts, executes the appropriate differential workflow, annotates peaks, and produces summary plots.

---

## Example workflows

The `example/` directory orchestrates a matched ENCODE MYC ChIP-seq dataset (hg38) containing downsampled BAMs.

### End-to-end 2 vs 2 analysis
```bash
bash example/run_pipeline.sh
```
Generates results under `example/results/`, including `consensus_peaks.bed` that can be reused with `--consensus-peaks` in later runs.

### No-replicate 1 vs 1 analysis
```bash
bash example/run_example_1v1.sh
```
Exercises the MARS branch via `metadata_1v1.tsv`, producing outputs in `example/results_1v1/`.

### Custom inputs with existing peaks
```bash
bash example/run_example2.sh \
  --condition-a K562 \
  --a-bams example/data/K562_rep1.bam example/data/K562_rep2.bam \
  --a-peaks example/results/2v2/peaks/K562_rep1_summits.bed example/results/2v2/peaks/K562_rep2_summits.bed \
  --condition-b HepG2 \
  --b-bams example/data/HepG2_rep1.bam example/data/HepG2_rep2.bam \
  --b-peaks example/results/2v2/peaks/HepG2_rep1_summits.bed example/results/2v2/peaks/HepG2_rep2_summits.bed \
  --consensus-peaks example/results/2v2/consensus_peaks.bed
```
Skips MACS2 when peaks are provided and reuses consensus intervals for consistent comparisons.

### Prior-aware demo
```bash
bash example/run_with_prior.sh
```
Bootstraps priors, runs `tsvmode` with `--prior-manifest`, and performs peak-shape profiling with `--prior-shape`.

---

## Reusing consensus peaks
Any completed PeakForge analysis writes `consensus_peaks.bed` inside the output directory. Passing that file to either `tsvmode` or `runmode` through `--consensus-peaks` preserves genomic intervals across follow-up contrasts, keeping fold-change estimates directly comparable.

---

## Troubleshooting

- **Missing optional dependencies** – PyDESeq2, gseapy, and pyBigWig are only required for specific features. If the CLI warns about a missing module you can either install it (`pip install pydeseq2 gseapy pybigwig`) or run the pipeline without that capability.
- **External tool failures** – PeakForge wraps MACS2, deepTools, and samtools. Check their versions with `macs2 --version`, `multiBamSummary --version`, and `samtools --version` if a subprocess error occurs.
- **Empty consensus sets** – ensure that MACS2 produced peaks or provide existing peak files via the metadata sheet. Lowering `--min-overlap` can help when combining sparse datasets.

---

## Command reference
Run `./peakforge --help` (or `python chipdiff.py --help`) to inspect CLI options including peak-calling parameters, threading controls, annotation, enrichment, and prior configuration.

---

## References
- [PyDESeq2][pydeseq2]
- [MACS2 peak caller][macs2]
- [deepTools suite][deeptools]
- [MARS test (DEGseq)][degseq]
- [gseapy and Enrichr API][gseapy]

[pydeseq2]: https://github.com/owkin/PyDESeq2
[macs2]: https://github.com/macs3-project/MACS
[deeptools]: https://deeptools.readthedocs.io/en/latest/
[degseq]: https://academic.oup.com/bioinformatics/article/26/1/136/199566
[gseapy]: https://gseapy.readthedocs.io/
[samtools]: http://www.htslib.org/
