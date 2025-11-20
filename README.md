# PeakForge

PeakForge is a Python-native, DiffBind-style toolkit for end-to-end ATAC-seq, CUT&Tag, and ChIP-seq differential analysis. It ingests BAM files or pre-called MACS2 peaks, builds consensus intervals, quantifies coverage, and performs replicate-aware or no-replicate differential testing. Single-sample contrasts remain supported through an implementation of the MARS (MA-plot-based Random Sampling) test so motif and regulon shifts remain discoverable even when only rich individual libraries are available.

---

## Key capabilities

### Flexible inputs
- Accepts paired-end or single-end BAM files.
- Consumes existing MACS2 peak files (`narrowPeak` or `broadPeak`).
- Automatically launches MACS2 when only BAMs are supplied, with support for narrow or broad peak modes.

### Consensus peak management
- Enforces a configurable minimum overlap between samples.
- Expands narrow peaks symmetrically (default ±250 bp via `--peak-extension`) while leaving broad calls intact.
- Optionally reuses an existing consensus BED to guarantee identical genomic intervals between runs.

### Counting and quantification
- Uses deepTools `multiBamSummary BED-file` to build a counts matrix that is exported as TSV alongside the `.npz` archive.
- Calculates library sizes via `samtools idxstats` for single-sample MARS (MA-plot-based Random Sampling) testing.

### Differential analysis
- Automatically chooses between the PyDESeq2 workflow (replicates present) and the MARS (MA-plot-based Random Sampling) test (no replicates).
- Supports contrasts with multiple replicates per condition as well as 1 vs 1 comparisons.
- Emits `differential_results.tsv` plus optional annotations and enrichment tables.

### Optional prior integration
PeakForge can blend public resources into the analysis to stabilise peak statistics or down-weight outliers that deviate from
well-characterised regulatory regions.  The `PriorRegistry` orchestrates this behaviour and is shared by both the main pipeline
and the `peakforge peakshape` subcommand.

#### Building priors in practice
Many users asked how to *produce* and *ship* priors rather than only consume them.  A pragmatic workflow is:

1. **Curate a BED catalogue** – Start from high-quality replicates or public consortia data (ENCODE, BLUEPRINT, etc.).  Merge
   overlapping intervals (for example with `bedtools merge -d 50`) so the catalogue reflects the regulatory territory you
   trust.  This file becomes the `--prior-bed` input.
2. **Summarise intensity with bigWigs** – Generate signal tracks that reflect the typical coverage over those regions.  A common
   pattern is `bamCoverage --normalizeUsing RPGC ... --outFileName reference.bw` followed by `multiBigwigSummary BED-file --bwfiles reference.bw --outRawCounts reference_prior.tsv --outFileName /tmp/ignore.npz`.  The resulting bigWig is supplied through
   `--prior-bigwig` so PeakForge can compute expected read densities per prior interval.
3. **Capture shape statistics** – Run `peakforge peakshape` (or `python peak_shape.py`) on the same training samples to export
   per-peak metrics.  The command `peakforge peakshape --peaks curated.bed --bam reference.bam --prior-shape stats.tsv` will emit
   a TSV that records distributions for width, summit sharpness, shoulder ratios, etc.  These statistics inform z-score scaling
   and novelty penalties when passed to `--prior-shape` / `--prior-stats`.
4. **Describe everything with a manifest** – Place the BED, bigWig, and stats file in a directory alongside a `prior_manifest.json`:

   ```json
   {
     "prior_bed": "prior_catalogue.bed",
     "prior_bigwig": "reference.bw",
     "prior_stats": "stats.tsv",
     "prior_weight": 0.4
   }
   ```

   Shipping this folder with your project makes it trivial for collaborators to reproduce the same prior-aware analysis via
   `--prior-manifest`.

The `example/run_with_prior.sh` script bootstraps a miniature manifest that demonstrates the whole workflow end to end.

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

#### Why priors matter for 1 vs 1 contrasts
Single replicate contrasts (1 treatment vs 1 control) lean on the MARS (MA-plot-based Random Sampling) test, which has limited power to distinguish noise from signal without additional information.  Priors counteract this by providing reference behaviour for well-characterised loci:

- **Overlap-driven shrinkage** – When an observed peak intersects the prior catalogue, its log fold-change is retained.  Novel
  peaks have their `log2FC` (and, when available, `log2FC_shrunk`) multiplied by `1 - w · (1 - overlap_fraction)`, where `w` is
  `--prior-weight`.  This tempers aggressive fold-changes from noisy loci while preserving validated regions.
- **Penalty-aware p-values** – Novel peaks also inherit a `penalty_factor = 1 + w · (1 - overlap_fraction)` that inflates
  p-values and standard errors.  When the prior weight is non-zero this makes the MARS output less overconfident in unexplored
  regions, which is especially helpful with only two libraries.
- **Shape sanity checks** – Feeding shape priors into `peakforge peakshape` highlights peaks with atypical summit profiles
  (`PriorShapeCategory`), flagging candidates that may deserve manual inspection before being called differential hits.

When tuning for 1v1 work:

1. Start with a moderate `--prior-weight` such as `0.3–0.5`.  Inspect `differential_results_prior.tsv` to ensure familiar loci
   remain high-confidence while novel peaks are merely down-weighted rather than discarded.
2. Review the overlap tables (`peaks_prior_overlap.tsv`, `peaks_prior_novel.tsv`) to gauge how much of your catalogue is
   supported by the data.
3. Use the generated `prior_vs_observed.png` plot to confirm that observed signal distributions roughly align with the priors;
   large discrepancies indicate that the weight should be reduced or the catalogue refreshed.

Following these steps gives the sparse 1 vs 1 branch enough context to stabilise fold-change estimates without hiding truly novel
biology.

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
4. Automatically select PyDESeq2 (replicates present) or the MARS (MA-plot-based Random Sampling) test (no replicates) for differential analysis.
5. Generate plots, summary tables, and optional annotations/enrichment reports under the output directory.

## Installation

### Requirements
- Python ≥ 3.10 (tested with Python 3.10).
- Python libraries: `numpy`, `pandas`, `scipy`, `statsmodels`, `matplotlib`, `seaborn`, `pyranges`, `gseapy`, `pydeseq2`.
- External tools: [MACS2][macs2], [deepTools][deeptools], and [samtools][samtools].

### Install commands
```bash
pip install numpy pandas scipy statsmodels matplotlib seaborn pyranges gseapy pydeseq2 macs2 deeptools
conda install -c bioconda samtools
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

The sheet can be tab- or comma-delimited and must include `sample`, `condition`, and `bam`. Optional columns let you supply pre-called peaks as MACS2 `narrowPeak` or `broadPeak` files.

#### 1 vs 1 

| sample   | condition | bam                         | peaks                                   | peak_type |
|----------|-----------|-----------------------------|-----------------------------------------|-----------|
| K562_rep1 | treated   | example/data/K562_rep1.bam  | example/results/narrow/K562_rep1_peaks.narrowPeak | narrow    |
| HepG2_rep1 | control | example/data/HepG2_rep1.bam | example/results/narrow/HepG2_rep1_peaks.narrowPeak | narrow    |

#### 2 vs 2 

| sample     | condition | bam                         | peaks                                        | peak_type |
|------------|-----------|-----------------------------|----------------------------------------------|-----------|
| K562_rep1  | K562      | example/data/K562_rep1.bam  | example/results/broad/K562_rep1_peaks.broadPeak | broad     |
| K562_rep2  | K562      | example/data/K562_rep2.bam  | example/results/broad/K562_rep2_peaks.broadPeak | broad     |
| HepG2_rep1 | HepG2     | example/data/HepG2_rep1.bam | example/results/broad/HepG2_rep1_peaks.broadPeak | broad     |
| HepG2_rep2 | HepG2     | example/data/HepG2_rep2.bam | example/results/broad/HepG2_rep2_peaks.broadPeak | broad     |

### Example command
```bash
./peakforge tsvmode samples.tsv \
  --output-dir results \
  --peak-dir peaks \
  --peak-type narrow \
  --peak-extension 250 \
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
Generates results under `example/results/`, including `consensus_peaks.bed` that can be reused with `--consensus-peaks` in later runs. Internally the script executes the following PeakForge command (with `THREADS` defaulting to 16):

```bash
./peakforge tsvmode example/data/metadata.tsv \
  --output-dir example/results/2v2 \
  --peak-dir example/results/2v2/peaks \
  --peak-type narrow \
  --peak-extension 250 \
  --min-overlap 2 \
  --macs2-genome hs \
  --threads 16
```

### No-replicate 1 vs 1 analysis
```bash
bash example/run_example_1v1.sh
```
Exercises the MARS branch via `metadata_1v1.tsv`, producing outputs in `example/results_1v1/`. The wrapped CLI call is:

```bash
./peakforge tsvmode example/data/metadata_1v1.tsv \
  --output-dir example/results_1v1 \
  --peak-dir example/results_1v1/peaks \
  --peak-type narrow \
  --peak-extension 250 \
  --min-overlap 1 \
  --macs2-genome hs \
  --threads 16
```

### Custom inputs with existing peaks
```bash
bash example/run_example2.sh \
  --condition-a K562 \
  --a-bams example/data/K562_rep1.bam example/data/K562_rep2.bam \
  --a-peaks example/results/2v2/peaks/K562_rep1_peaks.narrowPeak example/results/2v2/peaks/K562_rep2_peaks.narrowPeak \
  --condition-b HepG2 \
  --b-bams example/data/HepG2_rep1.bam example/data/HepG2_rep2.bam \
  --b-peaks example/results/2v2/peaks/HepG2_rep1_peaks.narrowPeak example/results/2v2/peaks/HepG2_rep2_peaks.narrowPeak \
  --consensus-peaks example/results/2v2/consensus_peaks.bed
```
Skips MACS2 when peaks are provided and reuses consensus intervals for consistent comparisons.
The helper ultimately invokes `python chipdiff.py runmode` with the supplied arguments plus `--peak-type narrow --peak-extension 250 --min-overlap 2 --macs2-genome hs --threads 16`, so you can copy/paste the expanded command if you prefer to run PeakForge directly.

### Prior-aware demo
```bash
bash example/run_with_prior.sh
```
Bootstraps priors, runs `tsvmode` with `--prior-manifest`, and performs peak-shape profiling with `--prior-shape`. The manifest contains `prior_bed`, `prior_stats`, and `prior_weight` entries; the script then calls:

```bash
./peakforge tsvmode example/data/metadata_1v1.tsv \
  --output-dir results/prior_demo \
  --peak-dir results/prior_demo/peaks \
  --prior-manifest results/prior_demo/prior_inputs/prior_manifest.json \
  --prior-weight 0.4

./peakforge peakshape \
  --bigwig-a example/peak_shape/data/demo_sample_A.bw \
  --bigwig-b example/peak_shape/data/demo_sample_B.bw \
  --bed example/peak_shape/data/demo_regions.bed \
  --core 400 \
  --flank 800 2000 \
  --out results/prior_demo/shape \
  --prior-shape results/prior_demo/prior_inputs/demo_shape.json \
  --prior-weight 0.4
```

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
