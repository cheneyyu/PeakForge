# PeakForge

PeakForge is a Python-native, DiffBind-style toolkit for end-to-end ATAC-seq, CUT&Tag, and ChIP-seq differential
analysis. It ingests BAM files or pre-called MACS2 peaks, builds consensus intervals, quantifies coverage, and performs
replicate-aware or no-replicate differential testing. Single-sample contrasts remain supported through an implementation
of the MARS (MA-plot-based Random Sampling) test so motif and regulon shifts remain discoverable even when only rich
individual libraries are available.

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

### Outputs and visualisation
- Volcano plots, MA plots, sample correlation heatmaps, and top-peak heatmaps summarise differential results.
- JSON metadata captures run configuration, library sizes, and provenance of priors for reproducibility.
- The `peakforge peakshape` subcommand emits peak-shape delta metrics and comparison plots for signal profiling.

### Pipeline summary
1. Validate inputs and (optionally) call MACS2 to ensure each sample has peaks (paired-end BAMs automatically use `-f BAMPE`).
2. Merge peak calls into a consensus catalogue that respects the requested overlap threshold.
3. Quantify read counts per consensus interval with deepTools and estimate library sizes via `samtools idxstats`.
4. Automatically select PyDESeq2 (replicates present) or the MARS (MA-plot-based Random Sampling) test (no replicates) for differential analysis.
5. Generate plots, summary tables, and optional annotations/enrichment reports under the output directory.

## Installation

### Requirements
- Python ≥ 3.10 (tested with Python 3.10).
- Python libraries: `numpy`, `pandas`, `scipy`, `statsmodels`, `matplotlib`, `seaborn`, `pyranges`, `gseapy`, `pydeseq2`.
- External tools: [MACS2/3][macs2], [deepTools][deeptools] (including `multiBamSummary` and `bamCoverage`), and [samtools][samtools].

### Install commands
```bash
pip install numpy pandas scipy statsmodels matplotlib seaborn pyranges gseapy pydeseq2 macs3 deeptools
conda install -c bioconda samtools
```
PeakForge automatically prefers `macs2` when available and falls back to `macs3` if needed.

### Verify installation
After installing the dependencies, confirm that the CLI is reachable and that optional tooling is on your `PATH`:

```bash
./peakforge --help
which macs2 macs3 samtools multiBamSummary
```

If any of the external tools are missing you can re-run the `conda install` command or add them to an existing
environment.

---

## Quick start

1. Prepare a metadata file (`samples.tsv`) describing your libraries.
2. Run the pipeline: `./peakforge tsvmode samples.tsv --output-dir results` (or `python chipdiff.py tsvmode ...`).
3. Review tables, plots, and metadata under `results/`.

### Easier metadata generation (new)

If hand-writing CSV/TSV feels tedious, you can now generate it directly from run arguments:

```bash
./peakforge makesheet \
  --condition-a K562 \
  --a-bams example/data/K562_rep1.bam example/data/K562_rep2.bam \
  --a-peaks example/results/broad/K562_rep1_peaks.broadPeak example/results/broad/K562_rep2_peaks.broadPeak \
  --condition-b HepG2 \
  --b-bams example/data/HepG2_rep1.bam example/data/HepG2_rep2.bam \
  --b-peaks example/results/broad/HepG2_rep1_peaks.broadPeak example/results/broad/HepG2_rep2_peaks.broadPeak \
  --output samples.tsv
```

Then run as usual:

```bash
./peakforge tsvmode samples.tsv --output-dir results
```

`makesheet` supports both `.tsv` and `.csv` output paths and auto-generates unique sample IDs from BAM filenames.

### Sample sheet format

The sheet can be tab- or comma-delimited and must include `sample`, `condition`, and `bam`. Optional columns let you
supply pre-called peaks as MACS2 `narrowPeak` or `broadPeak` files.

#### 1 vs 1

| sample     | condition | bam                         | peaks                                               | peak_type |
|------------|-----------|-----------------------------|-------------------------------------|-----------|
| K562_rep1  | treated   | example/data/K562_rep1.bam  | example/results/narrow/K562_rep1_peaks.narrowPeak   | narrow    |
| HepG2_rep1 | control   | example/data/HepG2_rep1.bam | example/results/narrow/HepG2_rep1_peaks.narrowPeak  | narrow    |

#### 2 vs 2

| sample     | condition | bam                         | peaks                                               | peak_type |
|------------|-----------|-----------------------------|-------------------------------------|-----------|
| K562_rep1  | K562      | example/data/K562_rep1.bam  | example/results/broad/K562_rep1_peaks.broadPeak     | broad     |
| K562_rep2  | K562      | example/data/K562_rep2.bam  | example/results/broad/K562_rep2_peaks.broadPeak     | broad     |
| HepG2_rep1 | HepG2     | example/data/HepG2_rep1.bam | example/results/broad/HepG2_rep1_peaks.broadPeak    | broad     |
| HepG2_rep2 | HepG2     | example/data/HepG2_rep2.bam | example/results/broad/HepG2_rep2_peaks.broadPeak    | broad     |

### Example command
```bash
./peakforge tsvmode samples.tsv \
  --output-dir results \
  --peak-dir peaks \
  --peak-type narrow \
  --peak-extension 250 \
  --min-overlap 2 \
  --threads 16 \
  --gtf annotations.gtf \
  --enrichr
```
This run performs peak calling (if necessary), constructs consensus intervals, quantifies counts, executes the
appropriate differential workflow, annotates peaks, and produces summary plots. `--peak-type` controls whether MACS2 is
invoked in narrow- or broad-peak mode when peaks are missing; narrow peaks are expanded symmetrically by
`--peak-extension` (default 250 bp) to build the consensus, while broad peaks use their full width without extra
padding. Paired-end BAMs are detected automatically and passed to MACS2 with `BAMPE` format. If you supply `narrowPeak`
or `broadPeak` files in the sample sheet, PeakForge infers the mode from the extension and you can omit `--peak-type`
entirely. `--threads` applies to deepTools counting and to any `samtools index` calls needed to compute library sizes.

---

## Example workflows

The `example/` directory orchestrates a matched ENCODE MYC ChIP-seq dataset (hg38) containing downsampled BAMs.

### End-to-end 2 vs 2 analysis
```bash
bash example/run_pipeline.sh
```
Generates results under `example/results/`, including `consensus_peaks.bed` that can be reused with `--consensus-peaks`
in later runs. Internally the script executes the following PeakForge command (with `THREADS` defaulting to 16):

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
Skips MACS2 when peaks are provided and reuses consensus intervals for consistent comparisons. The helper ultimately
invokes `python chipdiff.py runmode` with the supplied arguments plus `--peak-type narrow --peak-extension 250
--min-overlap 2 --macs2-genome hs --threads 16`, so you can copy/paste the expanded command if you prefer to run
PeakForge directly.

### Prior-aware demo
```bash
bash example/run_with_prior.sh
```
Bootstraps priors, runs `tsvmode` with `--prior-manifest`, and performs peak-shape profiling with `--prior-shape`. The
manifest contains `prior_bed`, `prior_stats`, and `prior_weight` entries; the script then calls:

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
Any completed PeakForge analysis writes `consensus_peaks.bed` inside the output directory. Passing that file to either
`tsvmode` or `runmode` through `--consensus-peaks` preserves genomic intervals across follow-up contrasts, keeping
fold-change estimates directly comparable.

---

## Troubleshooting

- **Missing optional dependencies** – PyDESeq2, gseapy, and pyBigWig are only required for specific features. If the CLI warns about a missing module you can either install it (`pip install pydeseq2 gseapy pybigwig`) or run the pipeline without that capability.
- **External tool failures** – PeakForge wraps MACS2, deepTools, and samtools. Check their versions with `macs2 --version`, `multiBamSummary --version`, and `samtools --version` if a subprocess error occurs.
- **Peak shape inputs** – `peakshape` uses `bamCoverage` to convert BAMs to bigWig tracks. Install deepTools or provide `--bigwig-a/--bigwig-b` inputs if `bamCoverage` is not available.
- **Empty consensus sets** – ensure that MACS2 produced peaks or provide existing peak files via the metadata sheet. Lowering `--min-overlap` can help when combining sparse datasets.

---

## Command reference
Run `./peakforge --help` (or `python chipdiff.py --help`) to inspect CLI options including peak-calling parameters,
threading controls, annotation, enrichment, priors, and metadata-sheet helpers (`makesheet`).

---

## Optional priors (current behaviour)

Priors are used as **external covariates** for weighted multiple-testing correction, not as replacements for core differential statistics.

### What priors do (and do not do)
- Priors **do not** modify raw `log2FoldChange`.
- Priors **do not** modify raw `pvalue`.
- Priors **do** affect weighted BH outputs (`prior_weight`, `p_weighted`, `q_weighted`), which can change ranking and significance calls.

### Prior inputs and their roles
- `--prior-bed`: provides overlap/matching information and width context, producing covariates such as `PriorWeight`, `NoveltyPenalty`, and `WidthZ`.
- `--prior-bigwig`: provides reference signal intensity context, yielding per-peak `PriorIntensity` and `IntZ`.
- `--shape-tsv`: provides external peak-shape covariates aligned to consensus peaks (for example `ShapeZ`).
- `--prior-manifest`: bundles prior resources (commonly `prior_bed`, `prior_bigwig`, `prior_stats`, `prior_weight`) into one reusable config.

### Minimal prior usage example
```bash
./peakforge tsvmode example/data/metadata_1v1.tsv \
  --output-dir results/prior_demo \
  --prior-bed results/prior_demo/prior_inputs/demo_prior.bed \
  --prior-bigwig results/prior_demo/prior_inputs/demo_prior.bw \
  --shape-tsv results/prior_demo/shape/peak_shape_summary.tsv \
  --prior-weight 0.4
```

### When should I use priors?
Priors are most useful when external evidence exists for where true regulatory signal is likely to occur, especially in low-replicate or borderline-significance settings. Suitable priors should be as close as possible to the current analysis in assay type, cell or tissue context, genome build, and biological state.

### What makes a good prior?
The best priors are assay-matched and context-matched reference regions, such as public ATAC-seq peaks for ATAC-seq analyses, or matched CUT&Tag / ChIP peaks for the same target. DNase-seq peaks can also be used as a weaker prior for ATAC-seq because both assays capture open chromatin, but assay-matched ATAC priors are preferred.

### Can CUT&Tag priors be used for ATAC-seq?
Sometimes. CUT&Tag priors are most reasonable when the profiled mark or factor is expected to co-localize with accessible regulatory elements, such as active promoter/enhancer marks. They should be used cautiously and interpreted as functional-regulatory priors rather than generic accessibility priors.

### Should priors come from treatment or control samples?
In most cases, priors should be neutral and not strongly biased toward one comparison arm. We recommend using pooled, consensus, or independent reference peak sets rather than priors derived exclusively from only treatment or only control samples.

### What priors do not do
Priors in PeakForge do not modify the raw log2 fold-change or the raw p-value. Instead, priors are used as external covariates for weighted multiple-testing correction, which can improve ranking and prioritization of peaks supported by prior knowledge.

To ship a reusable prior bundle, place files next to a `prior_manifest.json` with keys such as `prior_bed`, `prior_bigwig`, `prior_stats`, and `prior_weight`. `example/run_with_prior.sh` demonstrates this layout end-to-end.

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
