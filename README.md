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
- Incorporates public resources (e.g. ENCODE, Roadmap Epigenomics) to regularise peak width, intensity, and shape metrics.
- Accepts priors via `--prior-bed`, `--prior-bigwig`, `--prior-manifest`, or `--prior-shape` with tunable weights.
- Peak shape profiling can reuse priors through the standalone `peak_shape.py` module or the `peakforge peakshape` subcommand.

### Outputs and visualisation
- Volcano plots, MA plots, sample correlation heatmaps, and top-peak heatmaps.
- JSON metadata capturing run configuration and summary statistics.
- Peak-shape delta metrics and plots when the dedicated subcommand is invoked.

---

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
