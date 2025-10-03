# PeakForge

**PeakForge** is a publishable-level Python pipeline for end-to-end CUT&Tag / ChIP-seq differential analysis.  
It takes BAM or MACS2 peak files as input, builds consensus peaks, counts reads, runs differential analysis (DESeq2-like and MARS), and produces publication-quality figures.

---

## 🚀 Features

- **Input**
  - Accepts BAM files (paired-end or single-end).
  - Accepts existing MACS2 results: `summits.bed`, `narrowPeak`, or `broadPeak`.
  - If only BAMs are provided, calls peaks automatically via **MACS2** (supports narrow or broad peaks).

- **Consensus peaks**
  - Supports **minOverlap** (e.g. ≥2 samples required).
  - Summit/narrow peaks: extend ±250bp windows (adjustable).
  - Broad peaks: merged directly.

- **Counting**
  - Uses **deepTools** `multiBamSummary BED-file` for counts matrix.
  - Exports counts matrix as `.tsv`.

- **Differential analysis**
  - Automatically routes to the appropriate workflow based on replicate availability.
  - **DESeq2-like** (with replicates):
    - Median-ratio size factor normalization.
    - Dispersion trend fitting + empirical Bayes shrinkage.
    - Wald test with LFC shrinkage (ridge-EB).
  - **MARS (DEGseq, Likun Wang 2010)** (without replicates):
    - Supports 1 vs 1, or pooled multiple vs multiple samples.
    - MA-plot based exact binomial test.
  - Consolidated differential results (`differential_results.tsv`) for downstream interpretation.

- **Annotation & Enrichment (optional)**
  - Annotate consensus peaks with nearest genes (via GTF).
  - Run GO BP enrichment using **gseapy** (Enrichr).

- **Output**
  - Volcano plots (DESeq2-like and MARS).
  - MA plots.
  - Sample correlation heatmap.
  - Heatmap of top differential peaks.
  - Results tables (`.tsv`) and metadata (`.json`).

---

## 📦 Installation

Requirements:

- Python ≥ 3.9
- Libraries: `numpy pandas scipy statsmodels matplotlib seaborn pyranges gseapy`
- External tools:
  - [MACS2](https://github.com/macs3-project/MACS) (for peak calling)
  - [deepTools](https://deeptools.readthedocs.io/en/develop/) (for multiBamSummary)

Install dependencies:

```bash
pip install numpy pandas scipy statsmodels matplotlib seaborn pyranges gseapy
conda install -c bioconda macs2 deeptools
```

---

## 🧪 Quick start

1. Prepare a sample sheet (`samples.tsv`) describing your BAM files and optional peak calls.
2. Run the pipeline with `python chipdiff.py --metadata samples.tsv --output-dir results`.
3. Inspect the figures and result tables written to the `results/` directory.

### Sample sheet format

The sheet can be tab- or comma-delimited and must include the columns `sample`, `condition`, and `bam`. Optional columns allow you to reference pre-called peaks.

| sample | condition | bam             | peaks                 | peak_type |
|--------|-----------|-----------------|-----------------------|-----------|
| S1     | treated   | data/S1.bam     | data/S1_summits.bed   | summit    |
| S2     | treated   | data/S2.bam     | -                     | -         |
| C1     | control   | data/C1.bam     | data/C1_peaks.bed     | narrow    |

### Example command

```bash
python chipdiff.py \
  --metadata samples.tsv \
  --output-dir results \
  --peak-dir peaks \
  --min-overlap 2 \
  --gtf annotations.gtf \
  --enrichr
```

This run will call peaks (if needed), build consensus peaks across samples, compute counts, perform DESeq2-like differential analysis, optionally annotate peaks, and produce publication-ready plots.

### Output overview

Key files generated under `results/` include:

- `counts.tsv` – consensus peak counts matrix.
- `differential_results.tsv` – unified differential analysis table (method annotated per peak).
- `plots/` – volcano, MA, correlation, and heatmap figures.
- `plots/differential_summary.png` – clusterProfiler-style overview of significant peaks.
- `metadata.json` – run configuration and provenance metadata.

---

## 🔧 Command reference

Run `python chipdiff.py --help` to see all available options (peak calling parameters, threading, annotation, and enrichment settings).
