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

## 🧬 ENCODE hg38 2v2 example dataset

To help you get started quickly, the repository ships with an end-to-end example based on a released ENCODE MYC ChIP-seq dataset generated in the HepG2 cell line (GRCh38 assembly):

- **Experiment A (MYC ChIP)** – ENCFF975ETI, ENCFF380OWL
- **Experiment B (MYC ChIP)** – ENCFF315AUW, ENCFF987GJQ

The two experiments capture independent MYC ChIP-seq replicates on the same HepG2 background, making them ideal for demonstrating PeakForge's replicate-aware consensus building and reproducibility reporting.

The scripts in `example/` orchestrate downloading the public alignments, executing the PeakForge pipeline for the 2 vs 2 comparison, repeating all four possible 1 vs 1 contrasts, and benchmarking how closely the 1 vs 1 runs reproduce the 2 vs 2 signal.

1. **Prepare the inputs**
   ```bash
   # Prints the curl commands by default (DRY_RUN=1)
   bash example/download_encode.sh

   # Actually download the BAMs (~11 GB total)
   DRY_RUN=0 bash example/download_encode.sh
   ```
   Downloads are written to `example/data/` and a manifest (`encode_manifest.tsv`) is generated alongside the tab-delimited sample sheet (`metadata.tsv`).

2. **Run the complete analysis**
   ```bash
   bash example/run_pipeline.sh
   ```
   This executes the 2v2 workflow plus four one-vs-one runs, storing results in `example/results/`.

3. **Inspect reproducibility reports**
   `example/analyze_replicates.py` (invoked automatically by `run_pipeline.sh`) aggregates:
   - peak-level overlap precision/recall, F1, bp-wise Jaccard, and sign concordance;
   - top-N recovery of the most significant 2v2 peaks;
   - Spearman correlations of log2 fold-changes between the 2v2 run and each 1v1 replicate pairing;
   - a union log2FC matrix (`global_log2fc_matrix.tsv`) for downstream clustering/QC.

All scripts respect relative paths, so you can copy the `example/` directory into your own project and customize it as needed.

---

## 🔧 Command reference

Run `python chipdiff.py --help` to see all available options (peak calling parameters, threading, annotation, and enrichment settings).
