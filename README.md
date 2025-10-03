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
  - **DESeq2-like**:
    - Median-ratio size factor normalization.
    - Dispersion trend fitting + empirical Bayes shrinkage.
    - Wald test with LFC shrinkage (ridge-EB).
  - **MARS (DEGseq, Likun Wang 2010)**:
    - Supports 1 vs 1, or pooled multiple vs multiple samples.
    - MA-plot based exact binomial test.

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
