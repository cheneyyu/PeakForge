# PeakForge

**PeakForge** is a Python-native, DiffBind-style toolkit for end-to-end ATAC-seq, CUT&Tag, and ChIP-seq differential analysis.
It takes BAM or MACS2 peak files as input, builds consensus peaks, counts reads, and runs differential analysis (PyDESeq2 or MARS) whether or not biological replicates are available.
Two ATAC contrasts searching for differential motifs is one of the core scenarios the pipeline targets, ensuring motif/regulon shifts are detectable even when each condition is represented by a single rich sample.
Beyond standard replicate-aware testing, the pipeline supports single-sample vs single-sample contrasts for cases where cohesive motif/regulon signals make no-replicate comparisons informative.

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
  - Reuse an existing BED with `--consensus-peaks` to keep genomic regions identical across runs.

- **Counting**
  - Uses **deepTools** `multiBamSummary BED-file` for counts matrix.
  - Exports counts matrix as `.tsv`.

- **Differential analysis**
  - Automatically routes to the appropriate workflow based on replicate availability.
  - **DESeq2 (PyDESeq2)** (with replicates):
    - Median-ratio normalization, dispersion estimation, and Wald testing handled by [PyDESeq2](https://github.com/owkin/PyDESeq2).
    - Requires PyDESeq2 to be installed when replicate designs are detected.
  - **MARS (DEGseq, Likun Wang 2010)** (without replicates):
    - Supports 1 vs 1, or pooled multiple vs multiple samples.
    - Uses `samtools idxstats` to derive library sizes from the full BAM before applying the MA-plot random sampling test.
  - Consolidated differential results (`differential_results.tsv`) for downstream interpretation.

- **Annotation & Enrichment (optional)**
  - Annotate consensus peaks with nearest genes (via GTF).
  - Run GO BP enrichment using **gseapy** (Enrichr).

- **Output**
  - Volcano plots (PyDESeq2 and MARS).
  - MA plots.
  - Sample correlation heatmap.
  - Heatmap of top differential peaks.
  - Results tables (`.tsv`) and metadata (`.json`).
  - Standalone peak shape profiling via `peak_shape.py` for comparing two bigWig
    tracks over a BED of regions (with delta metrics and plots).

---

## 📦 Installation

Requirements:

- Python ≥ 3.9
- Libraries: `numpy pandas scipy statsmodels matplotlib seaborn pyranges gseapy pydeseq2`
- External tools:
  - [MACS2](https://github.com/macs3-project/MACS) (for peak calling)
  - [deepTools](https://deeptools.readthedocs.io/en/develop/) (for multiBamSummary)
  - [samtools](http://www.htslib.org/) (for BAM indexing and library size estimation)

Install dependencies:

```bash
pip install numpy pandas scipy statsmodels matplotlib seaborn pyranges gseapy pydeseq2
conda install -c bioconda macs2 deeptools samtools
```

---

## 🧪 Quick start

1. Prepare a sample sheet (`samples.tsv`) describing your BAM files and optional peak calls.
2. Run the pipeline with `./peakforge tsvmode samples.tsv --output-dir results` (or `python chipdiff.py tsvmode samples.tsv --output-dir results`).
3. Inspect the summary plots and result tables written to the `results/` directory.

### Sample sheet format

The sheet can be tab- or comma-delimited and must include the columns `sample`, `condition`, and `bam`. Optional columns allow you to reference pre-called peaks.

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

This run will call peaks (if needed), build consensus peaks across samples, compute counts, perform PyDESeq2-based differential analysis when replicates are present (falling back to MARS otherwise), optionally annotate peaks, and generate standard comparison plots.

`multiBamSummary` is invoked with `--numberOfProcessors 16` by default; override with `--threads` if you need a different level of parallelism.

### Output overview

Key files generated under `results/` include:

- `counts.tsv` – consensus peak counts matrix.
- `differential_results.tsv` – unified differential analysis table (method annotated per peak).
- `plots/` – volcano, MA, correlation, and heatmap figures.
- `plots/differential_summary.png` – clusterProfiler-style overview of significant peaks.
- `metadata.json` – run configuration and provenance metadata.

---

## 🔍 Peak shape analysis module

The repository also ships a standalone shape-profiling utility for comparing
two signal tracks over a shared set of genomic regions:

```bash
python peak_shape.py \
  --bigwig-a sampleA.bw \
  --bigwig-b sampleB.bw \
  --bed peaks.bed \
  --core 500 \
  --flank 1000 3000 \
  --out results/shape
```

For each interval the script normalises the signal, computes FWHM, core:flank
ratios, centroid shifts, and skewness, then records the per-sample values plus
their deltas in `peak_shape.tsv`.  Summary histograms, scatter plots, and a
heatmap of the top outliers are written to `results/shape/plots/`.

To try the module without downloading real datasets, use the synthetic example
under `example/peak_shape/`:

```bash
cd example/peak_shape
bash run_peak_shape.sh
```

The helper script builds toy bigWig tracks on demand so you can inspect the
workflow end-to-end.  Replace the generated files with your own bigWigs and BED
to run the same analysis on real data.

---

## 🧬 ENCODE hg38 2v2 example dataset

To help you get started quickly, the repository ships with an end-to-end example based on a matched ENCODE MYC ChIP-seq dataset (hg38):

- **K562 MYC** – replicates [`ENCFF975ETI.bam`, `ENCFF380OWL.bam`]
- **HepG2 MYC** – replicates [`ENCFF315AUW.bam`, `ENCFF987GJQ.bam`]

The scripts in `example/` orchestrate downloading the public alignments and executing the PeakForge pipeline for the 2 vs 2 comparison. The dataset comprises four compact (downsampled) BAM files; make sure all four are present before running `run_pipeline.sh` so the analysis has the expected inputs.

1. **Prepare the inputs**
   ```bash
   # Prints the curl commands by default (DRY_RUN=1)
   bash example/download_encode.sh

   # Actually download the BAMs (~1.5 GB total)
   DRY_RUN=0 bash example/download_encode.sh
   ```
   Downloads are written to `example/data/` and a manifest (`encode_manifest.tsv`) is generated alongside the tab-delimited sample sheet (`metadata.tsv`).
   When `samtools` is available on `PATH`, the script also creates `.bai` indices for each BAM so the downstream counting step can stream efficiently. (If `samtools` is missing a warning is emitted.)

2. **Run the complete analysis**
   ```bash
   bash example/run_pipeline.sh
   ```
   This executes the 2v2 workflow and stores results in `example/results/`.

All scripts respect relative paths, so you can copy the `example/` directory into your own project and customize it as needed.

### Example 1: 1v1 quick start

To exercise the no-replicate branch of the pipeline, run the bundled 1 vs 1 script:

```bash
bash example/run_example_1v1.sh
```

This uses the `metadata_1v1.tsv` sheet to compare `K562_rep1` against `HepG2_rep1` and produces MARS differential results under `example/results_1v1/`.

### Example 2: 1v1 with existing BAM/peak files

If you already have a pair of BAMs (and optionally MACS2 peak calls) on disk,
the helper script `example/run_example2.sh` lets you drive the 1 vs 1 branch
directly without re-running the download/index step:

```bash
bash example/run_example2.sh \
  --condition-a K562 \
  --a-bams example/data/K562_rep1.bam \
  --a-peaks example/results/2v2/peaks/K562_rep1_summits.bed \
  --condition-b HepG2 \
  --b-bams example/data/HepG2_rep1.bam \
  --b-peaks example/results/2v2/peaks/HepG2_rep1_summits.bed \
  --consensus-peaks example/results/2v2/consensus_peaks.bed
```

The script calls `peakforge runmode` with the provided paths and sensible
defaults (including `--threads 16`, which maps to `multiBamSummary
--numberOfProcessors`). Provide peak files to skip MACS2 entirely; omit them if
you want the pipeline to call peaks from your BAMs on the fly. Supplying
`--consensus-peaks` reuses the same peak set produced by the 2 vs 2 workflow so
the MARS comparison stays on the identical genomic intervals.

### Example 3: 2v2 quick start

Once the ENCODE dataset is downloaded, you can launch the bundled 2 vs 2
workflow via:

```bash
bash example/run_pipeline.sh
```

This script reads `example/data/metadata.tsv`, runs the full PeakForge pipeline
with a default of 16 `multiBamSummary` threads, and stores results under
`example/results/`. The generated consensus (`example/results/2v2/consensus_peaks.bed`)
can be passed to future `peakforge` runs with `--consensus-peaks` for consistent
peak definitions.

### Example 4: 2v2 with existing BAM/peak files

When your own paired conditions (with replicates) are already aligned, reuse
`example/run_example2.sh` to spin up the full differential analysis without
taking the download shortcut:

```bash
bash example/run_example2.sh \
  --condition-a K562 \
  --a-bams example/data/K562_rep1.bam example/data/K562_rep2.bam \
  --a-peaks example/results/2v2/peaks/K562_rep1_summits.bed example/results/2v2/peaks/K562_rep2_summits.bed \
  --condition-b HepG2 \
  --b-bams example/data/HepG2_rep1.bam example/data/HepG2_rep2.bam \
  --b-peaks example/results/2v2/peaks/HepG2_rep1_summits.bed example/results/2v2/peaks/HepG2_rep2_summits.bed
```

Supplying peak calls for each replicate skips MACS2 entirely; otherwise the
script will trigger peak calling for the provided BAMs before consensus/DE
analysis.

### Reusing consensus peaks across runs

Any completed PeakForge analysis produces `consensus_peaks.bed` inside its
output directory. Pass that file to either `tsvmode` or `runmode` with
`--consensus-peaks` to reuse the exact same genomic intervals in follow-up
contrasts (for example, when you want to compare a single replicate against the
2 vs 2 consensus). Doing so avoids re-running MACS2 and keeps fold-change
estimates directly comparable across analyses.

## 📝 TODO

- Add a gallery of publication-level figures that showcases typical PeakForge outputs.

---

## 🔧 Command reference

Run `./peakforge --help` (or `python chipdiff.py --help`) to see all available options (peak calling parameters, threading, annotation, and enrichment settings).
