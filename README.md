# PeakForge

PeakForge is a modern, easy-to-use toolkit for differential analysis of ATAC-seq, CUT&Tag, and ChIP-seq data. It handles everything from peak calling (wrapping MACS2/3) to consensus building, counting, and differential testing (using PyDESeq2 or MARS).

Designed to be as simple to use as `macs3`, PeakForge offers a unified CLI for all your analysis needs.

---

## Installation

### Requirements
- Python ≥ 3.10
- Standard bioinformatics tools: `samtools`, `macs3` (or `macs2`), `multiBamSummary` (from deepTools).

### Install
```bash
pip install numpy pandas scipy statsmodels matplotlib seaborn pyranges gseapy pydeseq2 macs3 deeptools typer rich
conda install -c bioconda samtools
```

### Verify
```bash
peakforge --help
```

---

## Quick Start

### 1. Run using a sample sheet (Recommended)
The easiest way to run PeakForge is with a simple metadata file (`samples.tsv`).

**samples.tsv:**
```tsv
sample      condition   bam
K562_rep1   K562        data/K562_rep1.bam
K562_rep2   K562        data/K562_rep2.bam
HepG2_rep1  HepG2       data/HepG2_rep1.bam
HepG2_rep2  HepG2       data/HepG2_rep2.bam
```

**Command:**
```bash
peakforge run samples.tsv --output-dir results
```
This single command will:
- Call peaks for each sample (if not provided).
- Build a consensus peak set.
- Count reads.
- Run differential analysis (DESeq2 for replicates, MARS for 1v1).
- Generate interactive reports and plots.

### 2. Run with direct arguments (Ad-hoc)
Quickly compare two conditions without creating a sample sheet.

```bash
peakforge diff \
  --condition-a K562  --a-bams data/K562_rep1.bam data/K562_rep2.bam \
  --condition-b HepG2 --b-bams data/HepG2_rep1.bam data/HepG2_rep2.bam \
  --output-dir results
```

---

## Commands

### `peakforge run`
Run the full pipeline from a sample sheet.

**Usage:**
```bash
peakforge run [OPTIONS] METADATA
```
**Key Options:**
- `--output-dir PATH`: Directory for results (default: `results`).
- `--peak-type [narrow|broad]`: Peak calling mode (default: `narrow`).
- `--threads INT`: Number of threads for counting and indexing (default: 16).
- `--enrichr`: Run GO enrichment analysis on top differential peaks.
- `--gtf PATH`: Gene annotation file (GTF) for peak annotation.

### `peakforge diff`
Run differential analysis by explicitly listing BAM files.

**Usage:**
```bash
peakforge diff --condition-a NAME --a-bams BAMS... --condition-b NAME --b-bams BAMS... [OPTIONS]
```
Same options as `run` are available for customization.

### `peakforge shape`
Analyze and compare peak shapes (metrics like FWHM, core:flank ratio) between two conditions.

**Usage:**
```bash
peakforge shape --bed regions.bed --bigwig-a A.bw --bigwig-b B.bw --out results/shape
```

---

## Output

All results are saved in the `--output-dir`:

- **`differential_results.tsv`**: The main result table with log2FC and adjusted p-values.
- **`plots/`**: Visualization including:
    - `volcano.png`: Volcano plot of significant changes.
    - `ma_plot.png`: MA plot (log2FC vs mean expression).
    - `sample_correlation.png`: Heatmap of sample clustering.
    - `top_peaks_heatmap.png`: Expression heatmap of top differential peaks.
- **`consensus_peaks.bed`**: The merged consensus peak set used for analysis.
- **`counts/`**: Raw read counts per peak.

---

## Advanced Usage

### Reusing Consensus Peaks
To ensure identical genomic intervals across multiple runs, generate `consensus_peaks.bed` once and reuse it:
```bash
peakforge run new_samples.tsv --consensus-peaks results/consensus_peaks.bed
```

### Prior-Informed Analysis
PeakForge supports integrating prior knowledge (e.g., ATAC-seq peaks or known binding sites) to improve sensitivity.
```bash
peakforge run samples.tsv --prior-bed priors.bed --prior-weight 0.4
```

---

## Troubleshooting

- **"Command not found: macs2/macs3"**: Ensure MACS is installed (`pip install macs3`).
- **"Command not found: multiBamSummary"**: Ensure deepTools is installed (`pip install deeptools`).
- **Missing Replicates**: PeakForge automatically switches to the MARS statistical test if no replicates are found (1 vs 1 comparison).

---

## References
- [PyDESeq2](https://github.com/owkin/PyDESeq2)
- [MACS3](https://github.com/macs3-project/MACS)
- [deepTools](https://deeptools.readthedocs.io/)
