# PeakForge

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/cheneyyu/PeakForge/blob/main/colab/PeakForge_Colab_Quickstart.ipynb)

PeakForge is a Python-native toolkit for differential analysis of ATAC-seq, CUT&Tag, and ChIP-seq peak data. It supports both replicate-aware comparisons and data-limited `1 vs 1` analyses, starting from BAM files or pre-called peaks and producing consensus peaks, count matrices, differential results, and standard diagnostic plots.

This directory is a GitHub-ready distribution of the current PeakForge codebase. Large BAM example files are intentionally excluded so the repository stays lightweight and uploadable. Example scripts and the Colab notebook download the public ENCODE MYC ChIP-seq BAMs together with the matched input-control BAMs on demand.

## Google Colab

For most users, the fastest way to try PeakForge is the Colab notebook:

```text
https://colab.research.google.com/github/cheneyyu/PeakForge/blob/main/colab/PeakForge_Colab_Quickstart.ipynb
```

The notebook reproduces the main user-facing workflow in a lightweight environment:

- installs PeakForge with `samtools` and `macs3`
- verifies the CLI
- downloads six ENCODE MYC ChIP-seq BAMs plus two matched input-control BAMs
- runs the standard input-aware `1 vs 1` and `2 vs 2` ChIP-seq analyses
- reproduces the overlap-based ROC comparison used in the manuscript

Colab is suitable for demos, tutorials, and small analyses. The setup cell usually takes about `3-8 min`. For longer runs or larger BAM collections, use a local environment instead.

## Local Installation

For a full local setup, clone the repository and run the installation commands together:

```bash
git clone https://github.com/cheneyyu/PeakForge.git
cd PeakForge
curl -LsSf https://astral.sh/uv/install.sh | sh   # only if uv is not installed yet
sudo apt-get update -y
sudo apt-get install -y samtools
uv sync --extra macs3
uv run peakforge --help
samtools --version
macs3 --version
```

PeakForge requires `samtools` plus a working `macs2` or `macs3` on `PATH`, and it prefers `macs3` when both are available. `uv sync --extra macs3` installs the Python package together with a Python-managed `macs3`, which is the recommended path for most users.

If you already manage bioinformatics tools through Conda or a system package manager, that is also fine as long as `samtools` and `macs2` or `macs3` run successfully from the shell. Optional downstream tools such as `HOMER` are only needed for motif analysis outside the core pipeline.

## What Is Included

- `chipdiff.py`, `io_utils.py`, `motif_ranking.py`, `peak_shape.py`
- `pyproject.toml` for `uv`-based environment management
- `example/` scripts for:
  - ENCODE `2 vs 2` replicate-aware benchmark
  - ENCODE `1 vs 1` no-replicate benchmark
  - optional `3-fold` held-out validation
  - peak-shape demo

## Quick start

### Show the CLI

```bash
uv run peakforge --help
```

### Run the ENCODE `2 vs 2` example

```bash
DRY_RUN=0 bash example/download_encode.sh
bash example/run_pipeline.sh
```

Outputs will be written to:

```text
example/results/2v2/
```

### Run the ENCODE `1 vs 1` example

```bash
DRY_RUN=0 bash example/download_encode.sh
bash example/run_example_1v1.sh
```

Outputs will be written to:

```text
example/results_1v1/
```

### Run the held-out `3-fold` validation

Download the extra third replicate pair as well:

```bash
DRY_RUN=0 INCLUDE_THIRD_REPLICATES=1 bash example/download_encode.sh
bash example/run_lopo_validation.sh
```

Outputs will be written to:

```text
example/results_3v3/
```

## Main commands

### Metadata-driven mode

```bash
uv run peakforge tsvmode samples.tsv --output-dir results
```

### Direct argument mode

```bash
uv run peakforge runmode \
  --condition-a K562 \
  --a-bams path/to/K562_rep1.bam path/to/K562_rep2.bam \
  --a-controls path/to/K562_input.bam \
  --condition-b HepG2 \
  --b-bams path/to/HepG2_rep1.bam path/to/HepG2_rep2.bam \
  --b-controls path/to/HepG2_input.bam \
  --output-dir results
```

For ATAC-seq or CUT&Tag, omit the control arguments.

### Generate a sample sheet

```bash
uv run peakforge makesheet \
  --condition-a K562 \
  --a-bams path/to/K562_rep1.bam path/to/K562_rep2.bam \
  --a-controls path/to/K562_input.bam \
  --condition-b HepG2 \
  --b-bams path/to/HepG2_rep1.bam path/to/HepG2_rep2.bam \
  --b-controls path/to/HepG2_input.bam \
  --output samples.tsv
```

### Peak-shape profiling

```bash
uv run peakforge peakshape --help
```

## Repository layout

```text
.
├── colab/
├── chipdiff.py
├── io_utils.py
├── motif_ranking.py
├── peak_shape.py
├── pyproject.toml
└── example/
    ├── README.md
    ├── download_encode.sh
    ├── run_pipeline.sh
    ├── run_example_1v1.sh
    ├── run_example2.sh
    ├── run_lopo_validation.sh
    └── peak_shape/
```

## Notes

- The example BAM files are not committed here because they exceed normal GitHub-friendly repository size.
- `example/download_encode.sh` fetches the public ENCODE MYC ChIP-seq BAMs together with the matched K562 and HepG2 input-control BAMs required by the standard ChIP-seq examples.
- `1 vs 1` mode is supported and useful for ranking and exploratory follow-up, but replicate-supported analysis remains the primary validation setting.

## License

PeakForge is distributed under the license provided in [LICENSE](LICENSE).
