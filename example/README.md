# Example workflows

This directory contains public-data examples for the GitHub-ready PeakForge distribution.

## What is included

- `download_encode.sh`: download the public ENCODE MYC ChIP-seq BAMs plus matched input-control BAMs
- `run_pipeline.sh`: replicate-aware `2 vs 2` benchmark
- `run_example_1v1.sh`: no-replicate `1 vs 1` benchmark
- `run_example2.sh`: direct `runmode` example using BAMs and optional peak files
- `run_lopo_validation.sh`: optional `3-fold` held-out validation using a third replicate pair
- `peak_shape/`: synthetic peak-shape demo

## Before running

From the repository root:

```bash
uv sync --extra macs3
```

Also make sure `samtools` is available:

```bash
samtools --version
```

## Example 1: download the public ENCODE files

Preview the commands first:

```bash
bash example/download_encode.sh
```

Actually download the standard input-aware ChIP-seq example files:

```bash
DRY_RUN=0 bash example/download_encode.sh
```

To fetch the additional third replicate pair used by the held-out validation:

```bash
DRY_RUN=0 INCLUDE_THIRD_REPLICATES=1 bash example/download_encode.sh
```

## Example 2: run the replicate-aware benchmark

```bash
bash example/run_pipeline.sh
```

This is the standard input-aware `2 vs 2` ChIP-seq example. The script reads `control_bam` from `example/data/metadata.tsv`.

Expected output:

```text
example/results/2v2/
```

## Example 3: run the `1 vs 1` fallback

```bash
bash example/run_example_1v1.sh
```

This example is also input-aware by default and reads `control_bam` from `example/data/metadata_1v1.tsv`.

Expected output:

```text
example/results_1v1/
```

## Example 4: run the held-out validation

```bash
DRY_RUN=0 INCLUDE_THIRD_REPLICATES=1 bash example/download_encode.sh
bash example/run_lopo_validation.sh
```

Expected output:

```text
example/results_3v3/
```

The manuscript's no-input sensitivity check is not part of the default GitHub or Colab workflow.

## Example 5: direct `runmode`

```bash
bash example/run_example2.sh \
  --condition-a K562 \
  --a-bams example/data/K562_rep1.bam example/data/K562_rep2.bam \
  --a-controls example/data/K562_input.bam \
  --condition-b HepG2 \
  --b-bams example/data/HepG2_rep1.bam example/data/HepG2_rep2.bam \
  --b-controls example/data/HepG2_input.bam
```

## Data layout

- `example/data/` contains the main `2 vs 2` benchmark inputs, matched input controls, and metadata
- `example/data2/` is reserved for the extra third replicate pair used in `3-fold` validation

The large BAMs are not stored in git. They are downloaded on demand.
