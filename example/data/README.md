# Example data

This directory stores metadata files and, after download, the ENCODE MYC ChIP-seq BAMs and matched input-control BAMs for the core PeakForge examples.

Tracked files:

- `metadata.tsv`
- `metadata_1v1.tsv`
- `metadata_3v3.tsv`
- `encode_manifest.tsv`
- `encode_manifest_3v3.tsv`

Not tracked in git:

- `*.bam`
- `*.bai`

To download the public BAMs:

```bash
DRY_RUN=0 bash ../download_encode.sh
```

The generated metadata tables use the `control_bam` column for the standard input-aware ChIP-seq examples.

To also download the third replicate pair used by LOPO validation:

```bash
DRY_RUN=0 INCLUDE_THIRD_REPLICATES=1 bash ../download_encode.sh
```
