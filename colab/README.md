# Google Colab Support

PeakForge can run in Google Colab for lightweight tutorials, small user datasets, and smoke tests. The most practical Colab target is:

- `runmode` on a small number of BAM files
- analyses that reuse pre-called peaks or an existing consensus BED
- peak-shape profiling on small inputs

Colab is not the best environment for the full public ENCODE benchmark because:

- the example BAM files are large
- Colab VM disk and session lifetime are limited
- repeated indexing and peak calling can be slow on ephemeral machines

The Colab notebook intentionally spells out the ENCODE downloads directly in notebook cells, rather than delegating them to `example/download_encode.sh`, so the workflow stays transparent for Colab users.

## Recommended Colab workflow

Clone the repository and run the setup commands from the repository root:

```bash
!git clone https://github.com/cheneyyu/PeakForge.git
%cd PeakForge
!apt-get update -y
!apt-get install -y samtools
!python3 -m pip install --prefer-binary -e '.[macs3]'
```

Or open the notebook directly in Colab:

```text
colab/PeakForge_Colab_Quickstart.ipynb
```

This installs:

- `samtools` with `apt`
- `macs3` through the Python environment
- the PeakForge Python package and its declared dependencies with `pip`
- without upgrading Colab's preinstalled `pip`, `setuptools`, or `wheel` unless you explicitly request it

On a typical Colab runtime, this setup step usually takes about `3-8 min`.

This Colab path intentionally avoids Ubuntu's `macs` package because recent Colab images can expose a broken `macs2` launcher.

If the installation fails because of a packaging issue in a future Colab runtime, retry with:

```bash
!python3 -m pip install --upgrade pip setuptools wheel
!python3 -m pip install --prefer-binary -e '.[macs3]'
```

After setup, verify the environment:

```bash
!peakforge --help
!samtools --version
!macs3 --version
```

The default notebook follows the standard input-aware ChIP-seq path: it downloads six ENCODE MYC BAMs plus the two matched input-control BAMs and uses those controls in the example `1 vs 1` and `2 vs 2` runs. The manuscript's no-input sensitivity rerun is not part of the default Colab workflow.

## Minimal Colab example

If you already have BAM and peak files in the Colab working directory or on Google Drive, a small direct run looks like:

```bash
!peakforge runmode \
  --condition-a Tumor \
  --a-bams data/tumor.bam \
  --a-controls data/tumor_input.bam \
  --a-peaks data/tumor_peaks.narrowPeak \
  --condition-b Nearby \
  --b-bams data/nearby.bam \
  --b-controls data/nearby_input.bam \
  --b-peaks data/nearby_peaks.narrowPeak \
  --output-dir colab_results \
  --peak-type narrow \
  --threads 2
```

For ATAC-seq or CUT&Tag, omit `--a-controls` and `--b-controls`.

If you already have a consensus BED, you can reuse it:

```bash
!peakforge runmode \
  --condition-a Tumor \
  --a-bams data/tumor.bam \
  --condition-b Nearby \
  --b-bams data/nearby.bam \
  --output-dir colab_results \
  --consensus-peaks data/consensus_peaks.bed \
  --threads 2
```

## Practical recommendations

- Start with `--threads 2` or `--threads 4`; Colab does not benefit much from the desktop default of `16`.
- On Colab, prefer `pip install -e '.[macs3]'` and ignore `macs2`.
- Prefer pre-called peaks or a prebuilt consensus set when possible.
- Keep BAM files on Google Drive if you need to preserve them across sessions.
- Avoid downloading the full example benchmark unless you specifically want to test the end-to-end workflow on a transient VM.
