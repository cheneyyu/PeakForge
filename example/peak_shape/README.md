# Peak shape demo

This directory provides a self-contained walkthrough for the
`peak_shape.py` module.  It builds miniature synthetic bigWig tracks so the
workflow can be exercised without downloading large public datasets.

## Contents

- `make_demo_tracks.py` – generates toy bigWig files and a matching BED file.
- `run_peak_shape.sh` – wrapper script that ensures the demo data exist and runs
  the analysis with reasonable defaults.
- `data/` – created on demand; contains the generated bigWigs and BED regions.
- `results/` – output directory written by `peak_shape.py`.

## Usage

1. **Create the demo tracks (optional)**

   ```bash
   python make_demo_tracks.py
   ```

   The script is idempotent and can be re-run to refresh the inputs.  The
generated files land in `data/` under this directory.

2. **Execute the peak shape analysis**

   ```bash
   bash run_peak_shape.sh
   ```

   Set the `THREADS` environment variable to control the worker pool size:

   ```bash
   THREADS=4 bash run_peak_shape.sh
   ```

   Results (TSV and plots) are written to `results/`.

## Customising the example

Replace `data/demo_sample_*.bw` and `data/demo_regions.bed` with your own tracks
and regions to visualise real datasets.  The helper script simply ensures the
input files exist; if you supply your own data, skip the generation step and run
`run_peak_shape.sh` directly.
