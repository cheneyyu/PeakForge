#!/usr/bin/env bash
set -euo pipefail

if [[ ! -f "pyproject.toml" ]]; then
  echo "Run this script from the PeakForge repository root." >&2
  exit 1
fi

export DEBIAN_FRONTEND=noninteractive
: "${UPGRADE_PIP_TOOLS:=0}"

apt-get update -y
apt-get install -y samtools

if [[ "${UPGRADE_PIP_TOOLS}" == "1" ]]; then
  python3 -m pip install --upgrade pip setuptools wheel
fi

python3 -m pip install --prefer-binary -e '.[macs3]'

echo
echo "PeakForge Colab setup completed."
echo "Try: peakforge --help"
echo "Verify with: macs3 --version"
if [[ "${UPGRADE_PIP_TOOLS}" != "1" ]]; then
  echo "Using Colab's preinstalled pip/setuptools/wheel."
  echo "If installation fails, rerun with: UPGRADE_PIP_TOOLS=1 bash colab/setup_colab.sh"
fi
