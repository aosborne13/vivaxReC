#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="vivaxReC"
YML_FILE="env/environment.yml"
R_INSTALLER="env/install_R_packages.R"

echo "[vivaxReC] Setting up environment: ${ENV_NAME}"

if command -v micromamba >/dev/null 2>&1; then
  echo "[vivaxReC] Using micromamba"

  # Enable micromamba activation in non-interactive scripts
  eval "$(micromamba shell hook -s bash)"

  echo "[vivaxReC] Creating environment from ${YML_FILE} ..."
  micromamba create -y -f "${YML_FILE}" || echo "[vivaxReC] Environment may already exist; continuing."

  echo "[vivaxReC] Activating environment..."
  micromamba activate "${ENV_NAME}"

elif command -v mamba >/dev/null 2>&1; then
  echo "[vivaxReC] Using mamba"

  echo "[vivaxReC] Creating environment from ${YML_FILE} ..."
  mamba env create -f "${YML_FILE}" || echo "[vivaxReC] Environment may already exist; continuing."

  echo "[vivaxReC] Activating environment..."
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "${ENV_NAME}"

elif command -v conda >/dev/null 2>&1; then
  echo "[vivaxReC] Using conda"

  echo "[vivaxReC] Creating environment from ${YML_FILE} ..."
  conda env create -f "${YML_FILE}" || echo "[vivaxReC] Environment may already exist; continuing."

  echo "[vivaxReC] Activating environment..."
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "${ENV_NAME}"

else
  echo "ERROR: micromamba, mamba, or conda not found."
  echo "Install micromamba (recommended) or conda/mamba, then retry."
  exit 1
fi

echo "[vivaxReC] Installing R packages from GitHub..."
Rscript "${R_INSTALLER}"

echo "[vivaxReC] Installing vivaxrec CLI..."
python -m pip install -e .

echo
echo "✅ vivaxReC install complete."
echo "Test with:"
echo "  vivaxrec --help"

