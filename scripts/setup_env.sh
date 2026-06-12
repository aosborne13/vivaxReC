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

  echo "[vivaxReC] Installing R packages..."
  Rscript "${R_INSTALLER}"

  echo "[vivaxReC] Installing vivaxrec Python package..."
  pip install -e .

  echo "[vivaxReC] Setup complete. Run: micromamba activate ${ENV_NAME}"

elif command -v conda >/dev/null 2>&1; then
  echo "[vivaxReC] Using conda"
  conda env create -f "${YML_FILE}" || conda env update -f "${YML_FILE}"
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "${ENV_NAME}"

  echo "[vivaxReC] Installing R packages..."
  Rscript "${R_INSTALLER}"

  echo "[vivaxReC] Installing vivaxrec Python package..."
  pip install -e .

  echo "[vivaxReC] Setup complete. Run: conda activate ${ENV_NAME}"

else
  echo "[vivaxReC] ERROR: Neither micromamba nor conda found."
  echo "           Install micromamba: curl -L micro.mamba.pm/install.sh | bash"
  exit 1
fi
