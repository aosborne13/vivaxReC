#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="vivaxReC"
YML_FILE="env/environment.yml"
R_INSTALLER="env/install_R_packages.R"

echo "[vivaxReC] Setting up environment: ${ENV_NAME}"

# Pick solver
if command -v mamba >/dev/null 2>&1; then
  SOLVER="mamba"
elif command -v conda >/dev/null 2>&1; then
  SOLVER="conda"
else
  echo "ERROR: conda or mamba not found. Install Miniconda/Mambaforge/Micromamba first."
  exit 1
fi

# Create env (if it already exists, conda/mamba will error; we handle that)
echo "[vivaxReC] Creating conda environment from ${YML_FILE} ..."
${SOLVER} env create -f "${YML_FILE}" || echo "[vivaxReC] Environment may already exist; continuing."

# Activate env (works for conda/mamba)
echo "[vivaxReC] Activating environment..."
source "$(${SOLVER} info --base)/etc/profile.d/conda.sh"
conda activate "${ENV_NAME}"

# Install R packages (Pv3Rs, asterTES, dcifer)
echo "[vivaxReC] Installing R packages from GitHub..."
Rscript "${R_INSTALLER}"

# Install Python CLI (editable)
echo "[vivaxReC] Installing vivaxrec CLI..."
python -m pip install -e .

echo
echo "✅ vivaxReC install complete."
echo "Test with:"
echo "  vivaxrec --help"
