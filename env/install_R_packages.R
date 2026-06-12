#!/usr/bin/env Rscript
# =============================================================================
# install_R_packages.R
# Installs all R dependencies required by the vivaxReC pipeline.
# Called automatically by setup_env.sh.
# =============================================================================
message("[vivaxReC] Installing CRAN packages...")
cran_pkgs <- c(
  "readr", "dplyr", "stringr", "tidyr", "purrr", "rlang",
  "ggplot2", "patchwork", "scales", "forcats", "ggrepel",
  "tibble", "remotes"
)
install.packages(cran_pkgs, repos = "https://cloud.r-project.org", quiet = TRUE)

message("[vivaxReC] Installing Bioconductor packages...")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

message("[vivaxReC] Installing dcifer from CRAN...")
install.packages("dcifer", repos = "https://cloud.r-project.org", quiet = TRUE)

message("[vivaxReC] Installing Pv3Rs from GitHub...")
remotes::install_github("aimeertaylor/Pv3Rs", quiet = TRUE)

message("[vivaxReC] R package installation complete.")
