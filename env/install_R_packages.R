options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

required <- c("igraph", "partitions", "gmp")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing) > 0) {
  stop(
    paste0(
      "Missing R dependencies: ", paste(missing, collapse = ", "),
      "\n\nIf you are using the vivaxReC conda environment, update it with:\n",
      "  micromamba install -y -c conda-forge r-igraph r-partitions r-gmp\n",
      "or recreate the environment from env/environment.yml.\n"
    )
  )
}

remotes::install_github("EPPIcenter/asterTES")
remotes::install_github("EPPIcenter/dcifer")
remotes::install_github("aimeertaylor/Pv3Rs")
