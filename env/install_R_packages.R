options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# We manage dependencies via conda (env/environment.yml), so avoid source builds here
remotes::install_github("EPPIcenter/asterTES", dependencies = FALSE)
remotes::install_github("EPPIcenter/dcifer",  dependencies = FALSE)
remotes::install_github("aimeertaylor/Pv3Rs", dependencies = FALSE)
