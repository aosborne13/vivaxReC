options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("EPPIcenter/asterTES")
remotes::install_github("EPPIcenter/dcifer")
remotes::install_github("aimeertaylor/Pv3Rs")
