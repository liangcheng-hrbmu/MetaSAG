#!/usr/bin/env Rscript

# Install R packages that are not reliably available from conda.
#
# Run after creating and activating the MetaSAG conda environment:
#   conda activate metasag
#   Rscript setup_r.R

message("R executable: ", file.path(R.home("bin"), "R"))
message("R version: ", as.character(getRversion()))
message("R library paths:")
message(paste("  -", .libPaths(), collapse = "\n"))

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org")
}

if (!requireNamespace("sankeyD3", quietly = TRUE)) {
  message("Installing sankeyD3 from GitHub...")
  remotes::install_github("fbreitwieser/sankeyD3", upgrade = "never")
}

suppressPackageStartupMessages({
  library(sankeyD3)
})

message("sankeyD3 installation and loading test passed.")
