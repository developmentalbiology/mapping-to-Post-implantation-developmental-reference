#!/usr/bin/env Rscript

# ========================================================
# Install all R packages required for spatial annotation analysis
# Including packages from CRAN and GitHub
# ========================================================

cat("🚀 Starting installation of R packages required for spatial annotation analysis...\n")

# Set CRAN mirror (recommended to use official or domestic mirror)
options(repos = c(CRAN = "https://cloud.r-project.org/"))
options(timeout = 300) 
# List of packages from CRAN
cran_packages <- c(
  "Seurat",
  "readr",
  "tidyverse",
  "Matrix",
  "data.table",
  "doParallel",
  "parallelly",      # Note: you wrote "parallelly", not "parallel" ("parallel" is a base R package)
  "quadprog",
  "yaml",            # Note: the CRAN package name for r-yaml is simply "yaml"
  "ggplot2",
  "RColorBrewer",
  "viridis",
  "patchwork",
  "remotes"          # Must install "remotes" first to install GitHub packages
)

# Install CRAN packages
cat("📦 Installing CRAN packages...\n")
install.packages(cran_packages, dependencies = TRUE)

# Install GitHub packages
cat("🐙 Installing additional packages from GitHub...\n")

# Install argparse (trevorld/r-argparse)
remotes::install_github("trevorld/r-argparse")

# Install spacexr (dmcable/spacexr)
remotes::install_github("dmcable/spacexr", build_vignettes = FALSE)

cat("\n🎉 All packages have been installed successfully!\n")
cat("You can now run the spatial annotation analysis workflow.\n")
