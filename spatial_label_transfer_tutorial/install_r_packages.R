#!/usr/bin/env Rscript

# Install R packages required for the spatial annotation workflow
# 
# This script installs all the required R packages for the workflow.
# Run this script before using the workflow for the first time.

cat("Installing R packages for Spatial Annotation Workflow...\n")

# List of required packages
required_packages <- c(
  "Seurat",
  "Matrix", 
  "readr",
  "data.table",
  "tidyverse",
  "spacexr",
  "doParallel",
  "parallel",
  "quadprog",
  "argparse",
  "yaml",
  "ggplot2",
  "RColorBrewer",
  "viridis",
  "patchwork"
)

# Function to install packages if not already installed
install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat("Installing", pkg, "...\n")
      install.packages(pkg, dependencies = TRUE)
      
      # Check if installation was successful
      if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        cat("Failed to install", pkg, "\n")
      } else {
        cat("Successfully installed", pkg, "\n")
      }
    } else {
      cat(pkg, "is already installed\n")
    }
  }
}

# Install CRAN packages
cat("Installing CRAN packages...\n")
install_if_missing(required_packages)

# Install spacexr from GitHub if not available from CRAN
if (!require("spacexr", character.only = TRUE, quietly = TRUE)) {
  cat("Installing spacexr from GitHub...\n")
  if (!require("devtools", character.only = TRUE, quietly = TRUE)) {
    install.packages("devtools")
  }
  devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
}

cat("\n=== Installation Summary ===\n")
cat("Checking all required packages...\n")

all_installed <- TRUE
for (pkg in required_packages) {
  if (require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("✓", pkg, "\n")
  } else {
    cat("✗", pkg, "- FAILED\n")
    all_installed <- FALSE
  }
}

if (all_installed) {
  cat("\n✓ All packages installed successfully!\n")
  cat("You can now run the spatial annotation workflow.\n")
} else {
  cat("\n✗ Some packages failed to install.\n")
  cat("Please check the error messages above and install missing packages manually.\n")
}

cat("\nFor more information, see the README.md file.\n")

