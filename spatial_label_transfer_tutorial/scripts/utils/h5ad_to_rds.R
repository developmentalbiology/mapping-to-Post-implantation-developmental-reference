#!/usr/bin/env Rscript
# Convert H5AD Export to Seurat RDS Object
#
# This script converts the exported h5ad data to a Seurat object and saves as RDS
# 
# Required inputs (from export_h5ad.py):
# - matrix.mtx: Expression matrix
# - barcodes.tsv: Cell barcodes
# - features.tsv: Gene features  
# - obs_metadata.csv: Cell metadata
# - umap_coordinates.csv: UMAP coordinates (optional)
#
# Author: Spatial Annotation Workflow
# Version: 2.0

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(readr)
  library(data.table)
  library(argparse)
})

#' Convert exported h5ad data to Seurat RDS object
#' 
#' @param input_dir Directory containing exported h5ad data
#' @param output_file Path for output RDS file
#' @param project_name Name for the Seurat project
convert_to_seurat <- function(input_dir, output_file, project_name = "spatial_data") {
  
  cat("Converting h5ad export to Seurat object...\n")
  cat("Input directory:", input_dir, "\n")
  cat("Output file:", output_file, "\n")
  
  # Check required files
  required_files <- c("matrix.mtx", "barcodes.tsv", "features.tsv", "obs_metadata.csv")
  missing_files <- c()
  
  for (file in required_files) {
    file_path <- file.path(input_dir, file)
    if (!file.exists(file_path)) {
      missing_files <- c(missing_files, file)
    }
  }
  
  if (length(missing_files) > 0) {
    stop("Missing required files: ", paste(missing_files, collapse = ", "))
  }
  
  # Read expression matrix
  cat("Reading expression matrix...\n")
  mtx_path <- file.path(input_dir, "matrix.mtx")
  mtx <- readMM(mtx_path)
  
  # Read barcodes and features
  cat("Reading barcodes and features...\n")
  bc_path <- file.path(input_dir, "barcodes.tsv")
  fs_path <- file.path(input_dir, "features.tsv")
  
  bc <- read_tsv(bc_path, col_names = FALSE, show_col_types = FALSE)
  fs <- read_tsv(fs_path, col_names = FALSE, show_col_types = FALSE)
  
  # Check matrix dimensions
  cat("Matrix dimensions:\n")
  cat("  Cells (columns):", ncol(mtx), "\n")
  cat("  Genes (rows):", nrow(mtx), "\n")
  cat("  Barcodes:", nrow(bc), "\n")
  cat("  Features:", nrow(fs), "\n")
  
  # Assign row and column names
  colnames(mtx) <- bc$X1
  rownames(mtx) <- fs$X1
  
  # Read metadata
  cat("Reading metadata...\n")
  meta_path <- file.path(input_dir, "obs_metadata.csv")
  meta <- read.csv(meta_path, row.names = 1)
  
  # Create Seurat object
  cat("Creating Seurat object...\n")
  seurat_obj <- CreateSeuratObject(
    counts = mtx, 
    project = project_name,
    meta.data = meta
  )
  
  cat("Seurat object created with", ncol(seurat_obj), "cells and", nrow(seurat_obj), "genes\n")
  
  # Add UMAP coordinates if available
  umap_path <- file.path(input_dir, "umap_coordinates.csv")
  if (file.exists(umap_path)) {
    cat("Adding UMAP coordinates...\n")
    umap_coords <- read.csv(umap_path, row.names = 1)
    
    # Check cell name compatibility
    seurat_cells <- colnames(seurat_obj)
    umap_cells <- rownames(umap_coords)
    
    # Find missing cells
    missing_in_umap <- setdiff(seurat_cells, umap_cells)
    missing_in_seurat <- setdiff(umap_cells, seurat_cells)
    
    if (length(missing_in_umap) > 0) {
      cat("Warning: Cells in Seurat object but missing in UMAP:", length(missing_in_umap), "\n")
    }
    
    if (length(missing_in_seurat) > 0) {
      cat("Warning: Cells in UMAP but missing in Seurat object:", length(missing_in_seurat), "\n")
    }
    
    # Add UMAP coordinates to Seurat object
    if (nrow(umap_coords) > 0) {
      seurat_obj[['umap']] <- CreateDimReducObject(
        embeddings = as.matrix(umap_coords), 
        key = "UMAP_", 
        assay = DefaultAssay(seurat_obj)
      )
      cat("UMAP coordinates added successfully\n")
    }
  } else {
    cat("Warning: UMAP coordinates file not found, skipping...\n")
  }
  
  # Save RDS file
  cat("Saving Seurat object to RDS...\n")
  saveRDS(seurat_obj, file = output_file)
  cat("Seurat object saved to:", output_file, "\n")
  
  # Print summary
  cat("\n=== Conversion Summary ===\n")
  cat("Project name:", project_name, "\n")
  cat("Cells:", ncol(seurat_obj), "\n")
  cat("Genes:", nrow(seurat_obj), "\n")
  cat("Metadata columns:", ncol(seurat_obj@meta.data), "\n")
  cat("Reductions:", length(seurat_obj@reductions), "\n")
  
  return(seurat_obj)
}

# Main execution
main <- function() {
  # Parse command line arguments
  parser <- ArgumentParser(description = 'Convert h5ad export to Seurat RDS object')
  parser$add_argument('--input', '-i', required = TRUE,
                     help = 'Input directory containing exported h5ad data')
  parser$add_argument('--output', '-o', required = TRUE,
                     help = 'Output RDS file path')
  parser$add_argument('--project', '-p', default = 'spatial_data',
                     help = 'Project name for Seurat object (default: spatial_data)')
  
  args <- parser$parse_args()
  
  # Validate input directory
  if (!dir.exists(args$input)) {
    stop("Input directory does not exist: ", args$input)
  }
  
  # Create output directory if needed
  output_dir <- dirname(args$output)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Convert to Seurat
  convert_to_seurat(args$input, args$output, args$project)
}

# Run main function if script is executed directly
if (!interactive()) {
  main()
}

