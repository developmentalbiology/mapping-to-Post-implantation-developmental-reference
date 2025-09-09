#!/usr/bin/env Rscript
# Run RCTD Analysis with Human Reference Data
#
# This script demonstrates how to run RCTD analysis using the provided
# human reference dataset and your spatial transcriptomics data.
#
# Prerequisites:
# - Download human_ref.rds from the provided Google Drive link
# - Have your spatial data in RDS format
#
# Author: Spatial Annotation Workflow
# Version: 2.0

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(spacexr)
  library(Matrix)
  library(doParallel)
  library(parallel)
  library(quadprog)
  library(argparse)
})

#' Run RCTD analysis with human reference data
#' 
#' @param human_ref_path Path to human reference RDS file
#' @param spatial_path Path to spatial data RDS file
#' @param output_path Path for output RCTD RDS file
#' @param cell_type_column Column name for cell types in reference (default: "cell_type")
#' @param spatial_coords_columns Vector of column names for spatial coordinates
#' @param max_cores Maximum number of cores to use
#' @param doublet_mode RCTD doublet mode
run_rctd_with_human_ref <- function(human_ref_path,
                                   spatial_path,
                                   output_path,
                                   cell_type_column = "cell_type",
                                   spatial_coords_columns = c("x", "y"),
                                   max_cores = 20,
                                   doublet_mode = "full") {
  
  cat("=== RCTD Analysis with Human Reference ===\n")
  cat("Human reference file:", human_ref_path, "\n")
  cat("Spatial data file:", spatial_path, "\n")
  cat("Output file:", output_path, "\n")
  
  # Load human reference data
  cat("Loading human reference data...\n")
  if (!file.exists(human_ref_path)) {
    stop("Human reference file does not exist: ", human_ref_path, "\n",
         "Please download from: https://drive.google.com/file/d/1lhVc_tZCgXcWecW8aV7iINdDz696_r7r/view?usp=sharing")
  }
  
  human_ref <- readRDS(human_ref_path)
  human_ref <- UpdateSeuratObject(human_ref)
  
  # Check available cell type columns
  cat("Available metadata columns in reference:\n")
  print(colnames(human_ref@meta.data))
  
  # Check if specified cell type column exists
  if (!cell_type_column %in% colnames(human_ref@meta.data)) {
    available_cols <- colnames(human_ref@meta.data)
    potential_cols <- available_cols[grepl("type|anno|cluster|label", available_cols, ignore.case = TRUE)]
    
    if (length(potential_cols) > 0) {
      cat("Cell type column '", cell_type_column, "' not found.\n")
      cat("Potential cell type columns found:\n")
      print(potential_cols)
      cat("Please specify the correct column name using --cell-type-column\n")
    }
    stop("Cell type column '", cell_type_column, "' not found in reference metadata")
  }
  
  # Set cell type identity
  human_ref[[cell_type_column]] <- as.character(human_ref[[cell_type_column]])
  Idents(human_ref) <- cell_type_column
  
  cat("Human reference dataset loaded:\n")
  cat("  Cells:", ncol(human_ref), "\n")
  cat("  Genes:", nrow(human_ref), "\n")
  cat("  Cell types:", length(unique(Idents(human_ref))), "\n")
  
  # Remove clusters with only one cell
  cell_counts <- table(Idents(human_ref))
  valid_clusters <- names(cell_counts[cell_counts > 1])
  
  if (length(valid_clusters) < length(cell_counts)) {
    removed_clusters <- names(cell_counts[cell_counts <= 1])
    cat("Removing clusters with ≤1 cells:", paste(removed_clusters, collapse = ", "), "\n")
    human_ref <- subset(human_ref, idents = valid_clusters)
  }
  
  cat("Final reference cell types:", length(unique(Idents(human_ref))), "\n")
  
  # Extract reference information for RCTD
  ref_counts <- human_ref[["RNA"]]$counts
  ref_cluster <- as.factor(Idents(human_ref))
  names(ref_cluster) <- colnames(human_ref)
  ref_nUMI <- human_ref$nCount_RNA
  names(ref_nUMI) <- colnames(human_ref)
  
  # Create RCTD Reference object
  cat("Creating RCTD Reference object...\n")
  reference <- Reference(ref_counts, ref_cluster, ref_nUMI)
  
  # Load spatial data
  cat("Loading spatial data...\n")
  if (!file.exists(spatial_path)) {
    stop("Spatial data file does not exist: ", spatial_path)
  }
  
  spatial_data <- readRDS(spatial_path)
  spatial_data <- UpdateSeuratObject(spatial_data)
  
  cat("Spatial dataset loaded:\n")
  cat("  Spots:", ncol(spatial_data), "\n")
  cat("  Genes:", nrow(spatial_data), "\n")
  
  # Check available coordinate columns
  cat("Available metadata columns in spatial data:\n")
  print(colnames(spatial_data@meta.data))
  
  # Extract spatial coordinates
  if (!all(spatial_coords_columns %in% colnames(spatial_data@meta.data))) {
    missing_cols <- spatial_coords_columns[!spatial_coords_columns %in% colnames(spatial_data@meta.data)]
    available_cols <- colnames(spatial_data@meta.data)
    potential_cols <- available_cols[grepl("x|y|coord", available_cols, ignore.case = TRUE)]
    
    cat("Spatial coordinate columns not found:", paste(missing_cols, collapse = ", "), "\n")
    if (length(potential_cols) > 0) {
      cat("Potential coordinate columns found:\n")
      print(potential_cols)
      cat("Please specify the correct column names using --spatial-coords\n")
    }
    stop("Spatial coordinate columns not found: ", paste(missing_cols, collapse = ", "))
  }
  
  coords <- spatial_data@meta.data[, spatial_coords_columns]
  coords <- as.data.frame(coords)
  colnames(coords) <- c("x", "y")
  
  cat("Spatial coordinates extracted from columns:", paste(spatial_coords_columns, collapse = ", "), "\n")
  
  # Extract spatial counts
  spatial_counts <- spatial_data[["RNA"]]$counts
  
  # Create RCTD SpatialRNA object
  cat("Creating RCTD SpatialRNA object...\n")
  query <- SpatialRNA(coords, spatial_counts, colSums(spatial_counts))
  
  # Create RCTD object
  cat("Creating RCTD object...\n")
  cat("Parameters:\n")
  cat("  Max cores:", max_cores, "\n")
  cat("  Doublet mode:", doublet_mode, "\n")
  
  RCTD <- create.RCTD(query, reference, max_cores = max_cores)
  
  # Run RCTD analysis
  cat("Running RCTD analysis...\n")
  start_time <- Sys.time()
  
  RCTD <- run.RCTD(RCTD, doublet_mode = doublet_mode)
  
  end_time <- Sys.time()
  elapsed_time <- difftime(end_time, start_time, units = "mins")
  cat("RCTD analysis completed in", round(elapsed_time, 2), "minutes\n")
  
  # Save results
  cat("Saving RCTD results...\n")
  output_dir <- dirname(output_path)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  saveRDS(RCTD, file = output_path)
  cat("RCTD results saved to:", output_path, "\n")
  
  # Print summary
  cat("\n=== RCTD Analysis Summary ===\n")
  cat("Reference cell types:", length(unique(reference@cell_types)), "\n")
  cat("Spatial spots analyzed:", nrow(RCTD@spatialRNA@coords), "\n")
  cat("Doublet mode:", doublet_mode, "\n")
  cat("Processing time:", round(elapsed_time, 2), "minutes\n")
  
  # Print cell type summary
  cat("\nReference cell types:\n")
  ref_summary <- table(reference@cell_types)
  for (i in 1:length(ref_summary)) {
    cat("  ", names(ref_summary)[i], ":", ref_summary[i], "cells\n")
  }
  
  return(RCTD)
}

# Main execution
main <- function() {
  # Parse command line arguments
  parser <- ArgumentParser(description = 'Run RCTD analysis with human reference data')
  parser$add_argument('--human-ref', '-r', required = TRUE,
                     help = 'Path to human reference RDS file')
  parser$add_argument('--spatial', '-s', required = TRUE,
                     help = 'Path to spatial RDS file')
  parser$add_argument('--output', '-o', required = TRUE,
                     help = 'Output path for RCTD results RDS file')
  parser$add_argument('--cell-type-column', default = 'cell_type',
                     help = 'Column name for cell types in reference (default: cell_type)')
  parser$add_argument('--spatial-coords', nargs = 2, default = c('x', 'y'),
                     help = 'Column names for spatial coordinates (default: x y)')
  parser$add_argument('--max-cores', type = 'integer', default = 20,
                     help = 'Maximum number of cores to use (default: 20)')
  parser$add_argument('--doublet-mode', default = 'full',
                     choices = c('full', 'doublet', 'single'),
                     help = 'RCTD doublet mode (default: full)')
  
  args <- parser$parse_args()
  
  # Run RCTD analysis
  run_rctd_with_human_ref(
    human_ref_path = args$human_ref,
    spatial_path = args$spatial,
    output_path = args$output,
    cell_type_column = args$cell_type_column,
    spatial_coords_columns = args$spatial_coords,
    max_cores = args$max_cores,
    doublet_mode = args$doublet_mode
  )
}

# Run main function if script is executed directly
if (!interactive()) {
  main()
}

