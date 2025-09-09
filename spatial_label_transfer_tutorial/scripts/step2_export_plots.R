#!/usr/bin/env Rscript
# Step 2: Export Weight Plots and Visualizations (Streamlined Version)
#
# This script generates essential visualizations of the RCTD weight scores
# with expensive operations removed for faster processing
#
# Required inputs:
# - RCTD results RDS file
# - Original spatial data RDS file
# - Configuration parameters
#
# Author: Spatial Annotation Workflow
# Version: 2.1 (Streamlined)

# Set library paths (add this at the very beginning)
.libPaths(c("/home/liuxiaodongLab/fanxueying/R_LIBS",
            "/home/liuxiaodongLab/fanxueying/miniconda3/envs/r_env/lib/R/library",
            "/soft/devtools/R/R-4.4.3_installation/lib64/R/library"))

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(spacexr)
  library(Matrix)
  library(ggplot2)
  library(RColorBrewer)
  library(viridis)
  library(patchwork)
  library(yaml)
})

#' Generate streamlined visualization plots for RCTD weight analysis
#' 
#' @param rctd_path Path to RCTD results RDS file
#' @param spatial_path Path to original spatial data RDS file
#' @param output_dir Directory to save output plots
#' @param config_path Path to configuration YAML file
#' @param plot_width Width of plots in inches
#' @param plot_height Height of plots in inches
#' @param plot_dpi DPI for plot resolution
generate_weight_plots <- function(rctd_path,
                                 spatial_path,
                                 output_dir,
                                 config_path = NULL,
                                 plot_width = 12,
                                 plot_height = 8,
                                 plot_dpi = 300) {
  
  cat("=== Generating RCTD Weight Analysis Plots (Streamlined) ===\n")
  cat("RCTD results file:", rctd_path, "\n")
  cat("Spatial data file:", spatial_path, "\n")
  cat("Output directory:", output_dir, "\n")
  
  # Load configuration if provided
  if (!is.null(config_path) && file.exists(config_path)) {
    config <- yaml::read_yaml(config_path)
    params <- config$plots
    cat("Configuration loaded from:", config_path, "\n")
  } else {
    # Default parameters
    params <- list(
      point_size = 0.5,
      plot_title_size = 14,
      axis_text_size = 10,
      legend_text_size = 8
    )
    cat("Using default configuration parameters\n")
  }
  
  # Create output directories
  plots_dir <- file.path(output_dir, "plots")
  if (!dir.exists(plots_dir)) {
    dir.create(plots_dir, recursive = TRUE)
  }
  
  # Load RCTD results
  cat("Loading RCTD results...\n")
  rctd_obj <- readRDS(rctd_path)
  weights <- rctd_obj@results$weights
  
  # Load spatial data
  cat("Loading spatial data...\n")
  spatial_data <- readRDS(spatial_path)
  
  # Get spatial coordinates
  coords <- spatial_data@meta.data[, c("x", "y")]
  
  # Match spots between weights and coordinates
  common_spots <- intersect(rownames(weights), rownames(coords))
  weights <- weights[common_spots, ]
  coords <- coords[common_spots, ]
  
  cat("Data summary:\n")
  cat("  Spots:", nrow(weights), "\n")
  cat("  Cell types:", ncol(weights), "\n")
  cat("  Matching spots:", length(common_spots), "\n")
  
  # Get cell type names
  cell_type_names <- colnames(weights)
  
  # 1. Generate basic weight summary statistics
  cat("Generating weight summary statistics...\n")
  
  weight_stats <- data.frame(
    cell_type = cell_type_names,
    mean_weight = colMeans(weights),
    median_weight = apply(weights, 2, median),
    max_weight = apply(weights, 2, max),
    min_weight = apply(weights, 2, min),
    spots_with_weight = colSums(weights > 0),
    spots_high_weight = colSums(weights > 0.5)
  )
  
  write.csv(weight_stats, file.path(output_dir, "weight_statistics.csv"), row.names = FALSE)
  
  # 2. Generate individual cell type weight plots (exactly matching original script style)
  cat("Generating individual cell type weight plots...\n")
  
  # Simplified output path - save directly to output_dir
  pdf_file <- file.path(output_dir, "all_cell_type_weights.pdf")
  
  # Custom plot function matching your original script exactly (including grid removal)
  plot_puck_wrapper <- function(puck_coords, plot_val, title = NULL) {
    plot_data <- data.frame(
      x = puck_coords[, 1],
      y = puck_coords[, 2], 
      weight = plot_val
    )
    
    p <- ggplot(plot_data, aes(x = x, y = y, color = weight)) +
      geom_point(size = 0.1, alpha = 0.8) +  # Small dots like original
      scale_color_gradientn(colors = c("blue", "yellow", "red"), limits = c(0, 1)) +  # Exact original colors
      labs(title = title, x = "X", y = "Y") +
      theme_minimal() +
      theme(
        panel.background = element_rect(fill = "black"),
        plot.background = element_rect(fill = "black"),
        plot.title = element_text(color = "white", size = 36),  # Same large title size as original
        # Remove grid lines and axis elements like in original script
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right"  # Keep legend visible
      ) +
      coord_fixed()
    
    return(p)
  }
  
  # Normalize weights like in original script
  norm_weights <- weights / rowSums(weights + 1e-10)  # Add small value to avoid division by zero
  
  # Create plots list like in original script
  plots <- vector(mode = "list", length = length(cell_type_names))
  
  for (i in 1:length(cell_type_names)) {
    cell_type <- cell_type_names[i]
    cat("  Processing:", cell_type, "...\n")
    
    plot_var <- norm_weights[, cell_type]
    names(plot_var) <- rownames(norm_weights)
    
    if (sum(norm_weights[, cell_type]) > 0) {
      plots[[i]] <- plot_puck_wrapper(coords, plot_var, title = cell_type)
    }
  }
  
  # Save all plots in single PDF exactly like original script
  pdf(pdf_file, width = 25, height = 10)  # Same dimensions as original
  invisible(lapply(plots, print))
  dev.off()
  
  cat("Generated single PDF with all", ncol(weights), "cell type weight plots:", pdf_file, "\n")
  
  # Skip expensive operations (spatial overview, detailed histograms, etc.)
  cat("Skipping computationally expensive plots (histograms, overview) for faster processing...\n")
  
  cat("\n=== Streamlined Weight Analysis Complete ===\n")
  cat("Output directory:", output_dir, "\n")
  cat("Generated files:\n")
  cat("  - all_cell_type_weights.pdf (all cell types in single PDF, clean format)\n")
  cat("  - weight_statistics.csv (summary statistics)\n")
  
}

# Simple argument parsing without argparse
parse_simple_args <- function(args) {
  result <- list()
  
  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--rctd-results" && i < length(args)) {
      result$rctd <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--spatial-data" && i < length(args)) {
      result$spatial <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--output-dir" && i < length(args)) {
      result$output <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--config" && i < length(args)) {
      result$config <- args[i + 1]
      i <- i + 2
    } else {
      i <- i + 1
    }
  }
  
  return(result)
}

# Main execution
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 6) {
    cat("Error: Missing required arguments\n")
    cat("Required: --rctd-results, --spatial-data, --output-dir\n")
    quit(status = 1)
  }
  
  parsed_args <- parse_simple_args(args)
  
  if (is.null(parsed_args$rctd) || is.null(parsed_args$spatial) || is.null(parsed_args$output)) {
    cat("Error: Missing required arguments\n")
    cat("Required: --rctd-results, --spatial-data, --output-dir\n")
    quit(status = 1)
  }
  
  # Run the analysis
  generate_weight_plots(
    rctd_path = parsed_args$rctd,
    spatial_path = parsed_args$spatial,
    output_dir = parsed_args$output,
    config_path = parsed_args$config
  )
}

# Run main function
if (!interactive()) {
  main()
}