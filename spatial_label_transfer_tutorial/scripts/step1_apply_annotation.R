#!/usr/bin/env Rscript
# Step 1: Apply Multi-Tier Assignment for Final Annotation
#
# This script applies the multi-tier assignment algorithm to RCTD results
# to generate final cell type annotations with enhanced accuracy
#
# Required inputs:
# - RCTD results RDS file (generated using spacexr package)
# - Spatial data RDS file (original spatial data used for RCTD)
# - Configuration parameters
#
# Author: Spatial Annotation Workflow
# Version: 2.0

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(spacexr)
  library(Matrix)
  library(argparse)
  library(yaml)
})

# Source the multi-tier assignment function
# Get script directory - robust method that works with both source() and Rscript
get_script_dir <- function() {
  # Try multiple methods to get script directory
  script_dir <- tryCatch({
    if (sys.nframe() == 0) {
      # Running with Rscript - try different approaches
      args <- commandArgs(trailingOnly = FALSE)
      file_arg <- grep("--file=", args, value = TRUE)
      if (length(file_arg) > 0) {
        script_path <- sub("--file=", "", file_arg)
        dirname(normalizePath(script_path))
      } else {
        # Alternative method
        dirname(normalizePath(sys.frame(1)$ofile))
      }
    } else {
      # Being sourced
      dirname(sys.frame(1)$ofile)
    }
  }, error = function(e) {
    # If all methods fail, assume we're in the scripts directory or main directory
    current_dir <- getwd()
    if (file.exists(file.path(current_dir, "multi_tier_integrated_final.R"))) {
      return(current_dir)
    } else if (file.exists(file.path(current_dir, "scripts", "multi_tier_integrated_final.R"))) {
      return(file.path(current_dir, "scripts"))
    } else {
      return(current_dir)
    }
  })
  
  # Ensure we have a valid directory and the target file exists
  if (is.na(script_dir) || is.null(script_dir) || script_dir == "") {
    current_dir <- getwd()
    if (file.exists(file.path(current_dir, "scripts", "multi_tier_integrated_final.R"))) {
      script_dir <- file.path(current_dir, "scripts")
    } else {
      script_dir <- current_dir
    }
  }
  
  return(script_dir)
}

script_dir <- get_script_dir()

# Check if the required file exists, if not try alternative paths
multi_tier_file <- file.path(script_dir, "multi_tier_integrated_final.R")
if (!file.exists(multi_tier_file)) {
  # Try in scripts subdirectory
  alt_path <- file.path(getwd(), "scripts", "multi_tier_integrated_final.R")
  if (file.exists(alt_path)) {
    multi_tier_file <- alt_path
  } else {
    stop("Cannot find multi_tier_integrated_final.R file. Please ensure it exists in the scripts directory.")
  }
}

source(multi_tier_file)

#' Apply multi-tier assignment to RCTD results
#' 
#' @param rctd_path Path to RCTD results RDS file
#' @param spatial_path Path to spatial data RDS file
#' @param output_path Path for output annotated spatial data
#' @param config_path Path to configuration YAML file
#' @param spatial_coords_columns Vector of column names for spatial coordinates
apply_multi_tier_annotation <- function(rctd_path,
                                       spatial_path,
                                       output_path,
                                       config_path = NULL,
                                       spatial_coords_columns = c("x", "y")) {
  
  cat("=== Multi-Tier Assignment for Spatial Annotation ===\n")
  cat("RCTD results file:", rctd_path, "\n")
  cat("Spatial data file:", spatial_path, "\n")
  cat("Output file:", output_path, "\n")
  
  # Default parameters for multi-tier assignment
  params <- list(
    # Basic parameters
    excluded_cell_types = NULL,
    protected_cell_types = NULL,
    max_proportion = 0.6,
    min_proportion = 0.01,
    proportion_correction = TRUE,
    correction_method = "iterative",
    k_neighbors = 5,
    spatial_weight = 0.6,
    protection_threshold = 0.2,
    dominance_threshold = 1.2,
    enable_tier0_protection = TRUE,
    
    # Signal balancing parameters
    suppressed_cell_types = NULL,
    max_suppressed_proportion = 0.25,
    suppression_strength = 0.8,
    weak_signal_protected_types = NULL,
    weak_signal_threshold = 0.08,
    weak_signal_dominance_threshold = 1.05,
    weak_signal_spatial_boost = 1.2,
    enable_signal_balancing = TRUE,
    
    # Performance parameters
    enable_optimizations = TRUE,
    verbose = TRUE
  )
  
  # Load configuration if provided
  if (!is.null(config_path) && file.exists(config_path)) {
    cat("Loading configuration from:", config_path, "\n")
    config <- yaml::read_yaml(config_path)
    
    # Update parameters with config values
    if ("annotation" %in% names(config)) {
      annotation_config <- config$annotation
      for (param_name in names(annotation_config)) {
        if (param_name %in% names(params)) {
          params[[param_name]] <- annotation_config[[param_name]]
        }
      }
    }
  }
  
  # Load RCTD results
  cat("Loading RCTD results...\n")
  if (!file.exists(rctd_path)) {
    stop("RCTD results file does not exist: ", rctd_path)
  }
  RCTD <- readRDS(rctd_path)
  
  # Load spatial data
  cat("Loading spatial data...\n")
  if (!file.exists(spatial_path)) {
    stop("Spatial data file does not exist: ", spatial_path)
  }
  spatial_data <- readRDS(spatial_path)
  
  # Extract RCTD results
  cat("Extracting RCTD weights...\n")
  results <- RCTD@results
  norm_weights <- normalize_weights(results$weights)
  cell_type_names <- RCTD@cell_type_info$info[[2]]
  spatialRNA <- RCTD@spatialRNA
  
  # Get spatial coordinates
  if (!all(spatial_coords_columns %in% colnames(spatial_data@meta.data))) {
    missing_cols <- spatial_coords_columns[!spatial_coords_columns %in% colnames(spatial_data@meta.data)]
    stop("Spatial coordinate columns not found: ", paste(missing_cols, collapse = ", "))
  }
  
  xy <- spatial_data@meta.data[, spatial_coords_columns]
  xy <- as.matrix(xy)
  
  # Ensure coordinates match RCTD spots
  xy <- xy[row.names(spatialRNA@coords), ]
  spatialRNA@coords <- as.data.frame(xy)
  
  # Add spatial reduction to Seurat object
  colnames(xy) <- c("spatial_1", "spatial_2")
  spatial_reduction <- CreateDimReducObject(
    embeddings = xy, 
    key = "spatial_", 
    assay = DefaultAssay(spatial_data)
  )
  spatial_data@reductions[["spatial"]] <- spatial_reduction
  
  # Convert weights matrix
  norm_weights <- as.matrix(norm_weights)
  
  cat("RCTD results summary:\n")
  cat("  Spots:", nrow(norm_weights), "\n")
  cat("  Cell types:", ncol(norm_weights), "\n")
  cat("  Cell type names:", paste(head(cell_type_names, 5), collapse = ", "), 
      if(length(cell_type_names) > 5) "..." else "", "\n")
  
  # Print parameter summary
  cat("\n=== Multi-Tier Assignment Parameters ===\n")
  if (!is.null(params$excluded_cell_types)) {
    cat("Excluded cell types:", paste(params$excluded_cell_types, collapse = ", "), "\n")
  }
  if (!is.null(params$protected_cell_types)) {
    cat("Protected cell types:", paste(params$protected_cell_types, collapse = ", "), "\n")
  }
  cat("Max proportion:", params$max_proportion, "\n")
  cat("Spatial weight:", params$spatial_weight, "\n")
  cat("K neighbors:", params$k_neighbors, "\n")
  
  if (params$enable_signal_balancing) {
    cat("Signal balancing: ENABLED\n")
    if (!is.null(params$suppressed_cell_types)) {
      cat("Suppressed cell types:", paste(params$suppressed_cell_types, collapse = ", "), "\n")
    }
    if (!is.null(params$weak_signal_protected_types)) {
      cat("Weak signal protected:", paste(params$weak_signal_protected_types, collapse = ", "), "\n")
    }
  }
  
  # Apply multi-tier assignment
  cat("\nApplying multi-tier assignment algorithm...\n")
  
  assignment_results <- do.call(multi_tier_assignment_integrated_final, c(
    list(
      weights_matrix = norm_weights,
      cell_type_names = cell_type_names,
      spatial_coords = spatialRNA@coords
    ),
    params
  ))
  
  # Add results to spatial data (using coordinate-based matching like original script)
  cat("Adding annotation results to spatial data...\n")
  
  # Direct assignment without spot name matching (coordinates already matched)
  # This follows the same approach as your original script
  spatial_data$final_annotation <- assignment_results$final_assignment
  spatial_data$annotation_confidence <- assignment_results$assignment_confidence
  spatial_data$assignment_tier <- assignment_results$assignment_tier
  spatial_data$tier0_protected <- assignment_results$tier0_protected
  spatial_data$signal_balancing_applied <- assignment_results$signal_balancing_applied
  spatial_data$weak_signal_protected <- assignment_results$weak_signal_protected
  
  # Add weight information (direct assignment like annotation columns)
  for (i in 1:length(cell_type_names)) {
    col_name <- paste0(cell_type_names[i], "_weight")
    if (col_name %in% colnames(assignment_results)) {
      spatial_data@meta.data[[col_name]] <- assignment_results[[col_name]]
    }
  }
  
  # Save annotated spatial data
  cat("Saving annotated spatial data...\n")
  output_dir <- dirname(output_path)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  saveRDS(spatial_data, file = output_path)
  cat("Annotated spatial data saved to:", output_path, "\n")
  
  # Generate final assignment UMAP plot
  cat("Generating final assignment UMAP plot...\n")
  
  # Define color palette (same as your previous script)
  godsnot_102 <- c(
    "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
    "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF",
    "#997D87", "#5A0007", "#809693", "#6A3A4C", "#1B4400", "#4FC601", "#3B5DFF",
    "#4A3B53", "#FF2F80", "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92",
    "#FF90C9", "#B903AA", "#D16100", "#DDEFFF", "#000035", "#7B4F4B", "#A1C299",
    "#300018", "#0AA6D8", "#013349", "#00846F", "#372101", "#FFB500", "#C2FFED",
    "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09", "#00489C", "#6F0062",
    "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66", "#885578",
    "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
    "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F",
    "#938A81", "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757",
    "#C8A1A1", "#1E6E00", "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C",
    "#772600", "#D790FF", "#9B9700", "#549E79", "#FFF69F", "#201625", "#72418F",
    "#BC23FF", "#99ADC0", "#3A2465", "#922329", "#5B4534", "#FDE8DC", "#404E55",
    "#0089A3", "#CB7E98", "#A4E804", "#324E72"
  )
  
  # Define ordered labels (same as your previous script)
  final_anno_labels <- c('TE', 'CTB_1','CTB_2', 'STB_1', 'STB_2', 'STB_3', 'EVT_1', 'EVT_2',
                        'Epiblast_1','Epiblast_2','Epiblast_3','Ectoderm',
                        'Amniontic.epi','Amniontic.ectoderm',
                        'PGC',
                        'Primitive.streak',
                        'Neuromesodermal.progenitor',
                        'Neural.crest', 'Neural.ectoderm.forebrain', 'Neural.ectoderm.hindbrain', 'Neural.ectoderm.midbrain','Spinal.cord',
                        'Paraxial.mesoderm','Emergent.mesoderm','Pre-somatic.mesoderm','Somite', 'Rostral.mesoderm', 'Lateral.plate.mesoderm_1',
                        'Lateral.plate.mesoderm_2','Lateral.plate.mesoderm_3','Cardiac.mesoderm','Amniotic.mesoderm','Exe.meso.progenitor','YS.mesoderm_1', 'YS.mesoderm_2',
                        'Hypoblast_1', 'Hypoblast_2', 'AVE', 'VE', 'YS.endoderm',
                        'DE','Gut',
                        'Notochord',
                        'Hemogenic.endothelial.progenitor','Endothelium','Erythroid','Primitive.megakaryocyte','Myeloid.progenitor')
  
  # Reorder and factor the final annotation column (check actual column name)
  anno_column <- if ("final_anno" %in% colnames(spatial_data@meta.data)) {
    "final_anno"
  } else if ("final_annotation" %in% colnames(spatial_data@meta.data)) {
    "final_annotation"
  } else {
    # Find any column with "final" and "anno" in the name
    possible_cols <- grep("final.*anno|anno.*final", colnames(spatial_data@meta.data), value = TRUE, ignore.case = TRUE)
    if (length(possible_cols) > 0) {
      possible_cols[1]
    } else {
      stop("Cannot find final annotation column. Available columns: ", paste(colnames(spatial_data@meta.data), collapse = ", "))
    }
  }
  
  cat("Using annotation column:", anno_column, "\n")
  
  # Factor the annotation column (same as your original script)
  spatial_data@meta.data[[anno_column]] <- factor(spatial_data@meta.data[[anno_column]], levels = final_anno_labels, ordered = TRUE)
  spatial_data@meta.data[[anno_column]] <- droplevels(spatial_data@meta.data[[anno_column]])
  
  # Create color palette for ALL final_anno_labels (not just existing ones) to maintain consistent colors
  # This matches your original script approach - lines 128 in your original
  color_palette <- setNames(godsnot_102[1:length(final_anno_labels)], final_anno_labels)
  
  # Generate the UMAP plot with small points like in original script
  p <- DimPlot(spatial_data, reduction = "spatial", group.by = anno_column) +
    scale_color_manual(values = color_palette) +
    ggtitle("Final Assignment - Spatial Plot")
  
  # Save the plot in the same directory as the output with correct naming
  output_dir <- dirname(output_path)
  output_base <- tools::file_path_sans_ext(basename(output_path))
  plot_path <- file.path(output_dir, paste0(output_base, "_final_assignment.pdf"))
  
  cat("Saving plot to:", plot_path, "\n")
  pdf(plot_path, width = 18, height = 15)
  print(p)
  dev.off()
  
  cat("Final assignment plot saved to:", plot_path, "\n")
  
  # Save detailed assignment results
  results_csv_path <- gsub("\\.rds$", "_assignment_results.csv", output_path)
  write.csv(assignment_results, results_csv_path, row.names = TRUE)
  cat("Detailed assignment results saved to:", results_csv_path, "\n")
  
  # Print summary
  cat("\n=== Final Annotation Summary ===\n")
  final_counts <- table(spatial_data$final_annotation)
  cat("Final cell type distribution:\n")
  for (i in 1:length(final_counts)) {
    cell_type <- names(final_counts)[i]
    count <- final_counts[i]
    percentage <- round(count / sum(final_counts) * 100, 1)
    cat("  ", cell_type, ":", count, "(", percentage, "%)\n")
  }
  
  # Tier distribution
  tier_counts <- table(spatial_data$assignment_tier)
  cat("\nTier distribution:\n")
  for (i in 1:length(tier_counts)) {
    tier <- names(tier_counts)[i]
    count <- tier_counts[i]
    percentage <- round(count / sum(tier_counts) * 100, 1)
    cat("  ", tier, ":", count, "(", percentage, "%)\n")
  }
  
  # Protection summary
  if (params$enable_tier0_protection) {
    protected_count <- sum(spatial_data$tier0_protected, na.rm = TRUE)
    cat("\nTier 0 protection:", protected_count, "spots (", 
        round(protected_count / ncol(spatial_data) * 100, 1), "%)\n")
  }
  
  if (params$enable_signal_balancing) {
    signal_balanced_count <- sum(spatial_data$signal_balancing_applied, na.rm = TRUE)
    weak_protected_count <- sum(spatial_data$weak_signal_protected, na.rm = TRUE)
    cat("Signal balancing applied:", signal_balanced_count, "spots\n")
    cat("Weak signal protection:", weak_protected_count, "spots\n")
  }
  
  return(spatial_data)
}

# Main execution
main <- function() {
  # Parse command line arguments
  parser <- ArgumentParser(description = 'Apply multi-tier assignment for final annotation')
  parser$add_argument('--rctd', '-r', required = TRUE,
                     help = 'Path to RCTD results RDS file')
  parser$add_argument('--spatial', '-s', required = TRUE,
                     help = 'Path to spatial data RDS file')
  parser$add_argument('--output', '-o', required = TRUE,
                     help = 'Output path for annotated spatial data RDS file')
  parser$add_argument('--config', '-c',
                     help = 'Path to configuration YAML file')
  parser$add_argument('--coords-columns', nargs = 2, default = c('x', 'y'),
                     help = 'Column names for spatial coordinates (default: x y)')
  
  args <- parser$parse_args()
  
  # Apply multi-tier annotation
  apply_multi_tier_annotation(
    rctd_path = args$rctd,
    spatial_path = args$spatial,
    output_path = args$output,
    config_path = args$config,
    spatial_coords_columns = args$coords_columns
  )
}

# Run main function if script is executed directly
if (!interactive()) {
  main()
}

