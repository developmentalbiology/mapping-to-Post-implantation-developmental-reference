#!/bin/bash

# Spatial Annotation Workflow - Main Execution Script
# 
# This script runs the simplified spatial annotation workflow:
# 1. Apply multi-tier annotation to RCTD results
# 2. Generate comprehensive visualization plots
#
# Prerequisites: RCTD analysis should be completed beforehand
# For RCTD analysis, please refer to: https://github.com/dmcable/spacexr
#
# Author: Spatial Annotation Workflow
# Version: 2.0

set -e  # Exit on any error

# Default parameters
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_DIR="$(dirname "$SCRIPT_DIR")"
CONFIG_FILE=""
VERBOSE=true

# Function to print usage
usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Spatial Annotation Workflow - Post-RCTD Analysis Pipeline

Required Arguments:
  --rctd-results PATH       Path to RCTD results RDS file
  --spatial-data PATH       Path to original spatial RDS file
  --output-dir PATH         Output directory for results

Optional Arguments:
  --config PATH             Configuration YAML file
  --spatial-coords X Y      Column names for spatial coordinates (default: x y)
  --skip-step STEP          Skip specific step (1-2, can be used multiple times)
  --verbose                 Enable verbose output (default: true)
  --quiet                   Disable verbose output
  --help                    Show this help message

Steps:
  1. Apply multi-tier assignment algorithm
  2. Generate visualization plots and weight maps

Prerequisites:
  - RCTD analysis completed (see https://github.com/dmcable/spacexr)
  - RCTD results saved as RDS file
  - Original spatial data available

Example:
  $0 --rctd-results rctd_results.rds \\
     --spatial-data spatial_data.rds \\
     --output-dir results/ \\
     --config config.yaml

EOF
}

# Function to log messages
log() {
    if [ "$VERBOSE" = true ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
    fi
}

# Function to check if command exists
check_command() {
    if ! command -v "$1" &> /dev/null; then
        echo "Error: $1 is not installed or not in PATH"
        exit 1
    fi
}

# Parse command line arguments
RCTD_RESULTS=""
SPATIAL_DATA=""
OUTPUT_DIR=""
SPATIAL_COORDS=("x" "y")
SKIP_STEPS=()

while [[ $# -gt 0 ]]; do
    case $1 in
        --rctd-results)
            RCTD_RESULTS="$2"
            shift 2
            ;;
        --spatial-data)
            SPATIAL_DATA="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --config)
            CONFIG_FILE="$2"
            shift 2
            ;;
        --spatial-coords)
            SPATIAL_COORDS=("$2" "$3")
            shift 3
            ;;
        --skip-step)
            SKIP_STEPS+=("$2")
            shift 2
            ;;
        --verbose)
            VERBOSE=true
            shift
            ;;
        --quiet)
            VERBOSE=false
            shift
            ;;
        --help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# Validate required arguments
if [ -z "$RCTD_RESULTS" ] || [ -z "$SPATIAL_DATA" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: Missing required arguments"
    usage
    exit 1
fi

# Check if input files exist
if [ ! -f "$RCTD_RESULTS" ]; then
    echo "Error: RCTD results file does not exist: $RCTD_RESULTS"
    echo "Please run RCTD analysis first. See: https://github.com/dmcable/spacexr"
    exit 1
fi

if [ ! -f "$SPATIAL_DATA" ]; then
    echo "Error: Spatial data file does not exist: $SPATIAL_DATA"
    exit 1
fi

# Check if config file exists (if provided)
if [ -n "$CONFIG_FILE" ] && [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: Configuration file does not exist: $CONFIG_FILE"
    exit 1
fi

# Check required commands
log "Checking required software..."
check_command "Rscript"

# Create output directory
mkdir -p "$OUTPUT_DIR"
log "Created output directory: $OUTPUT_DIR"

# Define output file paths
ANNOTATED_SPATIAL="$OUTPUT_DIR/annotated_spatial.rds"
PLOTS_DIR="$OUTPUT_DIR/plots"

# Function to check if step should be skipped
should_skip_step() {
    local step=$1
    for skip in "${SKIP_STEPS[@]}"; do
        if [ "$skip" = "$step" ]; then
            return 0
        fi
    done
    return 1
}

# Step 1: Apply multi-tier annotation
if ! should_skip_step "1"; then
    log "=== Step 1: Applying multi-tier annotation ==="
    
    CONFIG_ARG=""
    if [ -n "$CONFIG_FILE" ]; then
        CONFIG_ARG="--config $CONFIG_FILE"
    fi
    
    Rscript "$SCRIPT_DIR/step1_apply_annotation.R" \
        --rctd "$RCTD_RESULTS" \
        --spatial "$SPATIAL_DATA" \
        --output "$ANNOTATED_SPATIAL" \
        --coords-columns "${SPATIAL_COORDS[0]}" "${SPATIAL_COORDS[1]}" \
        $CONFIG_ARG
    log "Step 1 completed successfully"
else
    log "Skipping Step 1: Apply multi-tier annotation"
fi

# Step 2: Generate visualization plots
if ! should_skip_step "2"; then
    log "=== Step 2: Generating visualization plots ==="
    
    CONFIG_ARG=""
    if [ -n "$CONFIG_FILE" ]; then
        CONFIG_ARG="--config $CONFIG_FILE"
    fi
    
    Rscript "$SCRIPT_DIR/step2_export_plots.R" \
        --rctd-results "$RCTD_RESULTS" \
        --spatial-data "$SPATIAL_DATA" \
        --output-dir "$PLOTS_DIR" \
        $CONFIG_ARG
    log "Step 2 completed successfully"
else
    log "Skipping Step 2: Generate visualization plots"
fi

log "=== Workflow completed successfully ==="
log "Results saved to: $OUTPUT_DIR"
log ""
log "Output files:"
log "  - Annotated spatial data: $ANNOTATED_SPATIAL"
log "  - Assignment results CSV: ${ANNOTATED_SPATIAL%.*}_assignment_results.csv"
log "  - Visualization plots: $PLOTS_DIR"
log ""
log "Key visualizations:"
log "  - Cell type distribution: $PLOTS_DIR/cell_type_distribution.png"
log "  - Spatial annotation map: $PLOTS_DIR/spatial_annotation.png"
log "  - Confidence scores: $PLOTS_DIR/annotation_confidence.png"
log "  - Individual weight maps: $PLOTS_DIR/cell_type_weights/"
log ""
log "To view results, check the plots directory: $PLOTS_DIR"