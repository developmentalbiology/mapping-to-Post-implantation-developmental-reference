#!/bin/bash

# Example: Running Spatial Annotation Workflow for Human CS8 Embryo Data
# 
# This script demonstrates how to reproduce the exact results from the 
# original human CS8 spatial transcriptomics analysis using the provided
# configuration parameters that match the published research.

set -e

# Get the directory of this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_DIR="$(dirname "$SCRIPT_DIR")"

# Human CS8 embryo data paths - modify these paths to point to your actual data
RCTD_RESULTS="$WORKFLOW_DIR/data/rctd_results.rds"          # Your RCTD results file
SPATIAL_DATA="$WORKFLOW_DIR/data/spatial_data.rds"          # Your human CS8 spatial data file  
OUTPUT_DIR="$WORKFLOW_DIR/output/human_CS8_example"
CONFIG_FILE="$WORKFLOW_DIR/config/human_CS8_spatial.yaml"   # Updated config file name

echo "=== Human CS8 Embryo Spatial Annotation Workflow ==="
echo "This reproduces the exact analysis parameters from the original research"
echo ""
echo "RCTD results: $RCTD_RESULTS"
echo "Spatial data: $SPATIAL_DATA" 
echo "Output directory: $OUTPUT_DIR"
echo "Configuration: $CONFIG_FILE"
echo ""

# Check if input files exist
if [ ! -f "$RCTD_RESULTS" ]; then
    echo "Warning: RCTD results file not found: $RCTD_RESULTS"
    echo "Please place your human CS8 RCTD results file at this location."
    echo "See data/README.md for data acquisition instructions."
    echo ""
fi

if [ ! -f "$SPATIAL_DATA" ]; then
    echo "Warning: Spatial data file not found: $SPATIAL_DATA"
    echo "Please place your human CS8 spatial data file at this location."
    echo "See data/README.md for data acquisition instructions."
    echo ""
fi

# Run the workflow with human CS8 specific parameters
echo "Running human CS8 embryo spatial annotation workflow..."
echo "Using exact parameters from the original analysis to reproduce results..."
"$WORKFLOW_DIR/scripts/run_workflow.sh" \
    --rctd-results "$RCTD_RESULTS" \
    --spatial-data "$SPATIAL_DATA" \
    --output-dir "$OUTPUT_DIR" \
    --config "$CONFIG_FILE" \
    --spatial-coords "x" "y" \
    --verbose

echo ""
echo "=== Human CS8 Analysis Completed ==="
echo "Results saved to: $OUTPUT_DIR"
echo ""
echo "Key outputs:"
echo "  - Final annotations: $OUTPUT_DIR/annotated_spatial.rds"
echo "  - Final assignment plot: $OUTPUT_DIR/annotated_spatial_final_assignment.pdf"
echo "  - Assignment results CSV: $OUTPUT_DIR/annotated_spatial_assignment_results.csv"
echo "  - Cell type weight plots: $OUTPUT_DIR/plots/all_cell_type_weights.pdf"
echo "  - Weight statistics: $OUTPUT_DIR/plots/weight_statistics.csv"
echo ""
echo "The results should match the original human CS8 analysis exactly."