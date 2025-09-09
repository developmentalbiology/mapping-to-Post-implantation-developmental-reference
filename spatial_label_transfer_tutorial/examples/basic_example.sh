#!/bin/bash

# Basic Example: Running Spatial Annotation Workflow
# 
# This script demonstrates how to run the spatial annotation workflow
# with your own RCTD results and spatial data.

set -e

# Get the directory of this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_DIR="$(dirname "$SCRIPT_DIR")"

# Data paths - these point to the new flattened structure
RCTD_RESULTS="$WORKFLOW_DIR/data/rctd_results.rds"
SPATIAL_DATA="$WORKFLOW_DIR/data/spatial_data.rds" 
OUTPUT_DIR="$WORKFLOW_DIR/output/basic_example"
CONFIG_FILE="$WORKFLOW_DIR/config/default.yaml"

echo "=== Spatial Annotation Workflow Example ==="
echo "RCTD results: $RCTD_RESULTS"
echo "Spatial data: $SPATIAL_DATA"
echo "Output directory: $OUTPUT_DIR"
echo "Configuration: $CONFIG_FILE"
echo ""

# Check if input files exist
if [ ! -f "$RCTD_RESULTS" ]; then
    echo "Warning: RCTD results file not found: $RCTD_RESULTS"
    echo "Please place your RCTD results file at this location."
    echo "See data/README.md for data acquisition instructions."
    echo "To generate RCTD results, see: https://github.com/dmcable/spacexr"
    echo ""
fi

if [ ! -f "$SPATIAL_DATA" ]; then
    echo "Warning: Spatial data file not found: $SPATIAL_DATA"
    echo "Please place your spatial data file at this location."
    echo "See data/README.md for data acquisition instructions."
    echo ""
fi

# Run the workflow
echo "Running spatial annotation workflow..."
"$WORKFLOW_DIR/scripts/run_workflow.sh" \
    --rctd-results "$RCTD_RESULTS" \
    --spatial-data "$SPATIAL_DATA" \
    --output-dir "$OUTPUT_DIR" \
    --config "$CONFIG_FILE" \
    --verbose

echo ""
echo "=== Workflow completed ==="
echo "Results saved to: $OUTPUT_DIR"
echo ""
echo "Generated files:"
echo "  - Annotated spatial data: $OUTPUT_DIR/annotated_spatial.rds"
echo "  - Assignment results CSV: $OUTPUT_DIR/annotated_spatial_assignment_results.csv"
echo "  - Plots: $OUTPUT_DIR/plots/"
echo ""
echo "Key visualizations:"
echo "  - Cell type distribution: $OUTPUT_DIR/plots/cell_type_distribution.png"
echo "  - Spatial annotation map: $OUTPUT_DIR/plots/spatial_annotation.png"
echo "  - Confidence scores: $OUTPUT_DIR/plots/annotation_confidence.png"
echo "  - Individual weight maps: $OUTPUT_DIR/plots/cell_type_weights/"
echo ""
echo "To view the results, check the plots directory:"
echo "  ls $OUTPUT_DIR/plots/"