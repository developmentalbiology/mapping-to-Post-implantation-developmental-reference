#!/bin/bash

# =============================================================================
# Human CS8 Embryo Spatial Annotation Workflow (Flexible I/O)
#
# Allows full control over input and output paths.
# No more forcing users into the workflow directory!
# =============================================================================

set -euo pipefail

# -------------------------------
# Default values
# -------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_DIR="$(dirname "$SCRIPT_DIR")"

# ✅ 默认输出目录：当前目录下的 output/human_CS8_example
DEFAULT_OUTPUT_DIR="./output/human_CS8_example"

# 其他默认值（可被覆盖）
DEFAULT_CONFIG_FILE="$WORKFLOW_DIR/config/human_CS8_spatial.yaml"
RCTD_RESULTS=""
SPATIAL_DATA=""
OUTPUT_DIR="$DEFAULT_OUTPUT_DIR"
CONFIG_FILE="$DEFAULT_CONFIG_FILE"
SPATIAL_COORDS_X="x"
SPATIAL_COORDS_Y="y"
VERBOSE=true

# -------------------------------
# Help message
# -------------------------------
print_help() {
    cat << 'EOF'
Usage: human_embryo_example.sh [OPTIONS]

Run spatial annotation with flexible input/output paths.
You can specify your own data files and output location.

Options:
    --rctd-results <file>     Path to RCTD results .rds file (required)
    --spatial-data <file>     Path to spatial data .rds file (required)
    --output-dir <dir>        Output directory (default: ./output/human_CS8_example)
    --config <file>           Configuration YAML (default: internal config)
    --spatial-coords <x> <y>  Column names for spatial coordinates (default: "x" "y")
    --quiet                   Disable verbose output
    --help                    Show this help and exit

Example:
    ./human_embryo_example.sh \
        --rctd-results ~/data/rctd.rds \
        --spatial-data ~/data/spatial.rds \
        --output-dir ./results/cs8_analysis
EOF
}

# -------------------------------
# Parse arguments
# -------------------------------
while [[ $# -gt 0 ]]; do
    case $1 in
        --rctd-results)
            RCTD_RESULTS="$2"; shift; shift ;;
        --spatial-data)
            SPATIAL_DATA="$2"; shift; shift ;;
        --output-dir)
            OUTPUT_DIR="$2"; shift; shift ;;
        --config)
            CONFIG_FILE="$2"; shift; shift ;;
        --spatial-coords)
            SPATIAL_COORDS_X="$2"; SPATIAL_COORDS_Y="$3"; shift; shift; shift ;;
        --quiet)
            VERBOSE=false; shift ;;
        --help|-h)
            print_help; exit 0 ;;
        *)
            echo "Unknown option: $1" >&2; print_help; exit 1 ;;
    esac
done

# -------------------------------
# Validate required inputs
# -------------------------------
if [[ -z "$RCTD_RESULTS" ]]; then
    echo "ERROR: Missing --rctd-results" >&2
    print_help
    exit 1
fi

if [[ -z "$SPATIAL_DATA" ]]; then
    echo "ERROR: Missing --spatial-data" >&2
    print_help
    exit 1
fi

if [[ ! -f "$RCTD_RESULTS" ]]; then
    echo "ERROR: RCTD file not found: $RCTD_RESULTS" >&2
    exit 1
fi

if [[ ! -f "$SPATIAL_DATA" ]]; then
    echo "ERROR: Spatial data file not found: $SPATIAL_DATA" >&2
    exit 1
fi

if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "WARNING: Config file not found: $CONFIG_FILE" >&2
fi

# ✅ 创建输出目录（基于用户指定或默认）
mkdir -p "$OUTPUT_DIR"

# -------------------------------
# Run the workflow
# -------------------------------
echo "=== Human CS8 Embryo Spatial Annotation Workflow ==="
echo "Using flexible I/O paths."
echo ""
echo "RCTD Results:     $RCTD_RESULTS"
echo "Spatial Data:     $SPATIAL_DATA"
echo "Output Directory: $OUTPUT_DIR"
echo "Config File:      $CONFIG_FILE"
echo "Spatial Coords:   ($SPATIAL_COORDS_X, $SPATIAL_COORDS_Y)"
echo ""

# 检查 run_workflow.sh 是否可执行
WORKFLOW_SCRIPT="$WORKFLOW_DIR/scripts/run_workflow.sh"
if [[ ! -x "$WORKFLOW_SCRIPT" ]]; then
    echo "ERROR: Workflow script not executable: $WORKFLOW_SCRIPT" >&2
    echo "Run: chmod +x '$WORKFLOW_SCRIPT'" >&2
    exit 1
fi

echo "Running analysis..."
"$WORKFLOW_SCRIPT" \
    --rctd-results "$RCTD_RESULTS" \
    --spatial-data "$SPATIAL_DATA" \
    --output-dir "$OUTPUT_DIR" \
    --config "$CONFIG_FILE" \
    --spatial-coords "$SPATIAL_COORDS_X" "$SPATIAL_COORDS_Y" \
    ${VERBOSE:+--verbose}

echo ""
echo "✅ Analysis completed!"
echo "Results saved to: $(realpath "$OUTPUT_DIR")"
echo ""
