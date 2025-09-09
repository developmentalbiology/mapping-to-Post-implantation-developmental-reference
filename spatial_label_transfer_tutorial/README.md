# Spatial Data Label Transfer Tutorial

**Advanced spatial annotation using a 5-tier assignment algorithm for spatial transcriptomics data**

This tutorial provides a complete workflow for refining spatial transcriptomics cell type annotations using RCTD (Robust Cell Type Decomposition) weight scores and a multi-tier assignment algorithm.

## Installation

### Option 1: Clone from GitHub
```bash
git clone https://github.com/developmentalbiology/mapping-to-Post-implantation-developmental-reference.git
cd mapping-to-Post-implantation-developmental-reference/spatial_label_transfer_tutorial
```

### Option 2: Download ZIP
```bash
# Download and extract the repository
wget https://github.com/developmentalbiology/mapping-to-Post-implantation-developmental-reference/archive/main.zip
unzip main.zip
cd mapping-to-Post-implantation-developmental-reference-main/spatial_label_transfer_tutorial
```

### Install Dependencies
```bash
# Install R packages
./install_r_packages.R

# Install Python packages (if using H5AD conversion utilities)
pip install -r requirements.txt
```

## Quick Start

**Prerequisites:** RCTD results (RDS file) + original spatial data (RDS file)

```bash
# 1. Install dependencies
./install_r_packages.R

# 2. Run the workflow
./scripts/run_workflow.sh \
    --rctd-results your_rctd_results.rds \
    --spatial-data your_spatial_data.rds \
    --output-dir results/

# 3. View results in results/plots/
```

## Key Features

- **Weight Score-Based Refinement:** Uses RCTD weight patterns to guide parameter optimization
- **5-Tier Assignment Algorithm:** Sophisticated multi-tier approach for robust cell type assignments  
- **Comprehensive Visualization:** Spatial weight maps, annotation plots, and detailed statistics
- **Highly Configurable:** YAML-based configuration for fine-tuned control
- **Complete Pipeline:** Includes utilities for H5AD conversion and RCTD analysis

## Workflow Overview

The workflow operates in three main steps:
1. **Generate RCTD Weight Scores** - Initial cell type likelihood scores for each location
2. **Analyze Weight Patterns** - Examine weight distributions to guide parameter selection  
3. **Apply Refined Assignment** - Use multi-tier algorithm with optimized parameters

*Note: See [workflow_diagram_v2.puml](docs/workflow_diagram_v2.puml) for the detailed workflow diagram source.*

## Directory Structure

```
spatial_label_transfer_tutorial/
├── README.md                    # This file - quick start guide
├── TUTORIAL.md                  # Detailed methodology and parameter guide
├── LICENSE
├── install_r_packages.R         # R package installation
├── requirements.txt             # Python requirements
├── config/                      # Configuration files
│   ├── default.yaml
│   └── human_CS8_spatial.yaml
├── scripts/                     # Main workflow scripts  
│   ├── run_workflow.sh          # Main execution script
│   ├── step1_apply_annotation.R
│   ├── step2_export_plots.R
│   ├── multi_tier_integrated_final.R
│   └── utils/                   # Utility scripts
│       ├── export_h5ad.py
│       ├── h5ad_to_rds.R
│       └── run_rctd_with_human_ref.R
├── examples/                    # Example usage scripts
│   ├── basic_example.sh
│   └── human_embryo_example.sh   # Human CS8 embryo analysis example
├── docs/                        # Additional documentation
│   ├── TUTORIAL.md              # Detailed methodology guide
│   └── workflow_diagram_v2.puml # Workflow diagram source
├── data/                        # Data directory (see data/README.md)
└── output/                      # Output directory
```

## Usage Examples

### Basic Usage
```bash
./examples/basic_example.sh
```

### Human CS8 Embryo Analysis (Reproducible Example)
```bash
# Download example data first:
# Human CS8 spatial data: https://drive.google.com/file/d/1VoQqgAsqxF3mXwl_HsbbAu64Jp-duXwh/view?usp=drive_link
# Human reference data: https://drive.google.com/file/d/1lhVc_tZCgXcWecW8aV7iINdDz696_r7r/view?usp=drive_link

./examples/human_embryo_example.sh  # Reproduces human CS8 analysis exactly
```

### Custom Configuration
```bash
./scripts/run_workflow.sh \
    --rctd-results rctd_results.rds \
    --spatial-data spatial_data.rds \
    --output-dir results/ \
    --config config/human_CS8_spatial.yaml
```

## Data Preparation

If you need to convert H5AD files or run RCTD analysis first, see:
- [TUTORIAL.md](TUTORIAL.md) - Complete data preparation guide
- [data/README.md](data/README.md) - Data acquisition instructions

## Documentation

- **[TUTORIAL.md](TUTORIAL.md)** - Complete step-by-step tutorial with installation and usage
- **[docs/TUTORIAL.md](docs/TUTORIAL.md)** - Detailed methodology guide and technical details  
- **[data/README.md](data/README.md)** - Data acquisition and preparation instructions

## Output Files

The workflow generates:
- `annotated_spatial.rds` - Final annotated spatial data
- `annotated_spatial_assignment_results.csv` - Assignment results table
- `plots/` directory containing:
  - Spatial annotation maps
  - Individual cell type weight plots  
  - Summary statistics and distributions
  - Confidence score visualizations
