# Label Transfer Tutorials for Developmental Biology

**Advanced computational workflows for cell type annotation in developmental biology research**

This repository provides two complementary tutorials for cutting-edge cell type annotation and label transfer methods, specifically developed and optimized for developmental biology applications including embryonic and post-implantation developmental studies.

## 🧬 Repository Overview

### Our Methods

| Tutorial | Data Type | Description | Key Innovations |
|----------|-----------|-------------|-----------------|
| [**Single-Cell Label Transfer**](single-cell_data_label_transfer_tutorial/) | Single-cell RNA-seq | Advanced methods for single-cell data analysis and cell type annotation | Multi-reference integration, robust batch correction, confidence scoring |
| [**Spatial Label Transfer**](spatial_label_transfer_tutorial/) | Spatial transcriptomics | Novel 5-tier assignment algorithm for spatial transcriptomics annotation | RCTD-based refinement, spatial consensus, optimized weight scoring |

## 🚀 Quick Start

### Choose Your Analysis Path

**For Single-Cell Data:**
```bash
cd single-cell_data_label_transfer_tutorial/
./install_r_packages.R
# Follow tutorial for your specific use case
```

**For Spatial Transcriptomics Data:**
```bash
cd spatial_label_transfer_tutorial/
./install_r_packages.R
./examples/human_embryo_example.sh  # Reproduce published results
```

## 🔬 Applications in Developmental Biology

### Supported Research Areas
- **Embryonic development** - Early embryo cell fate mapping
- **Post-implantation development** - Gastrulation and organogenesis
- **Tissue morphogenesis** - Spatial organization analysis  
- **Cell fate transitions** - Developmental trajectory analysis
- **Cross-species comparisons** - Comparative developmental biology

### Example Use Cases
- Map cell types in human embryo spatial data (CS7, CS8)
- Integrate multiple developmental atlases
- Annotate new developmental datasets with established references
- Study spatial organization during morphogenesis
- Cross-reference single-cell and spatial measurements

## 📁 Repository Structure

```
mapping-to-Post-implantation-developmental-reference/
├── README.md                                    # This overview
├── single-cell_data_label_transfer_tutorial/   # Single-cell methods
│   ├── README.md                               # Single-cell quick start
│   ├── scripts/                                # Analysis pipelines
│   ├── config/                                 # Configuration files
│   └── examples/                               # Usage examples
└── spatial_label_transfer_tutorial/            # Spatial methods  
    ├── README.md                               # Spatial quick start
    ├── scripts/                                # Analysis pipelines
    ├── config/                                 # Parameter configurations
    ├── examples/                               # Human CS8 example
    └── data/                                   # Data download instructions
```

## 🛠️ Prerequisites

### System Requirements
- **R** (≥ 4.0.0) with development tools
- **Python** (≥ 3.7) - optional, for H5AD conversion utilities
- **Memory**: 8GB+ RAM recommended for large datasets
- **Storage**: Variable depending on dataset size

### Common Dependencies
Both tutorials share these core packages:
- `Seurat` - Single-cell and spatial data analysis
- `tidyverse` - Data manipulation and visualization  
- `Matrix` - Sparse matrix operations
- `ggplot2` - Advanced plotting

*See individual tutorial directories for complete dependency lists*

## 📖 Getting Started

### 1. Choose Your Tutorial
- **New to label transfer?** Start with the [single-cell tutorial](single-cell_data_label_transfer_tutorial/)
- **Have spatial data?** Go directly to the [spatial tutorial](spatial_label_transfer_tutorial/)
- **Want to compare methods?** Try both with the same biological system

### 2. Installation
```bash
# Clone the repository
git clone https://github.com/developmentalbiology/spatial_data_label_transfer_tutorial.git
cd spatial_data_label_transfer_tutorial

# Choose your tutorial directory
cd single-cell_data_label_transfer_tutorial/  # OR
cd spatial_label_transfer_tutorial/

# Install dependencies
./install_r_packages.R
```

### 3. Run Examples
Each tutorial includes working examples with real developmental biology datasets:
- **Single-cell**: Multi-atlas integration workflows
- **Spatial**: Human CS8 embryo analysis (exact parameter reproduction)

## 🎯 Reproducible Research

### Published Examples
- **Human CS8 Spatial Analysis** - Complete reproduction of published spatial transcriptomics analysis
- **Cross-species Comparisons** - Methods for comparative developmental studies  
- **Multi-timepoint Integration** - Developmental trajectory mapping

### Parameter Configurations
- Pre-optimized configurations for common developmental biology applications
- YAML-based parameter files for reproducibility
- Detailed parameter explanations and tuning guides

## 🤝 Community & Support

### Getting Help
1. **Check the documentation** - Each tutorial has comprehensive guides
2. **Review examples** - Working examples for common scenarios
3. **Open an issue** - For bugs, questions, or feature requests

### Contributing
We welcome contributions from the developmental biology community:
- Report issues or bugs
- Suggest new features or applications
- Share parameter configurations for new systems
- Contribute example datasets or use cases

---

**Get started with label transfer for your developmental biology research today!** 

Choose your data type above and follow the tutorial that matches your analysis needs. Both workflows are designed to be accessible to biologists while providing the depth needed for rigorous computational analysis.

*More details and additional content will be disclosed in the near future...*
