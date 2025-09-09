# Data Acquisition Instructions

This directory is where you should place your input data files for the spatial annotation workflow.

## Required Files

### For the Main Workflow
You need either:

**Option 1: If you already have RCTD results**
- `rctd_results.rds` - Output from your RCTD analysis
- `spatial_data.rds` - Original spatial transcriptomics data (RDS format)

**Option 2: If you need to run RCTD first**  
- `spatial_data.rds` - Your spatial transcriptomics data
- `reference.rds` - Single-cell reference data (or use our human reference)

## Getting Example Data

### Human CS8 Embryo Example Data

To reproduce the exact results from the published human CS8 spatial transcriptomics analysis, you'll need:

1. **Human CS8 spatial data** - The original human CS8 embryo spatial transcriptomics dataset
2. **RCTD results** - Pre-computed RCTD results for the CS8 data using the human reference
3. **Human reference data** - Single-cell reference for running RCTD (if needed)

The provided `config/human_CS8_spatial.yaml` contains the exact parameters used in the original analysis to ensure reproducible results.

**Download Example Data:**

```bash
# Download human CS8 spatial data
wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1VoQqgAsqxF3mXwl_HsbbAu64Jp-duXwh' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1VoQqgAsqxF3mXwl_HsbbAu64Jp-duXwh" -O cs8_human_embryo.rds && rm -rf /tmp/cookies.txt

# Or download manually from:
# https://drive.google.com/file/d/1VoQqgAsqxF3mXwl_HsbbAu64Jp-duXwh/view?usp=drive_link

# Download human reference data  
wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1lhVc_tZCgXcWecW8aV7iINdDz696_r7r' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1lhVc_tZCgXcWecW8aV7iINdDz696_r7r" -O human_ref.rds && rm -rf /tmp/cookies.txt

# Or download manually from:
# https://drive.google.com/file/d/1lhVc_tZCgXcWecW8aV7iINdDz696_r7r/view?usp=drive_link
```

After downloading both files, you can run the complete workflow:

```bash
# First run RCTD (if needed)
Rscript ../scripts/utils/run_rctd_with_human_ref.R \
    --spatial cs8_human_embryo.rds \
    --reference human_ref.rds \
    --output rctd_results.rds

# Then run the annotation workflow
../examples/human_embryo_example.sh
```

### Option 1: Download Our Human Reference Dataset

```bash
# Download pre-processed human reference for embryonic analysis
wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1lhVc_tZCgXcWecW8aV7iINdDz696_r7r' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1lhVc_tZCgXcWecW8aV7iINdDz696_r7r" -O human_ref.rds && rm -rf /tmp/cookies.txt

# Or download manually from:
# https://drive.google.com/file/d/1lhVc_tZCgXcWecW8aV7iINdDz696_r7r/view?usp=drive_link
```

### Option 2: Use Your Own Data

Place your files in this directory and ensure they follow the expected format:

**Spatial Data Requirements:**
- Must be a Seurat RDS object
- Should contain spatial coordinates in metadata (columns "x", "y" or specify with `--spatial-coords`)
- Must have raw counts in the RNA assay

**Reference Data Requirements:**
- Must be a Seurat RDS object  
- Should have cell type annotations in metadata (column "cell_type")
- Must have raw counts in the RNA assay

## Data Conversion

### From H5AD Format

If your data is in H5AD (AnnData) format:

```bash
# Step 1: Export H5AD data
python3 ../scripts/utils/export_h5ad.py \
    --input your_data.h5ad \
    --output exported_data/ \
    --hvg-count 2000

# Step 2: Convert to RDS
Rscript ../scripts/utils/h5ad_to_rds.R \
    --input exported_data/ \
    --output your_data.rds \
    --project "your_project_name"
```

### Data Validation

After placing or converting your data, validate it:

```R
# Load and check your data
data <- readRDS("spatial_data.rds")
print(data)
print(dim(data))
print(colnames(data@meta.data))

# Check for spatial coordinates
if (all(c("x", "y") %in% colnames(data@meta.data))) {
  print("✓ Spatial coordinates found")
} else {
  print("Available coordinate columns:")
  print(colnames(data@meta.data))
}

# For reference data, check cell type annotations
if ("cell_type" %in% colnames(data@meta.data)) {
  print("✓ Cell type annotations found")
  print(table(data$cell_type))
}
```

## File Naming Conventions

For the example scripts to work automatically, use these names:
- `rctd_results.rds` - RCTD analysis results
- `spatial_data.rds` - Spatial transcriptomics data  
- `human_ref.rds` - Human reference data (if using our reference)
- `reference.rds` - Your own reference data

Or modify the paths in the example scripts to match your file names.

## Directory Structure After Data Acquisition

```
data/
├── README.md                # This file
├── rctd_results.rds        # RCTD analysis output
├── spatial_data.rds        # Spatial transcriptomics data
└── human_ref.rds           # Reference data (optional)
```

## Troubleshooting

### Common Issues
1. **File not found errors**: Check that file paths match what's specified in scripts
2. **Memory issues**: Large datasets may require more RAM or processing in chunks
3. **Format issues**: Ensure RDS files contain Seurat objects with proper structure
4. **Coordinate columns**: Spatial coordinates must be named correctly ("x", "y") or specified with `--spatial-coords`

### Getting Help
- Check the main [TUTORIAL.md](../TUTORIAL.md) for detailed guidance
- Review example scripts in [examples/](../examples/) for proper usage patterns
- Ensure all required R packages are installed with `../install_r_packages.R`