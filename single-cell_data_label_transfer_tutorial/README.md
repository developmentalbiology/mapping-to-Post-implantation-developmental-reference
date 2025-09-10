# Tutorial: Single-cell Data Label Transfer using scPoli

---

## 1. Introduction

This tutorial will guide you through using the `label_transfer_integrated.py` script to perform cell type and lineage annotation on your own single-cell RNA sequencing (scRNA-seq) data using pre-trained scPoli models. 

**Key Advantage**: This version **only requires raw count data** and automatically handles all preprocessing steps internally.

### What You Need:
- **Input**: Raw count data in `.h5ad` format (just the basic count matrix)
- **Output**: Fully annotated data with cell type and lineage predictions and visualizations

By the end of this tutorial, you will be able to:

1. Set up the Python environment required to run the script
2. Download the necessary models and test data
3. Run the script with just raw count data
4. Understand and interpret the output files

---

## 2. Environment Setup

To successfully run the label transfer script, you need a properly configured Python environment with scArches (which includes scPoli). We recommend following the official scArches installation instructions for the most reliable setup.

### 2.1. System Requirements

- **Python**: 3.7 or 3.8 (as specified in the [scArches documentation](https://docs.scarches.org/en/latest/installation.html))
- **Operating System**: Linux, macOS, or Windows
- **Memory**: At least 8GB RAM recommended for typical datasets
- **GPU**: Optional but recommended for faster training (NVIDIA GPU with CUDA support)

### 2.2. Installation Methods

#### Method 1: Using pip (Recommended for most users)

The easiest way to install scArches (which includes scPoli) is through pip:

```bash
# Create a new virtual environment (recommended)
python -m venv scpoli_env
source scpoli_env/bin/activate  # On Windows: scpoli_env\Scripts\activate

# Install scArches with all dependencies
pip install -U scarches

# Install additional dependencies for our script
pip install matplotlib scikit-learn
```

#### Method 2: Using the official scArches conda environment

For a more controlled environment setup, you can use the official environment file from the scArches repository:

```bash
# Clone the scArches repository
git clone https://github.com/theislab/scarches
cd scarches

# Create conda environment from the official environment file
conda env create -f envs/scarches_linux.yaml
conda activate scarches
```

This method installs the exact package versions tested by the scArches developers, including:
- Python 3.10
- PyTorch with CUDA support
- scvi-tools
- All required dependencies

#### Method 3: Manual conda installation

If you prefer to set up your own conda environment:

```bash
# Create a new conda environment
conda create -n scpoli_env python=3.8
conda activate scpoli_env

# Install PyTorch (visit https://pytorch.org/get-started/locally/ for the latest commands)
# For CPU-only:
conda install pytorch torchvision torchaudio cpuonly -c pytorch

# For GPU (example - adjust CUDA version as needed):
# conda install pytorch torchvision torchaudio pytorch-cuda=11.8 -c pytorch -c nvidia

# Install scArches and additional dependencies
pip install scarches matplotlib scikit-learn
```

### 2.3. Verification

After installation, verify that everything is working correctly:

```bash
# Activate your environment
conda activate scpoli_env  # or: source scpoli_env/bin/activate

# Test the installation
python -c "
import scarches
from scarches.models.scpoli import scPoli
import scanpy as sc
import torch
print('✅ All packages imported successfully!')
print(f'scArches version: {scarches.__version__}')
print(f'PyTorch version: {torch.__version__}')
print(f'CUDA available: {torch.cuda.is_available()}')
"
```

You can also test the script's built-in environment check:

```bash
python label_transfer_improved_en.py --help
```

### 2.4. Additional Resources

- **scArches Documentation**: https://docs.scarches.org/
- **scArches GitHub Repository**: https://github.com/theislab/scarches
- **scPoli Paper**: [Nature Methods](https://www.nature.com/articles/s41592-023-02035-2)

Following these installation steps will ensure you have a properly configured environment for running the label transfer script with scPoli.


## 3. Model and Data Preparation

Before running label transfer, you need to download the script, pre-trained models, and reference datasets. You can also prepare your own data as query input.

### 3.1. Downloading the Script

First, you need to obtain the `label_transfer_integrated.py` script:

1. **Download the script**: [label_transfer_integrated.py](label_transfer_integrated.py) - Save this file to your project directory
2. **Make it executable** (optional on Linux/Mac):
   ```bash
   chmod +x label_transfer_integrated.py
   ```

### 3.2. Downloading Models and Reference Data

The script requires pre-trained scPoli models and corresponding reference AnnData files (`adata.h5ad`). These files are packaged together in the same directory.

1. **Access model download link**: 
   [https://drive.google.com/file/d/1rsNPoOKvBn8tOLM8Pr3yxUltp7AxT-GQ/view?usp=drive_link](https://drive.google.com/file/d/1rsNPoOKvBn8tOLM8Pr3yxUltp7AxT-GQ/view?usp=drive_link)

2. **Download model folders**: You will see two folders corresponding to two model types:
   - `enhanced_reference_model_lineage_2ndround`: For **lineage** annotation
   - `enhanced_reference_model_reanno_2ndround`: For **cell_type** annotation

   You can right-click on the folder you need and select "Download".

3. **Organize file structure**: After downloading, put the `models` folder, script, and query data in the same directory:

    ```
    your_project_directory/
    ├── label_transfer_integrated.py         # Downloaded script
    ├── models/                              # Downloaded models
    │   ├── enhanced_reference_model_lineage_2ndround/
    │   │   ├── adata.h5ad
    │   │   └── model_params.pt
    │   └── enhanced_reference_model_reanno_2ndround/
    │       ├── adata.h5ad
    │       └── model_params.pt
    └── my_raw_data.h5ad                             # Query data for annotation
    ```

    **Important Note**: Each model folder must contain both `adata.h5ad` (reference data) and `model_params.pt` (model weights).

### 3.3. Downloading Example Dataset (Optional)

To quickly test if the script works properly, you can download a pre-processed example dataset.

1. **Access example data link**:
   [https://drive.google.com/file/d/1w1rohJSr4xVm_b0o0b84pB_TUXS44XsS/view?usp=drive_link](https://drive.google.com/file/d/1w1rohJSr4xVm_b0o0b84pB_TUXS44XsS/view?usp=drive_link)

2. **Download and test**: Save the file to your `data` directory and test the script for 'lineage' and 'cell_type' annotation:

   ```bash
   python label_transfer_integrated.py \
       Chen_2025_GSE262081.h5ad \
       results/test_annotation 
   ```


### 3.4. Preparing Your Own Query Data

#### What You Need:
- An AnnData (`.h5ad`) file with raw count data in `adata.X`
- Non-negative values representing gene expression counts
- Basic gene and cell identifiers


#### Example: Loading Raw 10x Data

```python
import scanpy as sc

# Load raw 10x data
adata = sc.read_10x_mtx(
    'path/to/your/10x/data/',
    var_names='gene_symbols',
    cache=True
)

# Make variable names unique (if needed)
adata.var_names_unique()

# Save as h5ad file - that's it!
adata.write_h5ad('my_raw_data.h5ad')
```

#### Checking Your Data

You can verify your data is suitable:

```python
import scanpy as sc
import numpy as np

# Load your data
adata = sc.read_h5ad('my_raw_data.h5ad')

# Check basic properties
print(f"Data shape: {adata.shape}")
print(f"Data type: {adata.X.dtype}")
print(f"Min value: {adata.X.min()}")
print(f"Max value: {adata.X.max()}")

# Raw count data should be:
# - Non-negative (min >= 0)
# - Integer or float values
# - Reasonable range (typically 0-10000+ for highly expressed genes)
```

---

## 4. Running the Script

### 4.1. Command Line Arguments

- **`query_file`** (required): Path to your raw count `.h5ad` file
- **`output_folder`** (required): Folder to save all results
- **`--auto_preprocess`** (optional): Enable automatic preprocessing (default: `True`)
- **`--resolution`** (optional): Clustering resolution for preprocessing (default: `1.0`)
- **`--model_dir_lineage`** (optional): Custom model directory path for lineage
- **`--model_dir_celltype`** (optional): Custom model directory path for celltype
- **`--adata_path_lineage`** (optional): Custom reference data path for lineage
- **`--adata_path_celltype`** (optional): Custom reference data path for celltype
- **`--no_preprocess`** (optional): Disable preprocessing (only if data is already processed)

### 4.2. Basic Usage Examples

#### Example 1: annotation with provided dataset

```bash
python label_transfer_integrated.py \
    Chen_2025_GSE262081.h5ad \
    results/annotation  
```

#### Example 2: annotation with your own raw data

```bash
python label_transfer_integrated.py \
    my_data.h5ad \
    results/annotation
```


### 4.3. What Happens During Execution

The script automatically performs these steps:

1. **Environment Check**: Verifies all dependencies are installed
2. **Data Loading**: Loads your raw count data
3. **Data Validation**: Checks if data appears to be raw counts
4. **Automatic Preprocessing**:
   - Gene filtering (removes genes expressed in <1 cell)
   - Normalization (target sum = 10,000)
   - Log transformation
   - Highly variable gene identification (top 2000 genes)
   - PCA computation (30 components)
   - Neighborhood graph construction
   - UMAP embedding generation
   - Leiden clustering
5. **Model Loading**: Loads the specified scPoli model
6. **Gene Alignment**: Matches genes between your data and reference
7. **Label Transfer**: Performs the actual annotation
8. **Results Generation**: Creates output files and visualizations

---

## 5. Understanding the Output

After successful execution, you'll find these files in your output folder:

### 5.1. Main Output Files

- **`*_preprocessed.h5ad`**: Your data after preprocessing (includes all intermediate results)
- **`*_annotated.h5ad`**: Final annotated data with predictions
- **`preprocessing_umap_resolution_*.png`**: UMAP plot showing preprocessing clusters
- **`*_prediction.pdf`**: UMAP plot with predicted labels
- **`*_uncertainty.pdf`**: Uncertainty visualization

### 5.2. Key Data Columns Added

In the final annotated file, you'll find:

- **`lineage_pred` or `reanno_pred`**: Predicted cell type/lineage for each cell
- **`lineage_uncert` or `reanno_uncert`**: Prediction uncertainty (0-1, lower is better)
- **`leiden_r*`**: Clustering results from preprocessing
- **`obsm['X_umap_*']`**: UMAP coordinates for visualization

---

## 6. Advanced Usage

### 6.1. Custom Model Paths

```bash
python label_transfer_integrated.py \
    my_data.h5ad \
    results/annotation \
    --model_dir_lineage /path/to/custom/model/lineage \
    --model_dir_celltype /path/to/custom/model/celltype
```

### 6.2. Batch Processing

```bash
#!/bin/bash
# Process multiple files
for file in data/*.h5ad; do
    filename=$(basename "$file" .h5ad)
    python label_transfer_integrated.py \
        "$file" \
        "results/${filename}_annotation" 
done
```

---

## 7. Summary

This simplified workflow makes single-cell annotation much easier:

1. **Install scArches**: `pip install scarches`
2. **Download script**: Get `label_transfer_integrated.py` and place in your project directory
3. **Download models**: Get pre-trained models from the provided Google Drive link
4. **Prepare raw data**: Just need counts in `.h5ad` format - no preprocessing required!
5. **Run script**: `python label_transfer_integrated.py data.h5ad results/ `
6. **Get results**: Annotated data with cell type predictions and visualizations

---

**Author Information**: This tutorial was created by Manus AI to simplify single-cell data annotation using scPoli. For questions or improvements, please refer to the [scArches documentation](https://docs.scarches.org/) or [GitHub repository](https://github.com/theislab/scarches).

