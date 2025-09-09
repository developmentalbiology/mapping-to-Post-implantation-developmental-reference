# Tutorial: Advanced Spatial Annotation with the 5-Tier Assignment Algorithm

This tutorial provides a comprehensive guide on using the spatial annotation workflow, from installation to final visualization. The workflow is designed to start with RCTD results, but we also provide optional utilities for data conversion and RCTD analysis.

## Table of Contents

1. [Installation](#installation)
2. [Quick Start with Example Data](#quick-start-with-example-data)
3. [Prerequisites](#prerequisites)
4. [Optional: Data Preparation](#optional-data-preparation)
5. [Main Workflow: 5-Tier Assignment Algorithm](#main-workflow-5-tier-assignment-algorithm)
6. [Configuration Parameters](#configuration-parameters)
7. [Step-by-Step Guide](#step-by-step-guide)
8. [Troubleshooting](#troubleshooting)

## Installation

### Download the Tutorial Package

**Option 1: Clone from GitHub (Recommended)**
```bash
git clone https://github.com/developmentalbiology/mapping-to-Post-implantation-developmental-reference.git
cd mapping-to-Post-implantation-developmental-reference/spatial_label_transfer_tutorial
```

**Option 2: Download ZIP**
```bash
wget https://github.com/developmentalbiology/mapping-to-Post-implantation-developmental-reference/archive/main.zip
unzip main.zip
cd mapping-to-Post-implantation-developmental-reference-main/spatial_label_transfer_tutorial
```

### Install Dependencies

```bash
# Install R packages
./install_r_packages.R

# Install Python packages (optional, only needed for H5AD conversion)
pip install -r requirements.txt
```

## Quick Start with Example Data

To reproduce the published human CS8 embryo analysis results:

### Step 1: Download Example Data

```bash
cd data/

# Download human CS8 spatial data
wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1VoQqgAsqxF3mXwl_HsbbAu64Jp-duXwh' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1VoQqgAsqxF3mXwl_HsbbAu64Jp-duXwh" -O cs8_human_embryo.rds

# Download human reference data  
wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1lhVc_tZCgXcWecW8aV7iINdDz696_r7r' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1lhVc_tZCgXcWecW8aV7iINdDz696_r7r" -O human_ref.rds

cd ..
```

**Manual Download Alternative:**
- Human CS8 data: [https://drive.google.com/file/d/1VoQqgAsqxF3mXwl_HsbbAu64Jp-duXwh/view?usp=drive_link](https://drive.google.com/file/d/1VoQqgAsqxF3mXwl_HsbbAu64Jp-duXwh/view?usp=drive_link)
- Human reference: [https://drive.google.com/file/d/1lhVc_tZCgXcWecW8aV7iINdDz696_r7r/view?usp=drive_link](https://drive.google.com/file/d/1lhVc_tZCgXcWecW8aV7iINdDz696_r7r/view?usp=drive_link)

### Step 2: Run RCTD Analysis (if needed)

```bash
Rscript scripts/utils/run_rctd_with_human_ref.R \
    --spatial data/cs8_human_embryo.rds \
    --reference data/human_ref.rds \
    --output data/rctd_results.rds
```

### Step 3: Run the Complete Workflow

```bash
./examples/human_embryo_example.sh
```

This will reproduce the exact results from the original human CS8 analysis using the same parameters and generate identical visualizations.

## Prerequisites

**Required for Main Workflow:**
- **RCTD Results:** An RDS file from your completed RCTD analysis
- **Spatial Data:** The original spatial data (RDS file) used for RCTD
- **R Environment:** All required packages installed via `install_r_packages.R`

**If you already have these files, skip to [Main Workflow](#main-workflow-5-tier-assignment-algorithm)**

## Optional: Data Preparation

**⚠️ Skip this section if you already have RCTD results in RDS format**

### Converting H5AD Files to RDS Format

If your data is in H5AD format (AnnData), you need to convert it to RDS format before running RCTD or this workflow.

#### Method 1: Using Our Conversion Scripts

**Step 1: Export H5AD Data**

```bash
python3 scripts/utils/export_h5ad.py \
    --input your_reference.h5ad \
    --output exported_data/ \
    --hvg-count 2000
```

This script exports:
- Expression matrix (MTX format)
- Cell barcodes and gene features
- Cell and gene metadata
- Highly variable genes
- UMAP coordinates (if available)

**Step 2: Convert to RDS**

```bash
Rscript scripts/utils/h5ad_to_rds.R \
    --input exported_data/ \
    --output reference.rds \
    --project "reference_data"
```

#### Method 2: Direct Conversion in R

You can also convert H5AD to RDS directly in R using SeuratDisk:

```R
library(Seurat)
library(SeuratDisk)

# Install SeuratDisk if not already installed
# remotes::install_github("mojaveazure/seurat-disk")

# Convert h5ad to h5seurat
Convert("reference.h5ad", dest = "h5seurat", overwrite = TRUE)

# Load as Seurat object
reference <- LoadH5Seurat("reference.h5seurat")

# Save as RDS
saveRDS(reference, "reference.rds")
```

### Running RCTD Analysis

If you haven't run RCTD analysis yet, this section provides detailed guidance on running RCTD with our provided human reference data or your own reference.

#### Option 1: Using Our Human Reference Data

We provide a pre-processed human reference dataset optimized for embryonic and developmental spatial transcriptomics analysis.

**Download the Reference Data**

```bash
# Method 1: Using wget (recommended)
wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1lhVc_tZCgXcWecW8aV7iINdDz696_r7r' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1lhVc_tZCgXcWecW8aV7iINdDz696_r7r" -O human_ref.rds && rm -rf /tmp/cookies.txt

# Method 2: Manual download
# Visit: https://drive.google.com/file/d/1lhVc_tZCgXcWecW8aV7iINdDz696_r7r/view?usp=drive_link
# Download and save as human_ref.rds
```

**Run RCTD with Human Reference**

```bash
# Using our provided script
Rscript scripts/utils/run_rctd_with_human_ref.R \
    --human-ref human_ref.rds \
    --spatial your_spatial_data.rds \
    --output rctd_results.rds \
    --spatial-coords "x" "y" \
    --max-cores 10  
```

**Manual RCTD Analysis**

If you prefer to run RCTD manually, here's the R code:

```R
library(Seurat)
library(spacexr)

# Load the human reference data
human_ref <- readRDS("human_ref.rds")

# Load your spatial data
spatial_data <- readRDS("your_spatial_data.rds")

# Extract reference information
ref_counts <- human_ref[["RNA"]]$counts
ref_cluster <- as.factor(human_ref$cell_type)  # adjust column name if needed
names(ref_cluster) <- colnames(human_ref)
ref_nUMI <- human_ref$nCount_RNA
names(ref_nUMI) <- colnames(human_ref)

# Create RCTD Reference object
reference <- Reference(ref_counts, ref_cluster, ref_nUMI)

# Extract spatial information
spatial_coords <- spatial_data@meta.data[, c("x", "y")]  # adjust column names
spatial_counts <- spatial_data[["RNA"]]$counts

# Create RCTD SpatialRNA object
query <- SpatialRNA(spatial_coords, spatial_counts, colSums(spatial_counts))

# Create and run RCTD
RCTD <- create.RCTD(query, reference, max_cores = 20)
RCTD <- run.RCTD(RCTD, doublet_mode = "full")

# Save results
saveRDS(RCTD, "rctd_results.rds")
```

#### Option 2: Using Your Own Reference Data

For detailed instructions on running RCTD with your own reference data, please refer to the official **[spacexr GitHub repository](https://github.com/dmcable/spacexr)**.

## Main Workflow: Weight Score-Based Assignment Refinement

**Starting Point:** You should now have RCTD results (RDS file) and original spatial data (RDS file).

### Understanding the Weight-Based Approach

The core principle of this workflow is that **all cell type assignments are fundamentally based on RCTD weight scores**. The 5-tier assignment algorithm refines these assignments by:

1. **Analyzing weight patterns** across all cell types and spatial locations
2. **Identifying problematic assignments** (over-assigned, under-assigned, or inconsistent)
3. **Applying intelligent rules** to optimize assignments based on weight characteristics
4. **Preserving biologically meaningful signals** while suppressing artifacts

### Three-Step Workflow

#### Step 1: Generate RCTD Weight Scores

RCTD produces a weight matrix where each cell (spatial location) has weight scores for all reference cell types. These weights represent the likelihood of each cell type being present at that location.

```bash
# If not already completed
Rscript scripts/utils/run_rctd_with_human_ref.R \
    --human-ref human_ref.rds \
    --spatial spatial_data.rds \
    --output rctd_results.rds
```

#### Step 2: Export and Analyze Weight Scores

**Critical Step:** Before making any assignments, analyze the weight score patterns to understand your data characteristics.

```bash
# Generate weight score visualizations
./scripts/run_workflow.sh \
    --rctd-results rctd_results.rds \
    --spatial-data spatial_data.rds \
    --output-dir weight_analysis/ \
    --skip-step 1  # Only generate weight plots, skip assignment
```

**Weight Analysis Checklist:**

1. **Examine individual cell type weight maps** (`weight_analysis/plots/cell_type_weights/`)
   - Look for spatial coherence vs. noise
   - Identify cell types with strong, clear patterns
   - Spot cell types with weak or scattered signals

2. **Review weight distribution histograms**
   - Identify cell types with high maximum weights (potentially over-assigned)
   - Find cell types with consistently low weights (potentially under-assigned)
   - Note cell types with bimodal distributions (mixed signals)

3. **Analyze spatial patterns**
   - Check if weight patterns match expected biology
   - Identify regions with conflicting high weights for multiple cell types
   - Look for rare cell types with meaningful but weak signals

#### Step 3: Configure Parameters Based on Weight Analysis

Based on your weight analysis, create a configuration file to optimize assignments:

```yaml
# config/optimized_config.yaml
annotation:
  # Exclude cell types that shouldn't exist in embryo proper
  excluded_cell_types: 
    - "TE"          # Trophoblast - extraembryonic lineage
    - "STB"         # Syncytiotrophoblast - not in embryo proper
    - "EVT"         # Extravillous trophoblast - placental cell type
  
  # Protect rare but biologically important cell types
  weak_signal_protected_types:
    - "Neural_crest"    # Critical developmental cell type
    - "Primitive_streak" # Transient but important organizer
  
  # Suppress over-dominant cell types based on weight analysis
  suppressed_cell_types:
    - "Epiblast_3"      # Too broad, overwhelming specific subtypes
  
  # Adjust thresholds based on weight distributions
  max_proportion: 0.3           # Reduce if seeing over-assignment
  spatial_weight: 0.8           # Increase for spatially coherent weights
  weak_signal_threshold: 0.05   # Lower for datasets with weak signals
```

#### Step 4: Apply Final Assignment

Run the complete workflow with optimized parameters:

```bash
./scripts/run_workflow.sh \
    --rctd-results rctd_results.rds \
    --spatial-data spatial_data.rds \
    --output-dir final_results/ \
    --config config/optimized_config.yaml
```

## Configuration Parameters: Weight Score-Guided Optimization

The configuration parameters should be adjusted based on your weight score analysis. Each parameter targets specific patterns observed in RCTD weight distributions.

### Weight Analysis-Based Parameter Selection

#### Cell Type Management Parameters

**Based on Weight Pattern Analysis:**

-   `excluded_cell_types`: Exclude cell types that shouldn't exist in your tissue/developmental context
    -   **When to use:** Cell types that are biologically inappropriate for your tissue or developmental stage
    -   **Examples:** Trophoblast cells in embryo proper analysis, adult cell types in embryonic data, tissue-specific cells in wrong tissue
    -   **Biological rationale:** Remove cell types that would be artifacts or contamination rather than true assignments
    -   **Example:** `["TE", "STB", "EVT"]` for embryo proper (excluding extraembryonic lineages)

-   `protected_cell_types`: Protect cell types with strong, clear weight patterns
    -   **When to use:** Cell types with high confidence weights that should never be filtered
    -   **Weight indicators:** High max weights (>0.5), spatially coherent patterns, biologically expected
    -   **Example:** `["Neural_crest", "Definitive_endoderm"]`

-   `weak_signal_protected_types`: Protect rare cell types with meaningful but weak weights
    -   **When to use:** Biologically important cell types with low but consistent weight patterns
    -   **Weight indicators:** Low max weights (0.1-0.3) but spatially coherent, matches expected biology
    -   **Example:** `["Rare_progenitor", "Transient_state"]`

-   `suppressed_cell_types`: Reduce dominance of over-assigned cell types
    -   **When to use:** Cell types with excessively high weights across many locations
    -   **Weight indicators:** High weights (>0.7) in too many locations, overwhelming other cell types
    -   **Example:** `["Abundant_celltype", "Over_represented"]`

#### Weight-Based Threshold Parameters

**Adjust based on weight score distributions:**

-   `max_proportion`: Maximum allowed proportion for any single cell type
    -   **Weight guidance:** If weight analysis shows one cell type dominating (>60% of locations), reduce this value
    -   **Typical range:** 0.2-0.6
    -   **Example:** `0.3` for datasets with many competing cell types

-   `spatial_weight`: Weight given to spatial consensus vs. RCTD weights
    -   **Weight guidance:** Increase if weight patterns are spatially coherent; decrease if weights are noisy
    -   **Typical range:** 0.4-0.9
    -   **Example:** `0.8` for clean, spatially organized weight patterns

-   `k_neighbors`: Number of neighbors for spatial consensus
    -   **Weight guidance:** Increase if weight patterns are smooth; decrease if patterns are highly localized
    -   **Typical range:** 3-10
    -   **Example:** `7` for moderately smooth weight transitions

#### Signal Balancing Parameters

**Fine-tune based on weight characteristics:**

-   `protection_threshold`: Minimum weight required to protect a cell type
    -   **Weight guidance:** Set based on the minimum meaningful weight in your analysis
    -   **Typical range:** 0.1-0.4
    -   **Example:** `0.2` for datasets with moderate signal strength

-   `weak_signal_threshold`: Lower threshold for weak signal cell types
    -   **Weight guidance:** Set to capture the lowest biologically meaningful weights
    -   **Typical range:** 0.02-0.15
    -   **Example:** `0.05` for datasets with very weak rare cell type signals

-   `suppression_strength`: Factor to reduce weights of suppressed cell types
    -   **Weight guidance:** Increase (closer to 0) for more aggressive suppression of dominant cell types
    -   **Typical range:** 0.5-0.9
    -   **Example:** `0.7` for moderate suppression of over-assigned cell types

### Weight Pattern Examples and Parameter Responses

#### Example 1: Over-Assigned Cell Type
**Weight Pattern:** Cell type "Epiblast" has weights >0.8 in 70% of locations
**Parameter Response:**
```yaml
suppressed_cell_types: ["Epiblast"]
max_suppressed_proportion: 0.3
suppression_strength: 0.6
```

#### Example 2: Weak Rare Cell Type
**Weight Pattern:** Cell type "Neural_crest" has max weight 0.15, spatially coherent
**Parameter Response:**
```yaml
weak_signal_protected_types: ["Neural_crest"]
weak_signal_threshold: 0.08
weak_signal_spatial_boost: 1.5
```

#### Example 3: Noisy Cell Type
**Weight Pattern:** Cell type "Artifact" has scattered weights, no spatial coherence
**Parameter Response:**
```yaml
excluded_cell_types: ["Artifact"]
```

## Step-by-Step Guide: Weight-Guided Parameter Optimization

### 1. Initial Weight Analysis

Start by generating weight score visualizations without any parameter optimization:

```bash
./scripts/run_workflow.sh \
    --rctd-results rctd_results.rds \
    --spatial-data spatial_data.rds \
    --output-dir initial_analysis/ \
    --skip-step 1  # Only generate weight plots
```

### 2. Analyze Weight Patterns

Examine the generated plots to understand your data:

**Key Files to Review:**
- `initial_analysis/plots/cell_type_weights/`: Individual weight maps for each cell type
- `initial_analysis/plots/weight_summary.png`: Overall weight distribution summary
- `initial_analysis/plots/spatial_overview.png`: Spatial patterns overview

**Analysis Questions:**
1. Which cell types show strong, spatially coherent patterns?
2. Which cell types have weak but meaningful signals?
3. Which cell types appear over-assigned (high weights everywhere)?
4. Which cell types show noisy, inconsistent patterns?

### 3. Create Initial Configuration

Based on your weight analysis, create a configuration file:

```yaml
# config/iteration1.yaml
annotation:
  # Start with obvious biological exclusions
  excluded_cell_types: ["TE", "STB"]  # Extraembryonic cell types
  weak_signal_protected_types: ["Neural_crest"]  # Important rare cell type
  suppressed_cell_types: ["Epiblast_general"]  # Over-broad category
  
  # Conservative initial parameters
  max_proportion: 0.5
  spatial_weight: 0.6
  k_neighbors: 5
```

### 4. Run First Iteration

```bash
./scripts/run_workflow.sh \
    --rctd-results rctd_results.rds \
    --spatial-data spatial_data.rds \
    --output-dir iteration1/ \
    --config config/iteration1.yaml
```

### 5. Evaluate Results and Refine

Examine the assignment results:

**Key Evaluation Plots:**
- `iteration1/plots/spatial_annotation.png`: Final spatial assignments
- `iteration1/plots/cell_type_distribution.png`: Cell type proportions
- `iteration1/plots/tier_distribution.png`: Assignment tier usage
- `iteration1/plots/annotation_confidence.png`: Assignment confidence

**Refinement Strategy:**
- If a cell type is still over-assigned → Add to `suppressed_cell_types` or reduce `max_proportion`
- If rare cell types are missing → Add to `weak_signal_protected_types` or lower `weak_signal_threshold`
- If annotations are too noisy → Increase `spatial_weight` or `k_neighbors`
- If annotations are too smooth → Decrease `spatial_weight` or `k_neighbors`

### 6. Iterate Configuration

Create refined configuration based on results:

```yaml
# config/iteration2.yaml
annotation:
  excluded_cell_types: ["TE", "STB", "EVT"]  # All extraembryonic lineages
  weak_signal_protected_types: ["Neural_crest", "Primitive_streak"]  # Rare but critical
  suppressed_cell_types: ["Epiblast_general", "Mesoderm_broad"]  # Over-broad categories
  
  # Refined parameters based on iteration 1 results
  max_proportion: 0.3  # Reduced due to over-assignment
  spatial_weight: 0.8  # Increased for better spatial coherence
  k_neighbors: 7       # Increased for smoother transitions
  
  # Fine-tuned thresholds
  weak_signal_threshold: 0.05  # Lowered to capture more rare types
  suppression_strength: 0.7    # Increased suppression
```

### 7. Final Iteration

```bash
./scripts/run_workflow.sh \
    --rctd-results rctd_results.rds \
    --spatial-data spatial_data.rds \
    --output-dir final_results/ \
    --config config/iteration2.yaml
```

### 8. Validation

Validate your final results:

1. **Biological Plausibility:** Do the assignments match expected biology?
2. **Spatial Coherence:** Are similar cell types clustered appropriately?
3. **Rare Cell Type Recovery:** Are biologically important rare cell types preserved?
4. **Weight Consistency:** Do final assignments align with strong weight patterns?

### Example Iteration Log

**Iteration 1 Observations:**
- "Epiblast_3" assigned to 60% of locations (too high)
- "Neural_crest" completely missing (too strict filtering)
- Spatial patterns look noisy

**Iteration 2 Adjustments:**
- Added "Epiblast_3" to suppressed_cell_types
- Added "Neural_crest" to weak_signal_protected_types  
- Increased spatial_weight from 0.6 to 0.8
- Reduced max_proportion from 0.5 to 0.3

**Iteration 2 Results:**
- "Epiblast_3" now 25% of locations (improved)
- "Neural_crest" recovered in expected regions
- Spatial coherence much improved

By iteratively adjusting parameters based on weight score analysis and assignment results, you can achieve highly accurate and robust spatial annotations tailored to your specific dataset.


### Verifying Data Quality

After conversion, verify your RDS file:

```R
# Load and inspect the converted data
data <- readRDS("reference.rds")

# Check basic properties
print(data)
print(dim(data))
print(colnames(data@meta.data))

# Check for required columns
if ("cell_type" %in% colnames(data@meta.data)) {
  print("Cell type column found")
  print(table(data$cell_type))
} else {
  print("Available metadata columns:")
  print(colnames(data@meta.data))
}
```

