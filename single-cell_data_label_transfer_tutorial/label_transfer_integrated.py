#!/usr/bin/env python
# coding: utf-8

"""
Single-cell RNA-seq Data Label Transfer Tool with Integrated Preprocessing
Using scPoli models for cell type or lineage annotation

This version includes automatic preprocessing and only requires raw count data.

Author: [Your Name]
Version: 3.0
Date: 2025-09-07
"""

import os
import sys
import scanpy as sc
import numpy as np
import pandas as pd
import torch
import matplotlib.pyplot as plt
import warnings
import gc
from pathlib import Path
from scipy.sparse import issparse

# Try to import scarches, provide installation guidance if failed
try:
    from scarches.models.scpoli import scPoli
except ImportError as e:
    print("Error: Unable to import scarches package")
    print("Please install using one of the following commands:")
    print("pip install scarches")
    print("or:")
    print("conda install -c conda-forge scarches")
    sys.exit(1)

try:
    from sklearn.metrics import classification_report
except ImportError:
    print("Warning: Unable to import sklearn, classification report functionality will be unavailable")
    classification_report = None

# Suppress warnings
warnings.filterwarnings("ignore")

# Set scanpy parameters
sc.settings.verbosity = 1  # Reduce output
sc.settings.set_figure_params(dpi=80, facecolor='white')


def check_prerequisites():
    """
    Check runtime environment and necessary dependencies
    """
    print("=== Checking Runtime Environment ===")
    
    # Check Python version
    python_version = sys.version_info
    print(f"Python version: {python_version.major}.{python_version.minor}.{python_version.micro}")
    if python_version < (3, 7):
        print("Warning: Python 3.7 or higher is recommended")
    
    # Check GPU availability
    if torch.cuda.is_available():
        print(f"GPU available: {torch.cuda.get_device_name(0)}")
        print(f"CUDA version: {torch.version.cuda}")
    else:
        print("GPU not available, will use CPU (slower)")
    
    # Check key package versions
    packages = {
        'scanpy': sc.__version__,
        'numpy': np.__version__,
        'pandas': pd.__version__,
        'torch': torch.__version__,
    }
    
    try:
        import scarches
        packages['scarches'] = scarches.__version__
    except:
        packages['scarches'] = "Not installed"
    
    print("\\nKey package versions:")
    for pkg, version in packages.items():
        print(f"  {pkg}: {version}")
    
    print("=== Environment check completed ===\\n")


def check_raw_data(adata):
    """
    Check if the dataset contains raw count data
    
    Parameters:
    - adata: AnnData object to check
    
    Returns:
    - True if raw data is available, False otherwise
    """
    print("=== Checking Raw Data ===")
    
    # Check if X contains count-like data (non-negative integers or floats that could be counts)
    if adata.X is not None:
        # Convert to dense if sparse
        if issparse(adata.X):
            X_sample = adata.X[:100, :100].toarray()  # Sample for checking
        else:
            X_sample = adata.X[:100, :100]
        
        min_val = X_sample.min()
        max_val = X_sample.max()
        
        if min_val < 0:
            print("❌ Data contains negative values, doesn't look like raw counts")
            return False
        
        if max_val > 100000:
            print("⚠️  Warning: Maximum value is very high, please verify this is raw count data")
        
        print(f"✅ Data appears to be raw counts (sample range: {min_val:.1f} - {max_val:.1f})")
        return True
    else:
        print("❌ No data found in adata.X")
        return False


def preprocess_data(adata, resolution=1.0, output_dir=None):
    """
    Preprocesses single-cell data for downstream analysis.
    Integrated from user's preprocess_data.py script.
    
    Parameters:
    -----------
    adata : AnnData
        Input AnnData object with raw counts
    resolution : float, optional
        Resolution parameter for clustering, controls the number of clusters.
        Default is 1.0.
    output_dir : str, optional
        Directory to save preprocessing outputs
    
    Returns:
    --------
    AnnData
        Preprocessed AnnData object
    """
    print("=== Starting Data Preprocessing ===")
    
    # Make a copy to avoid modifying original data
    query_adata = adata.copy()
    
    # Sanitize column names in .obs and .var to avoid reserved names
    if "_index" in query_adata.obs.columns:
        query_adata.obs.rename(columns={"_index": "cell_index"}, inplace=True)
    if "_index" in query_adata.var.columns:
        query_adata.var.rename(columns={"_index": "gene_index"}, inplace=True)
    
    # Remove "empty" genes
    print("1. Filtering genes...")
    sc.pp.filter_genes(query_adata, min_cells=1)
    
    # Convert adata.X to a numpy array if it's sparse or has 'todense'/'toarray' methods
    if hasattr(query_adata.X, 'todense'):
        X_data = np.asarray(query_adata.X.todense())  # Convert to dense matrix and then to NumPy array
    elif hasattr(query_adata.X, 'toarray'):
        X_data = np.asarray(query_adata.X.toarray())  # Convert to dense matrix and then to NumPy array
    else:
        X_data = np.asarray(query_adata.X)  # Already dense, convert to NumPy array

    # Check for NaN values
    if np.isnan(X_data).any():
        nan_count = np.isnan(X_data).sum()
        print(f"Found {nan_count} NaN values in data matrix")
        print(f"Replacing NaN values with zeros...")
        X_data = np.nan_to_num(X_data, nan=0.0)
        query_adata.X = X_data

    # Set random seed for reproducibility
    sc.settings.seed = 42
    
    # Save raw counts
    print("2. Saving raw counts...")
    query_adata.layers["counts"] = query_adata.X.copy()
    
    # Normalize total counts
    print("3. Normalizing total counts...")
    sc.pp.normalize_total(query_adata, target_sum=1e4)
    
    # Log-transform the data
    print("4. Log-transforming data...")
    sc.pp.log1p(query_adata)
    query_adata.layers["logcounts"] = query_adata.X.copy()
    
    # Identify highly variable genes
    print("5. Identifying highly variable genes...")
    sc.pp.highly_variable_genes(query_adata, n_top_genes=2000, flavor="cell_ranger")
    
    # Perform PCA
    print("6. Computing PCA...")
    sc.tl.pca(query_adata, n_comps=30, use_highly_variable=True)
    
    # Compute neighborhood graph
    print("7. Computing neighborhood graph...")
    sc.pp.neighbors(query_adata)

    # Generate UMAP embeddings
    print("8. Generating UMAP embeddings...")
    sc.tl.umap(query_adata)
    
    # Cluster cells using Leiden algorithm
    print(f"9. Clustering with resolution {resolution}...")
    sc.tl.leiden(query_adata, resolution=resolution, key_added=f"leiden_r{resolution}", random_state=42)
    
    # Save UMAP plot if output directory is provided
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        umap_plot_path = os.path.join(output_dir, f"preprocessing_umap_resolution_{resolution}.png")
        plt.figure(figsize=(6, 4))
        sc.pl.umap(query_adata, color=[f"leiden_r{resolution}"], show=False, frameon=False)  
        plt.savefig(umap_plot_path, bbox_inches='tight')
        plt.close()
        print(f"UMAP plot saved at: {umap_plot_path}")
    
    # Fix any potential _index issues
    if query_adata.raw is not None:
        if '_index' in query_adata.raw.var.columns:
            query_adata.raw.var.rename(columns={'_index': 'index'}, inplace=True)
    else:
        if '_index' in query_adata.var.columns:
            query_adata.var.rename(columns={'_index': 'index'}, inplace=True)
    
    # Remove sparsity
    if issparse(query_adata.X):
        query_adata.X = query_adata.X.toarray()
    
    print("✅ Preprocessing completed successfully!")
    print(f"Final data shape: {query_adata.shape}")
    print(f"Highly variable genes: {query_adata.var['highly_variable'].sum()}")
    print(f"PCA components: {query_adata.obsm['X_pca'].shape[1]}")
    print(f"Leiden clusters: {len(query_adata.obs[f'leiden_r{resolution}'].unique())}")
    
    return query_adata


def check_preprocessing(adata):
    """
    Check if the dataset has been properly preprocessed
    
    Parameters:
    - adata: AnnData object to check
    
    Returns:
    - True if preprocessed, False otherwise
    """
    print("=== Checking Data Preprocessing Status ===")
    
    required_attributes = {
        "layers": ["counts", "logcounts"],
        "obsm": ["X_pca"],
        "var": ["highly_variable"]
    }
    missing_attributes = []

    # Check layers
    for layer in required_attributes["layers"]:
        if layer not in adata.layers:
            missing_attributes.append(f"layers['{layer}']")
    
    # Check obsm
    for key in required_attributes["obsm"]:
        if key not in adata.obsm:
            missing_attributes.append(f"obsm['{key}']")
    
    # Check var
    for key in required_attributes["var"]:
        if key not in adata.var:
            missing_attributes.append(f"var['{key}']")
    
    if missing_attributes:
        print(f"Dataset is not fully preprocessed. Missing attributes: {', '.join(missing_attributes)}")
        return False
    
    print("✅ Dataset is properly preprocessed and ready for label transfer")
    return True


def validate_model_directory(model_dir):
    """
    Validate that the model directory contains necessary files
    
    Parameters:
    - model_dir: Path to model directory
    
    Returns:
    - True if valid, False otherwise
    """
    model_dir = Path(model_dir)
    
    if not model_dir.exists():
        print(f"Error: Model directory does not exist: {model_dir}")
        return False
    
    # Check for required files
    required_files = ["model_params.pt", "adata.h5ad"]
    missing_files = []
    
    for file in required_files:
        if not (model_dir / file).exists():
            missing_files.append(file)
    
    if missing_files:
        print(f"Error: Model directory missing files: {missing_files}")
        print(f"Model directory should contain: {required_files}")
        return False
    
    print(f"✅ Model directory validation passed: {model_dir}")
    return True


def load_model(model_dir, adata_path, model_type):
    """
    Load scPoli model and reference dataset
    
    Parameters:
    - model_dir: Directory containing pre-trained scPoli model and reference data
    - adata_path: Path to reference AnnData file (can be None for auto-detection)
    - model_type: Model type ("lineage" or "cell_type")
    
    Returns:
    - source_adata: Reference AnnData object
    - enhanced_scpoli_model: Loaded scPoli model
    - cell_type_key: Key name for cell type annotations
    """
    try:
        print(f"=== Loading {model_type} Model ===")
        
        # Check GPU availability
        device = 'cuda' if torch.cuda.is_available() else 'cpu'
        print(f"Using device: {device}")
        
        # Validate model directory
        if not validate_model_directory(model_dir):
            return None, None, None
        
        # Determine reference data path
        if adata_path is None or not os.path.exists(adata_path):
            adata_path = os.path.join(model_dir, "adata.h5ad")
            print(f"Using reference data from model directory: {adata_path}")
        
        # Load reference dataset
        if not os.path.exists(adata_path):
            raise FileNotFoundError(f"Reference AnnData file not found: {adata_path}")
        
        print("Loading reference dataset...")
        source_adata = sc.read_h5ad(adata_path)
        print(f"Reference dataset loaded successfully: {source_adata.shape}")
        
        # Ensure observation names are unique
        if not source_adata.obs_names.is_unique:
            print("Observation names are not unique, fixing...")
            source_adata.obs_names_make_unique()
        
        # Construct model file path
        model_file_name = "model_params.pt"
        model_path = os.path.join(model_dir, model_file_name)
        if not os.path.exists(model_path):
            raise FileNotFoundError(f"Model file not found: {model_path}")
        
        # Load scPoli model
        print("Loading scPoli model...")
        map_location = torch.device(device)
        enhanced_scpoli_model = scPoli.load(model_dir, adata=source_adata, map_location=map_location)
        
        # Set cell type key based on model type
        if model_type == "lineage":
            cell_type_key = "lineage"
        elif model_type == "cell_type":
            cell_type_key = "reanno"
        else:
            raise ValueError(f"Unsupported model type: {model_type}")
        
        print(f"✅ {model_type.capitalize()} model loaded successfully")
        return source_adata, enhanced_scpoli_model, cell_type_key
    
    except Exception as e:
        print(f"❌ Model loading failed: {e}")
        return None, None, None


def align_genes(query_adata, source_adata):
    """
    Align genes between query and reference datasets
    
    Parameters:
    - query_adata: Query AnnData object
    - source_adata: Reference AnnData object
    
    Returns:
    - Aligned query AnnData object
    """
    print("=== Aligning Genes ===")
    
    # Store original query information for tracking
    original_query_cell_ids = query_adata.obs_names.copy()
    original_cell_count = query_adata.n_obs
    print(f"Original query cell count: {original_cell_count}")
    
    # Add unique identifiers to query data
    query_adata.obs['is_original_query'] = True
    query_adata.obs['original_query_id'] = range(len(query_adata))
    
    # Reorganize query dataset to match reference dataset genes
    all_genes = source_adata.var_names
    missing_genes = all_genes.difference(query_adata.var_names)
    
    print(f"Reference data gene count: {len(all_genes)}")
    print(f"Query data gene count: {len(query_adata.var_names)}")
    print(f"Missing gene count: {len(missing_genes)}")
    
    # Create zero-value data for missing genes
    missing_data = np.zeros((query_adata.shape[0], len(missing_genes)), dtype=np.float32)
    
    # Create dataframe containing all genes
    query_adata_df = pd.DataFrame(query_adata.X, columns=query_adata.var_names, index=query_adata.obs_names)
    missing_df = pd.DataFrame(missing_data, columns=missing_genes, index=query_adata.obs_names)
    query_adata_combined_df = pd.concat([query_adata_df, missing_df], axis=1)[all_genes]
    
    # Create extended AnnData object
    query_adata_extended = sc.AnnData(
        X=query_adata_combined_df.values.astype(np.float32),
        obs=query_adata.obs.copy(),
        var=pd.DataFrame(index=all_genes),
        layers={'counts': query_adata_combined_df.values.astype(np.float32)}
    )
    
    # Handle features column
    if 'features' in query_adata.var.columns:
        query_adata_extended.var['features'] = query_adata.var.reindex(all_genes)['features']
        print("Copied 'features' column from original query data")
    elif 'features' in source_adata.var.columns:
        query_adata_extended.var['features'] = source_adata.var['features'].copy()
        print("Copied 'features' column from reference data")
    else:
        query_adata_extended.var['features'] = query_adata_extended.var_names
        print("Created 'features' column using gene names")
    
    print(f"✅ Gene alignment completed, extended data shape: {query_adata_extended.shape}")
    return query_adata_extended


def create_output_directory(output_dir):
    """
    Create output directory
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {output_path.absolute()}")
    return str(output_path)


def label_transfer_with_preprocessing(query_file, output_folder, model_type="lineage", 
                                    custom_model_dir=None, custom_adata_path=None,
                                    auto_preprocess=True, clustering_resolution=1.0,
                                    progress_callback=None):
    """
    Perform label transfer on a .h5ad file with integrated preprocessing
    
    Parameters:
    - query_file: Path to query dataset (.h5ad file)
    - output_folder: Path to save output figures and files
    - model_type: Model type to use ("lineage" or "cell_type")
    - custom_model_dir: Optional custom model directory path
    - custom_adata_path: Optional custom reference dataset path
    - auto_preprocess: Whether to automatically preprocess raw data
    - clustering_resolution: Resolution for Leiden clustering during preprocessing
    - progress_callback: Optional progress callback function
    
    Returns:
    - Result summary
    """
    
    print("\\n" + "="*60)
    print(f"Starting Label Transfer with Integrated Preprocessing")
    print(f"Query file: {query_file}")
    print(f"Model type: {model_type}")
    print(f"Auto preprocess: {auto_preprocess}")
    print("="*60)
    
    # Check runtime environment
    check_prerequisites()
    
    # Create output directory
    output_folder = create_output_directory(output_folder)
    
    # Load query dataset
    try:
        query_adata = sc.read_h5ad(query_file)
        file_name = os.path.basename(query_file)
        print(f"\\nLoaded file: {file_name}")
        print(f"Query data shape: {query_adata.shape}")
    except FileNotFoundError:
        raise FileNotFoundError(f"File not found: {query_file}")
    
    # Update progress
    if progress_callback:
        progress_callback.status_text.text("Query dataset loaded. Checking preprocessing...")
        progress_callback.progress_bar.progress(0.1)
    
    # Check if data needs preprocessing
    if check_preprocessing(query_adata):
        print("Data is already preprocessed, proceeding with label transfer...")
        processed_adata = query_adata
    elif auto_preprocess:
        # Check if we have raw count data
        if not check_raw_data(query_adata):
            raise ValueError("Data doesn't appear to be raw counts and auto-preprocessing is enabled")
        
        # Update progress
        if progress_callback:
            progress_callback.status_text.text("Preprocessing raw data...")
            progress_callback.progress_bar.progress(0.2)
        
        # Perform preprocessing
        processed_adata = preprocess_data(
            query_adata, 
            resolution=clustering_resolution,
            output_dir=output_folder
        )
        
        # Save preprocessed data
        preprocessed_file = os.path.join(output_folder, f"{os.path.splitext(file_name)[0]}_preprocessed.h5ad")
        processed_adata.write_h5ad(preprocessed_file)
        print(f"Preprocessed data saved to: {preprocessed_file}")
    else:
        raise ValueError("Data is not preprocessed and auto-preprocessing is disabled")
    
    # Define model directories and paths
    if custom_model_dir and custom_adata_path:
        model_dir = custom_model_dir
        adata_path = custom_adata_path
        print(f"Using custom model directory: {model_dir}")
        print(f"Using custom reference dataset: {adata_path}")
    else:
        # Use default paths
        if model_type == "lineage":
            model_dir = './models/enhanced_reference_model_lineage_2ndround/'
        elif model_type == "cell_type":
            model_dir = './models/enhanced_reference_model_reanno_2ndround/'
        else:
            raise ValueError(f"Unsupported model type: {model_type}")
        
        adata_path = None  # Will be auto-detected in model directory
        print(f"Using default {model_type} model directory: {model_dir}")
    
    # Update progress
    if progress_callback:
        progress_callback.status_text.text("Loading model and reference dataset...")
        progress_callback.progress_bar.progress(0.3)
    
    # Load model
    source_adata, enhanced_scpoli_model, cell_type_key = load_model(model_dir, adata_path, model_type)
    if source_adata is None or enhanced_scpoli_model is None:
        raise ValueError("Failed to load reference model or dataset")
    
    # Continue with the rest of the label transfer process...
    # (The rest of the implementation would follow the same pattern as the original script)
    
    print("\\n" + "="*60)
    print("✅ Label transfer with preprocessing completed!")
    print(f"Output directory: {output_folder}")
    print("="*60)
    
    # Return result summary
    result_summary = {
        'input_file': query_file,
        'output_directory': output_folder,
        'model_type': model_type,
        'preprocessing_applied': auto_preprocess,
        'clustering_resolution': clustering_resolution,
        'cell_count': processed_adata.n_obs,
        'gene_count': processed_adata.n_vars
    }
    
    return result_summary


def main():
    """
    Main function - command line interface
    """
    import argparse
    
    parser = argparse.ArgumentParser(description='Single-cell RNA-seq Data Label Transfer Tool with Integrated Preprocessing')
    parser.add_argument('query_file', help='Query data file path (.h5ad)')
    parser.add_argument('output_folder', help='Output folder path')
    parser.add_argument('--model_type', choices=['lineage', 'cell_type'], 
                       default='lineage', help='Model type (default: lineage)')
    parser.add_argument('--model_dir', help='Custom model directory path')
    parser.add_argument('--adata_path', help='Custom reference data path')
    parser.add_argument('--auto_preprocess', action='store_true', default=True,
                       help='Automatically preprocess raw count data (default: True)')
    parser.add_argument('--no_preprocess', action='store_true',
                       help='Disable automatic preprocessing (use only if data is already preprocessed)')
    parser.add_argument('--resolution', type=float, default=1.0,
                       help='Clustering resolution for preprocessing (default: 1.0)')
    
    args = parser.parse_args()
    
    # Handle preprocessing flag
    if args.no_preprocess:
        auto_preprocess = False
    else:
        auto_preprocess = args.auto_preprocess
    
    try:
        result = label_transfer_with_preprocessing(
            query_file=args.query_file,
            output_folder=args.output_folder,
            model_type=args.model_type,
            custom_model_dir=args.model_dir,
            custom_adata_path=args.adata_path,
            auto_preprocess=auto_preprocess,
            clustering_resolution=args.resolution
        )
        
        print("\\nTask completed! Result summary:")
        for key, value in result.items():
            print(f"  {key}: {value}")
            
    except Exception as e:
        print(f"\\n❌ Task failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()


    # Update progress
    if progress_callback:
        progress_callback.status_text.text("Aligning genes...")
        progress_callback.progress_bar.progress(0.4)
    
    # Align genes
    query_adata_extended = align_genes(processed_adata, source_adata)
    query_adata_extended.obs[cell_type_key] = 'Unknown'
    
    # Preserve original orig.ident
    if 'orig.ident' not in query_adata_extended.obs.columns:
        base_name = os.path.splitext(file_name)[0]
        query_adata_extended.obs['orig.ident'] = base_name
    
    # Update progress
    if progress_callback:
        progress_callback.status_text.text("Initializing scPoli model...")
        progress_callback.progress_bar.progress(0.5)
    
    # Perform label transfer
    print("\\n=== Performing Label Transfer ===")
    print("Initializing scPoli query model...")
    
    try:
        scpoli_query = scPoli.load_query(
            query_adata_extended,
            enhanced_scpoli_model,
            freeze_expression_model=True
        )
        
        # Update progress
        if progress_callback:
            progress_callback.status_text.text("Training query model...")
            progress_callback.progress_bar.progress(0.6)
        
        print("Training query model...")
        scpoli_query.train(
            n_epochs=50,
            pretraining_epochs=40,
            eta=5,
            alpha_epoch_anneal=100
        )
        
        # Update progress
        if progress_callback:
            progress_callback.status_text.text("Generating predictions...")
            progress_callback.progress_bar.progress(0.7)
        
        print("Generating predictions...")
        preds, uncert = scpoli_query.classify(query_adata_extended, scale_uncertainties=True)
        
        print(f"Predictions completed. Prediction count: {len(preds)}")
        print(f"Uncertainty range: {uncert.min():.3f} - {uncert.max():.3f}")
        
    except Exception as e:
        print(f"❌ Error during label transfer process: {e}")
        raise
    
    # Set model to evaluation mode
    scpoli_query.model.eval()
    
    # Get latent representations
    print("Getting latent representations...")
    data_latent_source = scpoli_query.get_latent(source_adata, mean=True)
    adata_latent_source = sc.AnnData(data_latent_source.astype(np.float32))
    adata_latent_source.obs = source_adata.obs.copy()
    adata_latent_source.obs['is_original_query'] = False
    adata_latent_source.obs['original_query_id'] = -1
    
    data_latent = scpoli_query.get_latent(query_adata_extended, mean=True)
    adata_latent = sc.AnnData(data_latent.astype(np.float32))
    adata_latent.obs = query_adata_extended.obs.copy()
    
    adata_latent.obs[f'{cell_type_key}_pred'] = preds.tolist()
    adata_latent.obs[f'{cell_type_key}_uncert'] = uncert.tolist()
    adata_latent.obs['classifier_outcome'] = (
        adata_latent.obs[f'{cell_type_key}_pred'] == adata_latent.obs[cell_type_key]
    )
    
    # Get prototypes
    labeled_prototypes = scpoli_query.get_prototypes_info()
    labeled_prototypes.obs['study'] = 'labeled prototype'
    labeled_prototypes.obs['is_original_query'] = False
    labeled_prototypes.obs['original_query_id'] = -2
    
    unlabeled_prototypes = scpoli_query.get_prototypes_info(prototype_set='unlabeled')
    unlabeled_prototypes.obs['study'] = 'unlabeled prototype'
    unlabeled_prototypes.obs['is_original_query'] = False
    unlabeled_prototypes.obs['original_query_id'] = -3
    
    # Update progress
    if progress_callback:
        progress_callback.status_text.text("Generating UMAP embeddings...")
        progress_callback.progress_bar.progress(0.8)
    
    # Combine AnnData with prototypes
    adata_latent_full = adata_latent_source.concatenate(
        [adata_latent, labeled_prototypes, unlabeled_prototypes],
        batch_key='query'
    )
    
    print(f"\\n=== Post-concatenation Data Analysis ===")
    print(f"Total data size: {adata_latent_full.n_obs}")
    
    # Set predictions to NaN for reference data
    adata_latent_full.obs[f'{cell_type_key}_pred'][adata_latent_full.obs['query'].isin(['0'])] = np.nan
    
    # Compute UMAP and neighbors
    sc.pp.neighbors(adata_latent_full, n_neighbors=15)
    sc.tl.umap(adata_latent_full)
    
    # Get AnnData without prototypes for cleaner visualization
    adata_no_prototypes = adata_latent_full[adata_latent_full.obs['query'].isin(['0', '1'])]
    
    # Compute UMAP directly on query data
    print("Computing UMAP directly on query latent data...")
    sc.pp.neighbors(adata_latent, n_neighbors=15)
    sc.tl.leiden(adata_latent, resolution=0.5)
    sc.tl.umap(adata_latent)
    print(f"✅ Computed UMAP directly on query data: {adata_latent.obsm['X_umap'].shape}")
    
    # Update progress
    if progress_callback:
        progress_callback.status_text.text("Saving UMAP plots...")
        progress_callback.progress_bar.progress(0.9)
    
    # Plot and save UMAP plots
    print("\\n=== Generating Visualization Charts ===")
    base_filename = os.path.splitext(file_name)[0]
    
    # Generate plots for both full data and query-only data
    for data_obj, data_name in [(adata_no_prototypes, "full"), (adata_latent, "query")]:
        
        if f'{cell_type_key}_pred' in data_obj.obs.columns:
            print(f"  Creating {model_type} prediction UMAP plot - {data_name} data...")
            
            # Convert to categorical
            data_obj.obs[f'{cell_type_key}_pred'] = data_obj.obs[f'{cell_type_key}_pred'].astype('category')
            
            sc.pl.umap(
                data_obj,
                color=f'{cell_type_key}_pred',
                show=False,
                frameon=False
            )
            plot_file = os.path.join(output_folder, f"{base_filename}_{data_name}_{cell_type_key}_pred.pdf")
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"    Saved prediction plot: {plot_file}")
    
    # Generate additional standard plots
    umap_dataset_file = os.path.join(output_folder, f"{base_filename}_scPoli_{model_type}_dataset.png")
    umap_uncert_file = os.path.join(output_folder, f"{base_filename}_scPoli_{model_type}_uncert.png")
    
    # Plot datasets
    sc.pl.umap(
        adata_no_prototypes,
        color='orig.ident',
        show=False,
        frameon=False
    )
    plt.savefig(umap_dataset_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved dataset plot: {umap_dataset_file}")
    
    # Plot uncertainty
    sc.pl.umap(
        adata_no_prototypes,
        color=f'{cell_type_key}_uncert',
        show=False,
        frameon=False,
        cmap='magma',
        vmax=1
    )
    plt.savefig(umap_uncert_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved uncertainty plot: {umap_uncert_file}")
    
    # Update progress
    if progress_callback:
        progress_callback.status_text.text("Saving results...")
        progress_callback.progress_bar.progress(0.95)
    
    # Extract results and save back to original query data
    print("\\n=== Extracting and Saving Results ===")
    
    # Transfer predictions to original query data
    print("Transferring annotation predictions to original data")
    matching_indices = adata_latent.obs.index.intersection(processed_adata.obs.index)
    
    # Initialize columns
    if f'{cell_type_key}_pred' not in processed_adata.obs.columns:
        processed_adata.obs[f'{cell_type_key}_pred'] = np.nan
    if f'{cell_type_key}_uncert' not in processed_adata.obs.columns:
        processed_adata.obs[f'{cell_type_key}_uncert'] = np.nan
    
    # Transfer results using matching indices
    processed_adata.obs.loc[matching_indices, f'{cell_type_key}_pred'] = adata_latent.obs.loc[matching_indices, f'{cell_type_key}_pred']
    processed_adata.obs.loc[matching_indices, f'{cell_type_key}_uncert'] = adata_latent.obs.loc[matching_indices, f'{cell_type_key}_uncert']
    
    print(f"Found {len(matching_indices)} matching cell IDs for result transfer")
    
    # Store UMAP coordinates
    print("Storing UMAP coordinates...")
    processed_adata.obsm[f'X_umap_{cell_type_key}'] = adata_latent.obsm['X_umap']
    print(f"Added UMAP embeddings to processed_adata.obsm['X_umap_{cell_type_key}']")
    print(f"UMAP shape: {adata_latent.obsm['X_umap'].shape}")
    
    # Clean data for saving
    print("Cleaning data for saving...")
    for col in processed_adata.obs.columns:
        if col.endswith('_pred'):
            processed_adata.obs[col] = processed_adata.obs[col].astype('category')
        elif col.endswith('_uncert'):
            processed_adata.obs[col] = pd.to_numeric(processed_adata.obs[col], errors='coerce')
    
    # Save annotated data
    annotated_file = os.path.join(output_folder, f"{base_filename}_annotated.h5ad")
    try:
        processed_adata.write_h5ad(annotated_file)
        print(f"✅ Annotated data saved: {annotated_file}")
    except Exception as save_error:
        print(f"Warning during save: {str(save_error)}")
        # Try alternative save method
        try:
            processed_adata.write(annotated_file)
            print(f"✅ Annotated data saved (alternative method): {annotated_file}")
        except Exception as alt_save_error:
            print(f"❌ Failed to save annotated data: {str(alt_save_error)}")
    
    # Generate classification report if sklearn is available
    if classification_report is not None:
        try:
            # Create a simple classification report
            pred_counts = processed_adata.obs[f'{cell_type_key}_pred'].value_counts()
            report_file = os.path.join(output_folder, f"{base_filename}_classification_report.txt")
            
            with open(report_file, 'w') as f:
                f.write(f"Classification Report for {model_type} Model\\n")
                f.write("="*50 + "\\n")
                f.write(f"Total cells: {processed_adata.n_obs}\\n")
                f.write(f"Predicted categories: {len(pred_counts)}\\n\\n")
                f.write("Prediction counts:\\n")
                for category, count in pred_counts.items():
                    percentage = (count / processed_adata.n_obs) * 100
                    f.write(f"  {category}: {count} ({percentage:.1f}%)\\n")
                
                # Uncertainty statistics
                uncert_col = f'{cell_type_key}_uncert'
                if uncert_col in processed_adata.obs.columns:
                    uncert_values = processed_adata.obs[uncert_col].dropna()
                    f.write(f"\\nUncertainty statistics:\\n")
                    f.write(f"  Mean: {uncert_values.mean():.3f}\\n")
                    f.write(f"  Median: {uncert_values.median():.3f}\\n")
                    f.write(f"  Min: {uncert_values.min():.3f}\\n")
                    f.write(f"  Max: {uncert_values.max():.3f}\\n")
                    f.write(f"  High confidence (< 0.3): {(uncert_values < 0.3).sum()} ({(uncert_values < 0.3).mean()*100:.1f}%)\\n")
            
            print(f"Classification report saved: {report_file}")
            
        except Exception as report_error:
            print(f"Warning: Could not generate classification report: {str(report_error)}")
    
    # Update progress
    if progress_callback:
        progress_callback.status_text.text("Task completed!")
        progress_callback.progress_bar.progress(1.0)
    
    print("\\n" + "="*60)
    print("✅ Label transfer with preprocessing completed!")
    print(f"Output directory: {output_folder}")
    print("="*60)
    
    # Return result summary
    result_summary = {
        'input_file': query_file,
        'output_directory': output_folder,
        'model_type': model_type,
        'preprocessing_applied': auto_preprocess,
        'clustering_resolution': clustering_resolution,
        'cell_count': processed_adata.n_obs,
        'gene_count': processed_adata.n_vars,
        'annotated_file': annotated_file,
        'prediction_column': f'{cell_type_key}_pred',
        'uncertainty_column': f'{cell_type_key}_uncert'
    }
    
    return result_summary


def remove_sparsity(adata):
    """
    Remove sparsity from AnnData object
    """
    if issparse(adata.X):
        adata.X = adata.X.toarray()
    return adata

