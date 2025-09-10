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


def label_transfer_with_preprocessing_dual_model(
    query_file,
    output_folder,
    custom_model_dir_lineage=None,
    custom_model_dir_celltype=None,
    custom_adata_path_lineage=None,
    custom_adata_path_celltype=None,
    auto_preprocess=True,
    clustering_resolution=1.0,
    progress_callback=None
):
    """
    Perform dual-model label transfer using BOTH lineage and cell_type models.
    For lineage model, uses the UMAP layout from cell_type model for consistent visualization.
    All results are saved in a single .h5ad file.

    Parameters:
    - query_file: Path to query dataset (.h5ad file)
    - output_folder: Output directory
    - custom_model_dir_lineage: Optional custom lineage model directory
    - custom_model_dir_celltype: Optional custom cell_type model directory
    - custom_adata_path_lineage: Optional custom lineage reference adata
    - custom_adata_path_celltype: Optional custom cell_type reference adata
    - auto_preprocess: Whether to preprocess raw data
    - clustering_resolution: Leiden clustering resolution
    - progress_callback: Optional GUI progress callback

    Returns:
    - Result summary dict
    """
    import os
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import torch
    from umap import UMAP

    print("\n" + "="*60)
    print("Starting Dual-Model Label Transfer with Shared UMAP Layout")
    print(f"Query file: {query_file}")
    print("="*60)

    # Check prerequisites
    check_prerequisites()
    output_folder = create_output_directory(output_folder)

    # Load query data
    try:
        query_adata = sc.read_h5ad(query_file)
        file_name = os.path.basename(query_file)
        base_name = os.path.splitext(file_name)[0]
        print(f"Loaded: {file_name} | Shape: {query_adata.shape}")
    except FileNotFoundError:
        raise FileNotFoundError(f"File not found: {query_file}")

    if progress_callback and hasattr(progress_callback, 'status_text'):
        progress_callback.status_text.text("Query dataset loaded. Checking preprocessing...")
        if hasattr(progress_callback, 'progress_bar'):
            progress_callback.progress_bar.progress(0.05)

    # Preprocessing
    if check_preprocessing(query_adata):
        print("Data is already preprocessed.")
        processed_adata = query_adata
    elif auto_preprocess:
        if not check_raw_data(query_adata):
            raise ValueError("Data is not raw counts and auto-preprocessing is enabled")
        if progress_callback and hasattr(progress_callback, 'status_text'):
            progress_callback.status_text.text("Preprocessing raw data...")
            if hasattr(progress_callback, 'progress_bar'):
                progress_callback.progress_bar.progress(0.1)
        processed_adata = preprocess_data(
            query_adata,
            resolution=clustering_resolution,
            output_dir=output_folder
        )
        preprocessed_file = os.path.join(output_folder, f"{base_name}_preprocessed.h5ad")
        processed_adata.write_h5ad(preprocessed_file)
        print(f"Preprocessed data saved: {preprocessed_file}")
    else:
        raise ValueError("Data not preprocessed and auto_preprocess=False")

    if 'orig.ident' not in processed_adata.obs:
        processed_adata.obs['orig.ident'] = base_name

    # Model configurations: cell_type MUST run first
    model_configs = {
        "cell_type": {
            "model_dir": custom_model_dir_celltype or './models/enhanced_reference_model_reanno_2ndround/',
            "adata_path": custom_adata_path_celltype,
            "cell_type_key": "reanno"
        },
        "lineage": {
            "model_dir": custom_model_dir_lineage or './models/enhanced_reference_model_lineage_2ndround/',
            "adata_path": custom_adata_path_lineage,
            "cell_type_key": "lineage"
        }
    }

    all_results = {}
    shared_umap_coords = None  # 存储共享的UMAP坐标

    for model_name, config in model_configs.items():
        print(f"\n--- Processing {model_name.upper()} Model ---")
        
        # Load model
        source_adata, enhanced_scpoli_model, cell_type_key = load_model(
            model_dir=config["model_dir"],
            adata_path=config["adata_path"],
            model_type=model_name
        )
        if source_adata is None or enhanced_scpoli_model is None:
            raise ValueError(f"Failed to load {model_name} model")

        # Align genes
        query_adata_extended = align_genes(processed_adata, source_adata)
        query_adata_extended.obs[cell_type_key] = 'Unknown'
        if 'orig.ident' not in query_adata_extended.obs:
            query_adata_extended.obs['orig.ident'] = base_name

        # Load query model
        try:
            scpoli_query = scPoli.load_query_data(query_adata_extended, enhanced_scpoli_model, labeled_indices=[])
        except Exception as e:
            print(f"❌ Failed to load query data for {model_name}: {e}")
            raise

        # Train
        print(f"Training {model_name} query model...")
        scpoli_query.train(n_epochs=50, pretraining_epochs=40, eta=10)

        # Predict
        results = scpoli_query.classify(query_adata_extended, scale_uncertainties=True)
        preds = results[cell_type_key]['preds']
        uncert = results[cell_type_key]['uncert']

        # Get latent
        latent = scpoli_query.get_latent(query_adata_extended, mean=True)
        adata_latent = sc.AnnData(latent.astype(np.float32))
        adata_latent.obs = query_adata_extended.obs.copy()
        adata_latent.obs[f'{cell_type_key}_pred'] = preds.tolist()
        adata_latent.obs[f'{cell_type_key}_uncert'] = uncert.tolist()

        # --- 关键修改：确保两个模型使用相同的UMAP布局 ---
        if model_name == "cell_type":
            # 对于cell_type模型，计算UMAP
            sc.pp.neighbors(adata_latent)
            sc.tl.umap(adata_latent)
            shared_umap_coords = adata_latent.obsm['X_umap'].copy()
            print(f"✅ Computed UMAP for cell_type: {shared_umap_coords.shape}")
        else:
            # 对于lineage模型，使用cell_type的UMAP坐标
            if shared_umap_coords is not None:
                adata_latent.obsm['X_umap'] = shared_umap_coords
                print(f"✅ Using shared UMAP layout for lineage: {shared_umap_coords.shape}")
            else:
                # 如果cell_type还没运行，计算默认UMAP
                sc.pp.neighbors(adata_latent)
                sc.tl.umap(adata_latent)
                print(f"⚠️  Computed default UMAP for lineage (cell_type not run yet)")

        # Store result
        all_results[model_name] = {
            'adata_latent': adata_latent,
            'cell_type_key': cell_type_key,
            'source_adata': source_adata,
            'scpoli_query': scpoli_query
        }

    # === Save predictions to processed_adata ===
    print("\n=== Saving dual-model results to AnnData ===")
    
    # 确保两个模型都使用相同的UMAP坐标
    if shared_umap_coords is not None:
        # 将共享的UMAP坐标保存到主数据对象
        processed_adata.obsm['X_umap_anno'] = shared_umap_coords
        
        for model_name, result in all_results.items():
            cell_type_key = result['cell_type_key']
            pred_col = f'{cell_type_key}_pred'
            uncert_col = f'{cell_type_key}_uncert'

            latent_obs = result['adata_latent'].obs
            unique_preds = latent_obs[pred_col].dropna().astype(str).unique().tolist()
            if 'Unknown' in unique_preds:
                unique_preds.remove('Unknown')
            all_preds = ['Unknown'] + sorted(unique_preds)

            pred_map = dict(zip(latent_obs.index, latent_obs[pred_col].astype(str)))
            uncert_map = dict(zip(latent_obs.index, latent_obs[uncert_col]))

            processed_adata.obs[pred_col] = [pred_map.get(idx, 'Unknown') for idx in processed_adata.obs.index]
            processed_adata.obs[uncert_col] = [uncert_map.get(idx, np.nan) for idx in processed_adata.obs.index]

            # 两个模型都使用相同的UMAP坐标
            processed_adata.obsm[f'X_umap_{cell_type_key}'] = shared_umap_coords
            print(f"✅ Saved {pred_col} with shared UMAP layout")

    # Clean dtypes
    for col in processed_adata.obs.columns:
        if col.endswith('_pred'):
            processed_adata.obs[col] = processed_adata.obs[col].astype('category')

    # Save
    annotated_file = os.path.join(output_folder, f"{base_name}_annotated.h5ad")
    processed_adata.write_h5ad(annotated_file)
    print(f"✅ Combined result saved: {annotated_file}")

    # === Plotting ===
    # 使用共享的UMAP坐标进行绘图
    lineage_color = {
        "Amniotic_ecto": "#1f77b4", "Notochord": "#aa40fc", "Endoderm": "#ff7f0e",
        "PGC": "#8c564b", "ExE_endo": "#279e68", "Primitive.streak": "#e377c2",
        "NMP": "#d62728", "TE_TrB": "#b5bd61", "epi": "#17becf",
        "hemogenic": "#aec7e8", "meso_Exe.meso": "#ffbb78", "neural_ecto": "#98df8a"
    }
    lineage_ordered = list(lineage_color.keys())
    reanno_ordered = [
        'TE', 'CTB_1','CTB_2', 'STB_1', 'STB_2', 'STB_3', 'EVT_1', 'EVT_2',
        'Epiblast_1','Epiblast_2','Epiblast_3','Ectoderm',
        'Amniontic.epi','Amniontic.ectoderm', 'PGC', 'Primitive.streak',
        'Neuromesodermal.progenitor', 'Neural.crest', 'Neural.ectoderm.forebrain',
        'Neural.ectoderm.hindbrain', 'Neural.ectoderm.midbrain','Spinal.cord',
        'Paraxial.mesoderm','Emergent.mesoderm','Pre-somatic.mesoderm','Somite',
        'Rostral.mesoderm', 'Lateral.plate.mesoderm_1', 'Lateral.plate.mesoderm_2',
        'Lateral.plate.mesoderm_3','Cardiac.mesoderm','Amniotic.mesoderm',
        'Exe.meso.progenitor','YS.mesoderm_1', 'YS.mesoderm_2', 'Hypoblast_1',
        'Hypoblast_2', 'AVE', 'VE', 'YS.endoderm', 'DE','Gut', 'Notochord',
        'Hemogenic.endothelial.progenitor','Endothelium','Erythroid',
        'Primitive.megakaryocyte','Myeloid.progenitor'
    ]

    for model_name, result in all_results.items():
        adata_latent = result['adata_latent']
        cell_type_key = result['cell_type_key']
        pred_col = f'{cell_type_key}_pred'

        # 确保使用共享的UMAP坐标
        adata_latent.obsm['X_umap'] = shared_umap_coords

        adata_latent.obs[pred_col] = adata_latent.obs[pred_col].astype('category')

        if cell_type_key == "lineage":
            cats = [c for c in lineage_ordered if c in adata_latent.obs[pred_col].cat.categories]
            cats += [c for c in adata_latent.obs[pred_col].cat.categories if c not in cats]
            adata_latent.obs[pred_col] = adata_latent.obs[pred_col].cat.reorder_categories(cats)
            palette = [lineage_color.get(c, "#cccccc") for c in cats]
            sc.pl.umap(adata_latent, color=pred_col, palette=palette, show=False, frameon=False)
        else:
            all_cats = reanno_ordered + [c for c in adata_latent.obs[pred_col].cat.categories if c not in reanno_ordered]
            adata_latent.obs[pred_col] = adata_latent.obs[pred_col].cat.set_categories(all_cats)
            sc.pl.umap(adata_latent, color=pred_col, show=False, frameon=False)

        plt.savefig(os.path.join(output_folder, f"{base_name}_{model_name}_prediction.pdf"), dpi=300, bbox_inches='tight')
        plt.close()

        sc.pl.umap(adata_latent, color=f'{cell_type_key}_uncert', cmap='magma', vmax=1, show=False, frameon=False)
        plt.savefig(os.path.join(output_folder, f"{base_name}_{model_name}_uncertainty.pdf"), dpi=300, bbox_inches='tight')
        plt.close()

    # Report
    report_file = os.path.join(output_folder, f"{base_name}_dual_model_report.txt")
    with open(report_file, 'w') as f:
        f.write("Dual-Model Label Transfer Report\n")
        f.write("="*50 + "\n")
        f.write(f"Input: {query_file}\n")
        f.write(f"Cells: {processed_adata.n_obs}\n\n")
        for name, res in all_results.items():
            cnt = processed_adata.obs[f"{res['cell_type_key']}_pred"].value_counts()
            f.write(f"[{name.upper()}]\nTop 5 predictions:\n")
            for i, (cat, n) in enumerate(cnt.head(5).items()):
                f.write(f"  {cat}: {n}\n")
            f.write("\n")

    if progress_callback and hasattr(progress_callback, 'status_text'):
        progress_callback.status_text.text("Task completed!")
        if hasattr(progress_callback, 'progress_bar'):
            progress_callback.progress_bar.progress(1.0)

    print("\n" + "="*60)
    print("✅ Dual-model label transfer with shared UMAP layout completed!")
    print(f"Results saved to: {output_folder}")
    print("="*60)

    return {
        'input_file': query_file,
        'output_directory': output_folder,
        'cell_count': processed_adata.n_obs,
        'annotated_file': annotated_file,
        'prediction_columns': ['lineage_pred', 'reanno_pred'],
        'umap_keys': ['X_umap_lineage', 'X_umap_reanno']
    }


def remove_sparsity(adata):
    """
    Remove sparsity from AnnData object
    """
    if issparse(adata.X):
        adata.X = adata.X.toarray()
    return adata


def main():
    """
    Main function - command line interface
    """
    import argparse
    
    parser = argparse.ArgumentParser(description='Single-cell RNA-seq Data Label Transfer Tool with Integrated Preprocessing')
    parser.add_argument('query_file', help='Query data file path (.h5ad)')
    parser.add_argument('output_folder', help='Output folder path')
    parser.add_argument('--model_dir_lineage', help='Custom model directory path for lineage')
    parser.add_argument('--model_dir_celltype', help='Custom model directory path for celltype')
    parser.add_argument('--adata_path_lineage', help='Custom reference data path for lineage')
    parser.add_argument('--adata_path_celltype', help='Custom reference data path for celltype')
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
        result = label_transfer_with_preprocessing_dual_model(
           query_file=args.query_file,
           output_folder=args.output_folder,
           custom_model_dir_lineage=args.model_dir_lineage,
           custom_model_dir_celltype=args.model_dir_celltype,
           custom_adata_path_lineage=args.adata_path_lineage,
           custom_adata_path_celltype=args.adata_path_celltype,
           auto_preprocess=auto_preprocess,
           clustering_resolution=args.resolution
        )
            
        print("\nTask completed! Result summary:")
        for key, value in result.items():
            print(f"  {key}: {value}")
            
    except Exception as e:
        print(f"\n❌ Task failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()

