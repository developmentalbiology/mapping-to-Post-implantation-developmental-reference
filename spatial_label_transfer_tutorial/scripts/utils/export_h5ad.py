#!/usr/bin/env python3
"""
Export H5AD File Data for R/Seurat Conversion

This script exports data from an AnnData h5ad file to formats compatible with R/Seurat:
- Expression matrix (MTX format)
- Cell barcodes and gene features
- Metadata (obs and var)
- Highly variable genes
- UMAP coordinates

Author: Spatial Annotation Workflow
Version: 2.0
"""

import argparse
import os
import sys
import anndata
import scanpy as sc
import numpy as np
import pandas as pd
from scipy.io import mmwrite


def export_h5ad_data(h5ad_path, output_dir, hvg_count=2000):
    """
    Export h5ad file data to multiple formats for downstream analysis
    
    Parameters:
    -----------
    h5ad_path : str
        Path to the input h5ad file
    output_dir : str
        Directory to save output files
    hvg_count : int
        Number of highly variable genes to export (default: 2000)
    """
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"Loading h5ad file: {h5ad_path}")
    adata = anndata.read_h5ad(h5ad_path)
    
    print(f"Dataset shape: {adata.shape[0]} cells x {adata.shape[1]} genes")
    
    # Export highly variable genes
    if 'highly_variable' in adata.var.columns:
        hvg = adata.var[adata.var['highly_variable']].index.tolist()
        hvg_df = pd.DataFrame(hvg, columns=['highly_variable_genes'])
        hvg_output = os.path.join(output_dir, f'highly_variable_genes_{hvg_count}.csv')
        hvg_df.to_csv(hvg_output, index=False)
        print(f"Exported {len(hvg)} highly variable genes to {hvg_output}")
    else:
        print("Warning: No highly variable genes found in the dataset")
    
    # Export expression matrix in MTX format
    print("Exporting expression matrix...")
    adata_matrix = adata.X
    adata_transposed = np.transpose(adata_matrix)
    
    matrix_output = os.path.join(output_dir, "matrix.mtx")
    mmwrite(matrix_output, adata_transposed, field="integer")
    print(f"Expression matrix exported to {matrix_output}")
    
    # Export cell barcodes
    cell_names = adata.obs_names.tolist()
    barcodes_output = os.path.join(output_dir, "barcodes.tsv")
    with open(barcodes_output, "wt") as f:
        f.write("\n".join(cell_names))
    print(f"Cell barcodes exported to {barcodes_output}")
    
    # Export gene features
    gene_names = adata.var_names.tolist()
    features_output = os.path.join(output_dir, "features.tsv")
    with open(features_output, "wt") as f:
        f.write("\n".join(gene_names))
    print(f"Gene features exported to {features_output}")
    
    # Export metadata
    obs_output = os.path.join(output_dir, 'obs_metadata.csv')
    adata.obs.to_csv(obs_output)
    print(f"Cell metadata exported to {obs_output}")
    
    var_output = os.path.join(output_dir, 'var_metadata.csv')
    adata.var.to_csv(var_output)
    print(f"Gene metadata exported to {var_output}")
    
    # Export UMAP coordinates if available
    if 'X_umap' in adata.obsm.keys():
        umap_coordinates = adata.obsm['X_umap']
        umap_df = pd.DataFrame(umap_coordinates, 
                              index=adata.obs.index, 
                              columns=['UMAP1', 'UMAP2'])
        umap_output = os.path.join(output_dir, 'umap_coordinates.csv')
        umap_df.to_csv(umap_output)
        print(f"UMAP coordinates exported to {umap_output}")
    else:
        print("Warning: No UMAP coordinates found in the dataset")
    
    print(f"All files exported successfully to {output_dir}")


def main():
    parser = argparse.ArgumentParser(description='Export H5AD file data for R/Seurat conversion')
    parser.add_argument('--input', '-i', required=True, 
                       help='Path to input h5ad file')
    parser.add_argument('--output', '-o', required=True,
                       help='Output directory for exported files')
    parser.add_argument('--hvg-count', type=int, default=2000,
                       help='Number of highly variable genes (default: 2000)')
    
    args = parser.parse_args()
    
    # Validate input file
    if not os.path.exists(args.input):
        print(f"Error: Input file {args.input} does not exist")
        sys.exit(1)
    
    # Export data
    export_h5ad_data(args.input, args.output, args.hvg_count)


if __name__ == "__main__":
    main()

