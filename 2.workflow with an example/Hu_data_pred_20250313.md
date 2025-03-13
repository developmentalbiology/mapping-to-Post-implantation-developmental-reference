```python
## label transfer with scPoli
# In this notebook, we demonstrate a label transfer workflow using scPoli with an in-house embryo model dataset as the query. The model is available for download from Hugging Face.
```


```python
import os
import numpy as np
import anndata as ad
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import classification_report
from scarches.models.scpoli import scPoli
from scarches.dataset.trvae.data_handling import remove_sparsity

import warnings
warnings.filterwarnings('ignore')

```

    WARNING:root:In order to use the mouse gastrulation seqFISH datsets, please install squidpy (see https://github.com/scverse/squidpy).
    INFO:lightning_fabric.utilities.seed:[rank: 0] Global seed set to 0
    /home/liuxiaodongLab/fanxueying/miniconda3/envs/benchmarking/lib/python3.8/site-packages/flax/struct.py:132: FutureWarning: jax.tree_util.register_keypaths is deprecated, and will be removed in a future release. Please use `register_pytree_with_keys()` instead.
      jax.tree_util.register_keypaths(data_clz, keypaths)
    /home/liuxiaodongLab/fanxueying/miniconda3/envs/benchmarking/lib/python3.8/site-packages/flax/struct.py:132: FutureWarning: jax.tree_util.register_keypaths is deprecated, and will be removed in a future release. Please use `register_pytree_with_keys()` instead.
      jax.tree_util.register_keypaths(data_clz, keypaths)
    WARNING:root:mvTCR is not installed. To use mvTCR models, please install it first using "pip install mvtcr"
    WARNING:root:multigrate is not installed. To use multigrate models, please install it first using "pip install multigrate".



```python
sc.settings.set_figure_params(dpi=100, frameon=False)
sc.set_figure_params(dpi=100)
sc.set_figure_params(figsize=(3, 3))  
plt.rcParams['figure.dpi'] = 100
plt.rcParams['figure.figsize'] = (3, 3)
```


```python
## Data download
```


```python
# load query data
query_adata = sc.read_h5ad('/storage/liuxiaodongLab/fanxueying/mayanalysis/Garfield_run/Garfield_Bowen_data_20250102/bowen_data.h5ad')
query_adata
```




    AnnData object with n_obs × n_vars = 20494 × 38606
        obs: 'orig.ident', 'nCount_RNA', 'nFeature_RNA', 'sample', 'celltype'
        var: 'features'




```python
# Create a set of unique values
unique_values = set(query_adata.obs["orig.ident"])
# Ensure there is only one unique value
if len(unique_values) == 1:
    query_name = unique_values.pop()  # Extract the single value from the set
else:
    raise ValueError("Expected only one unique value in 'orig.ident', but found multiple or none.")
            
print(f"Query name: {query_name}")

# Preprocess query dataset
sc.settings.seed = 42
query_adata.layers["counts"] = query_adata.X.copy()
sc.pp.normalize_total(query_adata, target_sum=1e4)
sc.pp.log1p(query_adata)
query_adata.layers["logcounts"] = query_adata.X.copy()
sc.pp.highly_variable_genes(query_adata, n_top_genes=2000, flavor="cell_ranger")
sc.tl.pca(query_adata, n_comps=30, use_highly_variable=True)

from scipy.sparse import issparse
counts_matrix = query_adata.layers["counts"]
if issparse(counts_matrix):
    counts_matrix = counts_matrix.toarray()

query_adata = sc.AnnData(
    X=counts_matrix,
    obs=query_adata.obs.copy(),
    var=query_adata.var.copy(),
    layers={'counts': counts_matrix}
)
query_adata = remove_sparsity(query_adata)
```

    Query name: 0



```python
## Load lineage training model and perform label transfer of lineage
enhanced_model_dir='/storage/liuxiaodongLab/fanxueying/mayanalysis/scPoli/scpoli_optim/enhanced_reference_model_lineage/'
os.makedirs(enhanced_model_dir, exist_ok=True)  # Create directory for saving the model
adata_path = "/storage/liuxiaodongLab/fanxueying/mayanalysis/scPoli/scpoli_optim/enhanced_reference_model_lineage/combined_adata.h5ad"
source_adata = sc.read_h5ad(adata_path)
enhanced_scpoli_model = scPoli.load(enhanced_model_dir, adata=source_adata)

# load setting
condition_key = 'orig.ident'
cell_type_key = "final_lineage"
early_stopping_kwargs = {
    "early_stopping_metric": "val_prototype_loss",
    "mode": "min",
    "threshold": 0,
    "patience": 20,
    "reduce_lr": True,
    "lr_patience": 13,
    "lr_factor": 0.1,
}
```

    AnnData object with n_obs × n_vars = 40697 × 2000
        obs: 'orig.ident', 'nCount_RNA', 'nFeature_RNA', 'sample', 'stage', 'percent.mt', 'species', 'embryo', 'platform', 'ann_level_2', 'ann_level_3', 'ann_level_1', 'final_anno', 'final_lineage', 'conditions_combined'
        layers: 'counts'
    Embedding dictionary:
     	Num conditions: [15]
     	Embedding dim: [5]
    Encoder Architecture:
    	Input Layer in, out and cond: 2000 45 5
    	Mean/Var Layer in/out: 45 10
    Decoder Architecture:
    	First Layer in, out and cond:  10 45 5
    	Output Layer in/out:  45 2000 
    



```python
# Reorganize query dataset to match genes in the reference dataset
all_genes = source_adata.var_names
missing_genes = all_genes.difference(query_adata.var_names)
missing_data = np.zeros((query_adata.shape[0], len(missing_genes)))
query_adata_df = pd.DataFrame(query_adata.X, columns=query_adata.var_names, index=query_adata.obs_names)
missing_df = pd.DataFrame(missing_data, columns=missing_genes, index=query_adata.obs_names)
query_adata_combined_df = pd.concat([query_adata_df, missing_df], axis=1)[all_genes]
query_adata_extended = sc.AnnData(
    X=query_adata_combined_df.values,
    obs=query_adata.obs,
    var=pd.DataFrame(index=all_genes),
    layers={'counts': query_adata_combined_df.values}
)
query_adata_extended.var['features'] = query_adata.var.reindex(all_genes)['features']

```


```python
# Label transfer to query dataset
query_adata_extended.obs['final_lineage'] = 'Unknown'
query_adata_extended.obs['orig.ident'] = 'query'
```


```python
scpoli_query = scPoli.load_query_data(
    adata=query_adata_extended,
    reference_model=enhanced_scpoli_model,
    labeled_indices=[],
)
scpoli_query.train(
    n_epochs=50,
    pretraining_epochs=40,
    eta=10
)
query_adata_extended.X = query_adata_extended.X.astype(np.float32)

#Label transfer from reference to query
results_dict = scpoli_query.classify(query_adata_extended, scale_uncertainties=True)
```

    INFO:scarches.trainers.scpoli.trainer:GPU available: True, GPU used: True


    Embedding dictionary:
     	Num conditions: [16]
     	Embedding dim: [5]
    Encoder Architecture:
    	Input Layer in, out and cond: 2000 45 5
    	Mean/Var Layer in/out: 45 10
    Decoder Architecture:
    	First Layer in, out and cond:  10 45 5
    	Output Layer in/out:  45 2000 
    
    Warning: Labels in adata.obs[final_lineage] is not a subset of label-encoder!
    The missing labels are: {'Unknown'}
    Therefore integer value of those labels is set to -1
    Warning: Labels in adata.obs[final_lineage] is not a subset of label-encoder!
    The missing labels are: {'Unknown'}
    Therefore integer value of those labels is set to -1
    Warning: Labels in adata.obs[final_lineage] is not a subset of label-encoder!
    The missing labels are: {'Unknown'}
    Therefore integer value of those labels is set to -1
    Warning: Labels in adata.obs[final_lineage] is not a subset of label-encoder!
    The missing labels are: {'Unknown'}
    Therefore integer value of those labels is set to -1
    Initializing dataloaders
    Starting training
     |████████████████----| 80.0%  - val_loss: 1321.48 - val_cvae_loss: 1321.48
    Initializing unlabeled prototypes with Leiden with an unknown number of  clusters.
    Clustering succesful. Found 30 clusters.
     |████████████████████| 100.0%  - val_loss: 1242.25 - val_cvae_loss: 1242.25 - val_prototype_loss:    0.00 - val_unlabeled_loss:    0.18



```python
#get latent representation of reference data
scpoli_query.model.eval()
data_latent_source = scpoli_query.get_latent(
    source_adata,
    mean=True
)

adata_latent_source = sc.AnnData(data_latent_source)
adata_latent_source.obs = source_adata.obs.copy()

#get latent representation of query data
data_latent= scpoli_query.get_latent(
    query_adata_extended,
    mean=True
)

adata_latent = sc.AnnData(data_latent)
adata_latent.obs = query_adata_extended.obs.copy()

#get label annotations
adata_latent.obs['final_lineage_pred'] = results_dict['final_lineage']['preds'].tolist()
adata_latent.obs['final_lineage_uncert'] = results_dict['final_lineage']['uncert'].tolist()
adata_latent.obs['classifier_outcome'] = (
    adata_latent.obs['final_lineage_pred'] == adata_latent.obs['final_lineage']
)

#get prototypes
labeled_prototypes = scpoli_query.get_prototypes_info()
labeled_prototypes.obs['study'] = 'labeled prototype'
unlabeled_prototypes = scpoli_query.get_prototypes_info(prototype_set='unlabeled')
unlabeled_prototypes.obs['study'] = 'unlabeled prototype'

#join adatas
adata_latent_full = adata_latent_source.concatenate(
    [adata_latent, labeled_prototypes, unlabeled_prototypes],
    batch_key='query'
)
adata_latent_full.obs['final_lineage_pred'][adata_latent_full.obs['query'].isin(['0'])] = np.nan
sc.pp.neighbors(adata_latent_full, n_neighbors=15)
sc.tl.umap(adata_latent_full)
```


```python
#get adata without prototypes
adata_no_prototypes = adata_latent_full[adata_latent_full.obs['query'].isin(['0', '1'])]
```


```python
sc.pl.umap(
    adata_no_prototypes,
    color='final_lineage_pred',
    show=False,
    frameon=False,
)
sc.pl.umap(
    adata_no_prototypes,
    color='orig.ident',
    show=False,
    frameon=False,
)
sc.pl.umap(
    adata_no_prototypes,
    color='final_lineage_uncert',
    show=False,
    frameon=False,
    cmap='magma',
    vmax=1
)
```




    <AxesSubplot: title={'center': 'final_lineage_uncert'}, xlabel='UMAP1', ylabel='UMAP2'>




    
![png](Hu_data_pred_20250313_files/Hu_data_pred_20250313_12_1.png)
    



    
![png](Hu_data_pred_20250313_files/Hu_data_pred_20250313_12_2.png)
    



    
![png](Hu_data_pred_20250313_files/Hu_data_pred_20250313_12_3.png)
    



```python
sc.pp.neighbors(adata_latent)
sc.tl.leiden(adata_latent)
sc.tl.umap(adata_latent)
sc.pl.umap(
    adata_latent,
    color='final_lineage_pred',
    show=False,
    frameon=False,
)
```




    <AxesSubplot: title={'center': 'final_lineage_pred'}, xlabel='UMAP1', ylabel='UMAP2'>




    
![png](Hu_data_pred_20250313_files/Hu_data_pred_20250313_13_1.png)
    



```python
# Ensure indices match between the two AnnData objects
matching_indices = adata_latent.obs.index.intersection(query_adata.obs.index)

# Select relevant columns and subset based on matching indices
query_adata.obs.loc[matching_indices, ["final_lineage_pred", "final_lineage_uncert"]] = \
    adata_latent.obs.loc[matching_indices, ["final_lineage_pred", "final_lineage_uncert"]]
```


```python
query_adata.obsm['X_umap_lineage'] = adata_latent.obsm['X_umap']
```


```python
## Load anno training model and perform label transfer of lineage
enhanced_model_dir='/storage/liuxiaodongLab/fanxueying/mayanalysis/scPoli/scpoli_optim/enhanced_reference_model/'
os.makedirs(enhanced_model_dir, exist_ok=True)  # Create directory for saving the model
adata_path = "/storage/liuxiaodongLab/fanxueying/mayanalysis/scPoli/scpoli_optim/enhanced_reference_model/combined_adata.h5ad"
source_adata = sc.read_h5ad(adata_path)
enhanced_scpoli_model = scPoli.load(enhanced_model_dir, adata=source_adata)

# load setting
condition_key = 'orig.ident'
cell_type_key = "final_anno"
early_stopping_kwargs = {
    "early_stopping_metric": "val_prototype_loss",
    "mode": "min",
    "threshold": 0,
    "patience": 20,
    "reduce_lr": True,
    "lr_patience": 13,
    "lr_factor": 0.1,
}
```

    AnnData object with n_obs × n_vars = 36502 × 2000
        obs: 'orig.ident', 'nCount_RNA', 'nFeature_RNA', 'sample', 'stage', 'percent.mt', 'species', 'embryo', 'platform', 'ann_level_2', 'ann_level_3', 'ann_level_1', 'final_anno', 'final_lineage', 'conditions_combined'
        layers: 'counts'
    Embedding dictionary:
     	Num conditions: [15]
     	Embedding dim: [5]
    Encoder Architecture:
    	Input Layer in, out and cond: 2000 45 5
    	Mean/Var Layer in/out: 45 10
    Decoder Architecture:
    	First Layer in, out and cond:  10 45 5
    	Output Layer in/out:  45 2000 
    



```python
# Reorganize query dataset to match genes in the reference dataset
all_genes = source_adata.var_names
missing_genes = all_genes.difference(query_adata.var_names)
missing_data = np.zeros((query_adata.shape[0], len(missing_genes)))
query_adata_df = pd.DataFrame(query_adata.X, columns=query_adata.var_names, index=query_adata.obs_names)
missing_df = pd.DataFrame(missing_data, columns=missing_genes, index=query_adata.obs_names)
query_adata_combined_df = pd.concat([query_adata_df, missing_df], axis=1)[all_genes]
query_adata_extended = sc.AnnData(
    X=query_adata_combined_df.values,
    obs=query_adata.obs,
    var=pd.DataFrame(index=all_genes),
    layers={'counts': query_adata_combined_df.values}
)
query_adata_extended.var['features'] = query_adata.var.reindex(all_genes)['features']


```


```python
# Label transfer to query dataset
query_adata_extended.obs['final_anno'] = 'Unknown'
query_adata_extended.obs['orig.ident'] = 'query'
scpoli_query = scPoli.load_query_data(
    adata=query_adata_extended,
    reference_model=enhanced_scpoli_model,
    labeled_indices=[],
)
scpoli_query.train(
    n_epochs=50,
    pretraining_epochs=40,
    eta=10
)
query_adata_extended.X = query_adata_extended.X.astype(np.float32)

#Label transfer from reference to query
results_dict = scpoli_query.classify(query_adata_extended, scale_uncertainties=True)
```

    INFO:scarches.trainers.scpoli.trainer:GPU available: True, GPU used: True


    Embedding dictionary:
     	Num conditions: [16]
     	Embedding dim: [5]
    Encoder Architecture:
    	Input Layer in, out and cond: 2000 45 5
    	Mean/Var Layer in/out: 45 10
    Decoder Architecture:
    	First Layer in, out and cond:  10 45 5
    	Output Layer in/out:  45 2000 
    
    Warning: Labels in adata.obs[final_anno] is not a subset of label-encoder!
    The missing labels are: {'Unknown'}
    Therefore integer value of those labels is set to -1
    Warning: Labels in adata.obs[final_anno] is not a subset of label-encoder!
    The missing labels are: {'Unknown'}
    Therefore integer value of those labels is set to -1
    Warning: Labels in adata.obs[final_anno] is not a subset of label-encoder!
    The missing labels are: {'Unknown'}
    Therefore integer value of those labels is set to -1
    Warning: Labels in adata.obs[final_anno] is not a subset of label-encoder!
    The missing labels are: {'Unknown'}
    Therefore integer value of those labels is set to -1
    Initializing dataloaders
    Starting training
     |████████████████----| 80.0%  - val_loss: 1328.89 - val_cvae_loss: 1328.89
    Initializing unlabeled prototypes with Leiden with an unknown number of  clusters.
    Clustering succesful. Found 31 clusters.
     |████████████████████| 100.0%  - val_loss: 1247.94 - val_cvae_loss: 1247.94 - val_prototype_loss:    0.00 - val_unlabeled_loss:    0.23



```python
#get latent representation of reference data
scpoli_query.model.eval()
data_latent_source = scpoli_query.get_latent(
    source_adata,
    mean=True
)

adata_latent_source = sc.AnnData(data_latent_source)
adata_latent_source.obs = source_adata.obs.copy()

#get latent representation of query data
data_latent= scpoli_query.get_latent(
    query_adata_extended,
    mean=True
)

adata_latent = sc.AnnData(data_latent)
adata_latent.obs = query_adata_extended.obs.copy()


#get label annotations
adata_latent.obs['final_anno_pred'] = results_dict['final_anno']['preds'].tolist()
adata_latent.obs['final_anno_uncert'] = results_dict['final_anno']['uncert'].tolist()
adata_latent.obs['classifier_outcome'] = (
    adata_latent.obs['final_anno_pred'] == adata_latent.obs['final_anno']
)

#get prototypes
labeled_prototypes = scpoli_query.get_prototypes_info()
labeled_prototypes.obs['study'] = 'labeled prototype'
unlabeled_prototypes = scpoli_query.get_prototypes_info(prototype_set='unlabeled')
unlabeled_prototypes.obs['study'] = 'unlabeled prototype'

#join adatas
adata_latent_full = adata_latent_source.concatenate(
    [adata_latent, labeled_prototypes, unlabeled_prototypes],
    batch_key='query'
)
adata_latent_full.obs['final_anno_pred'][adata_latent_full.obs['query'].isin(['0'])] = np.nan
sc.pp.neighbors(adata_latent_full, n_neighbors=15)
sc.tl.umap(adata_latent_full)
```


```python
#get adata without prototypes
adata_no_prototypes = adata_latent_full[adata_latent_full.obs['query'].isin(['0', '1'])]
```


```python
sc.pl.umap(
    adata_no_prototypes,
    color='final_anno_pred',
    show=False,
    frameon=False,
)
sc.pl.umap(
    adata_no_prototypes,
    color='orig.ident',
    show=False,
    frameon=False,
)
sc.pl.umap(
    adata_no_prototypes,
    color='final_anno_uncert',
    show=False,
    frameon=False,
    cmap='magma',
    vmax=1
)
```




    <AxesSubplot: title={'center': 'final_anno_uncert'}, xlabel='UMAP1', ylabel='UMAP2'>




    
![png](Hu_data_pred_20250313_files/Hu_data_pred_20250313_21_1.png)
    



    
![png](Hu_data_pred_20250313_files/Hu_data_pred_20250313_21_2.png)
    



    
![png](Hu_data_pred_20250313_files/Hu_data_pred_20250313_21_3.png)
    



```python
sc.pp.neighbors(adata_latent)
sc.tl.leiden(adata_latent)
sc.tl.umap(adata_latent)
sc.pl.umap(
    adata_latent,
    color='final_anno_pred',
    show=False,
    frameon=False,
)
```




    <AxesSubplot: title={'center': 'final_anno_pred'}, xlabel='UMAP1', ylabel='UMAP2'>




    
![png](Hu_data_pred_20250313_files/Hu_data_pred_20250313_22_1.png)
    



```python
adata_latent
```




    AnnData object with n_obs × n_vars = 20494 × 10
        obs: 'orig.ident', 'nCount_RNA', 'nFeature_RNA', 'sample', 'celltype', 'final_lineage_pred', 'final_lineage_uncert', 'final_anno', 'conditions_combined', 'final_anno_pred', 'final_anno_uncert', 'classifier_outcome', 'leiden'
        uns: 'neighbors', 'leiden', 'umap', 'final_anno_pred_colors'
        obsm: 'X_umap'
        obsp: 'distances', 'connectivities'




```python
# Ensure indices match between the two AnnData objects
matching_indices = adata_latent.obs.index.intersection(query_adata.obs.index)

# Select relevant columns and subset based on matching indices
query_adata.obs.loc[matching_indices, ["final_anno_pred", "final_anno_uncert"]] = \
    adata_latent.obs.loc[matching_indices, ["final_anno_pred", "final_anno_uncert"]]
```


```python
query_adata.obsm['X_umap_anno'] = adata_latent.obsm['X_umap']
```


```python
query_adata
```




    AnnData object with n_obs × n_vars = 20494 × 38606
        obs: 'orig.ident', 'nCount_RNA', 'nFeature_RNA', 'sample', 'celltype', 'final_lineage_pred', 'final_lineage_uncert', 'final_anno_pred', 'final_anno_uncert'
        var: 'features', 'highly_variable', 'means', 'dispersions', 'dispersions_norm'
        obsm: 'X_umap_lineage', 'X_umap_anno'
        layers: 'counts'




```python
query_adata.obsm['X_umap'] = query_adata.obsm['X_umap_anno']
sc.pl.umap(query_adata, color=['final_lineage_pred'])
sc.pl.umap(query_adata, color=['final_anno_pred'])
sc.pl.umap(query_adata, color=['celltype'])
```


    
![png](Hu_data_pred_20250313_files/Hu_data_pred_20250313_27_0.png)
    



    
![png](Hu_data_pred_20250313_files/Hu_data_pred_20250313_27_1.png)
    



    
![png](Hu_data_pred_20250313_files/Hu_data_pred_20250313_27_2.png)
    



```python
# Define the list of genes you want to plot
genes_to_plot = ["POU5F1","NANOG",#epiblast
                 "SOX2", "TTYH1", #Neural_ectoderm
                 "GATA3","TFAP2A",
                 "TBXT",  "CDX1",
                 "PDGFRA","APOA2","FOXA2","NANOS3",
                 "PECAM1","HBZ","PTPRC",
                "GABRP","HEY1","COL6A1","COL6A2",
                "GATA6", "GATA4", "SOX17", "LHX1", "HHEX","LEFTY1", "LEFTY2", "EPCAM", "CLDN6"
         ]  # Replace "GENE2" and "GENE3" with the actual gene names

# Plot the UMAP embedding for the specified genes
sc.pl.embedding(
    query_adata,
    basis='X_umap',
    color=genes_to_plot,
    use_raw=False,
    cmap=sns.cubehelix_palette(dark=0, light=.9, as_cmap=True)
)
```


    
![png](Hu_data_pred_20250313_files/Hu_data_pred_20250313_28_0.png)
    



```python

```
