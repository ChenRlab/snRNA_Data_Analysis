# -*- coding: utf-8 -*-
"""
Created on Tue Dec 30 11:07:17 2025

@author: maihuanzhuo
This script performs:
- Data loading from CellRanger outputs
- Quality control and filtering
- Removal of hemoglobin genes
- Normalization and HVG selection
- Dimensionality reduction and batch correction
- Clustering and manual cell-type annotation
"""

import numpy as np
import pandas as pd
import scanpy as sc
import leidenalg
import bbknn
import anndata
import re
import matplotlib.pyplot as plt
import os

os.chdir('E:/Omics_Data/snRNA')
sc.logging.print_header()
# Basic settings
sc.settings.verbosity = 3 
sc.settings.set_figure_params(dpi=80, 
                              dpi_save=300, 
                              format='pdf', 
                              vector_friendly=True, 
                              facecolor="white") 

#
samples = [d for d in next(os.walk("./Cellranger_result"))[1]]
print(samples)
# 
sc_list = []
for sample in samples:
    folder = os.path.join('./Cellranger_result', sample, 'filtered_feature_bc_matrix')
    print(sample)
    print(folder)
    print(os.listdir(folder))
    adata = sc.read_10x_mtx(folder, 
                            var_names='gene_symbols', 
                            cache=True)
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    adata.obs['sample'] = sample
    sc_list.append(adata)
# 
adata = sc.concat(sc_list, label='sample', keys=samples, index_unique='-')
adata 
print(adata.obs["sample"].value_counts())
adata.raw = adata
adata.X.data[:20]

# mitochondrial genes, "MT-" for human, "Mt-" for mouse
adata.var["mt"] = adata.var_names.str.startswith("mt-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl"))
# hemoglobin genes
adata.var["hb"] = adata.var_names.str.contains("^Hb[^(p)]")
# 
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)
sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    groupby="sample",
    rotation=45,
    multi_panel=True,
    show= False,
)
plt.savefig("violin_raw.png", dpi=300, bbox_inches="tight")
# 
sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt", show=False)
plt.savefig("scatter_raw.png", dpi=300, bbox_inches="tight")
# 
adata = adata[
    (adata.obs['n_genes_by_counts'] > 300) &
    (adata.obs['n_genes_by_counts'] < 5000) &
    (adata.obs['pct_counts_mt'] < 40),
    :,
].copy()
print(adata.shape)
adata 
## 
genes_to_remove = ['Hbb-bs', 'Hbb-bt', 'Hba-a1', 'Hba-a2']
adata = adata[:, ~adata.var_names.isin(genes_to_remove)].copy()
print(adata.shape)

# # 
# sc.pp.scrublet(adata, batch_key="sample", random_state = 42,
#                expected_doublet_rate = 0.05)
# # 
# adata.layers["counts"] = adata.X.copy()
# adata.layers

# Normalizing to median total counts
sc.pp.normalize_total(adata, target_sum=1e4)
# Logarithmize the data
sc.pp.log1p(adata)
# 
sc.pp.highly_variable_genes(
    adata,
    batch_key='sample',
    n_top_genes = 2000, 
    min_mean= 0.0125,   
    max_mean= 3,    
    min_disp= 0.5,
    flavor="seurat_v3"
)
# sc.pp.scale(adata, max_value=10)

# 
sc.tl.pca(adata, n_comps=30) 
sc.pl.pca_variance_ratio(adata, n_pcs=30, log=True) 
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata)

# 
sc.pl.umap(adata, color='sample', wspace=0.5)
sc.pl.umap(
    adata,
    color=["predicted_doublet", "doublet_score"],
    wspace=0.5,
    show=False,
)
plt.savefig("predicted_doublet.pdf", dpi=300, bbox_inches="tight")
######### 
adata = adata[~adata.obs["predicted_doublet"]].copy()
print(adata.shape) # print(adata.shape) 

####### 
sc.external.pp.bbknn(adata, batch_key="sample")  # running bbknn 1.3.6
sc.tl.umap(adata)
sc.pl.umap(adata, color='sample', title = 'raw data', wspace=0.5, show=False)
plt.savefig("bbknn_before.pdf", dpi=300, bbox_inches="tight")

# Clustering，leiden
sc.tl.leiden(adata, flavor="igraph", n_iterations=2, resolution = 0.5)
sc.pl.umap(adata, color=["leiden"], legend_loc="on data", show=False)
plt.savefig("leiden_res_0.5.png", dpi=300, bbox_inches="tight")

# Manual cell-type annotation
for res in [0.02, 0.05, 0.1, 0.5]:
    sc.tl.leiden(adata, key_added=f"leiden_res_{res:4.2f}", resolution=res)
# 
sc.pl.umap(adata, color=["leiden_res_0.02", "leiden_res_0.05", "leiden_res_0.10", "leiden_res_0.50"],
           legend_loc="on data", ncols=2, show=False)
plt.savefig("multi_leiden_res.png", dpi=300, bbox_inches="tight")
#
sc.pl.umap(adata, color=["leiden_res_0.50"], legend_loc="on data", show=False)
plt.savefig("leiden_res_0.50.pdf", dpi=300, bbox_inches="tight")
# 
marker_genes = {
    "Cardiomyocytes": ["Ttn", "Myh6", "Tnnt2", "Actn2", "Cox6a2"],
    "Fibroblasts": ["Col3a1", "Pdgfra"],
    "Endothelial cells": ["Kdr", "Cdh5"],
    "Smooth mucle cells": ["Myh11", "Mylk", "Myo1b", "Rgs5"],
    "Epicardial cells": ["Msln", "Wt1"],
    "Immune cells": ["Lgals3", "Ctsb", "Tyrobp", "Cd79a", "Cd3d", "Cd14"],
    "Adipocytes": ["Plin1", "Adipoq"], 
    "Schwann cells" : ["S100b", "Mpz"]
}
sc.pl.dotplot(adata, marker_genes, groupby="leiden", standard_scale="var")
sc.pl.dotplot(adata, marker_genes, groupby="leiden_res_0.50", standard_scale="var", show=False)
plt.savefig("dotplot_leiden_res_0.50.pdf", dpi=300, bbox_inches="tight")

# 
# sc.tl.rank_genes_groups(adata, groupby="leiden_res_1.0", groups=['25'], method="wilcoxon")
# sc.pl.rank_genes_groups_dotplot(adata, groupby="leiden_res_0.50", standard_scale="var", n_genes=5)
# sc.get.rank_genes_groups_df(adata, group="8").head(20)
# sc.get.rank_genes_groups_df(adata, group="4").head(20)

# 
new_cluster = {
    'Cardiomyocytes': ['0','6','7','8','9','15','17'],
    'Fibroblasts': ['2','10','11'],
    'Endothelial cells' : ['1','3','12','14'],
    "Smooth mucle cells": ['4','16'],
    "Epicardial cells": ['11'],
    "Immune cells": ['5','13'],
    "Adipocytes": ['19'], 
    "Schwann cells" : ['18'],
}
adata.obs['cell_type'] = 'Unknown'
for celltype, clusters in new_cluster.items():
    adata.obs.loc[adata.obs['leiden_res_0.50'].isin(clusters), 'cell_type'] = celltype
adata.obs['cell_type'].value_counts()
# 
sc.pl.umap(adata, color='cell_type', show=False)
plt.savefig("annotation_umap.png", dpi=300, bbox_inches="tight")

adata.write("./sce.all.h5ad")
# 








