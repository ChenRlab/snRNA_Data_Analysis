# -*- coding: utf-8 -*-
"""
Created on Fri Dec 30 09:25:01 2025

@author: maihuanzhuo
"""

import scanpy as sc
import anndata
import pandas as pd
from sccoda.util import comp_ana as mod
from sccoda.util import cell_composition_data as dat
from sccoda.util import data_visualization as viz
import matplotlib.pyplot as plt
import pickle as pkl
import arviz as az 
import importlib
import warnings
warnings.filterwarnings("ignore")
import os
os.chdir('E:/Omics_Data/snRNA')
# 
adata = sc.read_h5ad("./ann.h5ad")
adata
adata.obs['cell_type'].value_counts()
# 
cov_df = pd.DataFrame({
    "sample": adata.obs["sample"].unique(),
})
cov_df["group"] = (
    cov_df["sample"].astype(str).str[0]
    .map({"C": "ang", "D": "BPF_ang"})
    .astype("category")
)
# 
cov_df = cov_df.set_index("sample")
# 
cov_df['group'] = pd.Categorical(cov_df['group'], 
                                categories=["ang", "BPF_ang"],
                                ordered=True)
cov_df
# scCODA object
scCODA_obj = dat.from_scanpy(
    adata,
    cell_type_identifier="cell_type",
    sample_identifier="sample",
    covariate_df=cov_df
)
scCODA_obj
# 
data_compare = scCODA_obj[scCODA_obj.obs["group"].isin(["ang", "BPF_ang"])]
print(data_compare.obs)
## 
viz.boxplots(data_compare, feature_name="group")
plt.savefig("boxplot_group.png", dpi=300, bbox_inches="tight")
plt.show()

# Stacked barplot for each sample
viz.stacked_barplot(data_compare, feature_name="samples")
plt.savefig("stacked_barplot_sample.png", dpi=300, bbox_inches="tight")
plt.show()

# Stacked barplot for the levels of "Condition"
viz.stacked_barplot(data_compare, feature_name="group")
plt.savefig("stacked_barplot_group.png", dpi=300, bbox_inches="tight")
plt.show()

viz.rel_abundance_dispersion_plot(
    data=data_compare,
    abundant_threshold=0.9
)
plt.savefig("rel_abundance_dispersion_plot.png", dpi=300, bbox_inches="tight")
plt.show()

## 
model_compare = mod.CompositionalAnalysis(data_compare, formula="group", reference_cell_type="Endothelial cells") 

# Run MCMC
compare_results = model_compare.sample_hmc(num_results=20000) 
# 
compare_results.summary()
# Intercepts:
#                     Final Parameter  Expected Sample
# Cell Type                                           
# Cardiomyocytes                4.792      6139.922369
# Endothelial cells             4.410      4190.466644
# Fibroblasts                   3.659      1977.457836
# Smooth mucle cells            2.718       771.678010
# Immune cells                  2.623       701.743096
# Epicardial cells              1.206       170.130770
# Schwann cells                 0.243        64.946729
# Adipocytes                   -0.223        40.754545

print(compare_results.credible_effects())
# 调整假阳性率
compare_results.set_fdr(est_fdr=0.05)
compare_results.summary()
# Effects:
#                                      Final Parameter  ...  log2-fold change
# Covariate        Cell Type                            ...                  
# group[T.BPF_ang] Cardiomyocytes            -0.306400  ...         -0.261030
#                  Endothelial cells          0.000000  ...          0.181012
#                  Fibroblasts                0.000000  ...          0.181012
#                  Smooth mucle cells         0.000000  ...          0.181012
#                  Immune cells               0.000000  ...          0.181012
#                  Epicardial cells          -0.191783  ...         -0.095672
#                  Schwann cells              0.019043  ...          0.208485
#                  Adipocytes                -0.241490  ...         -0.167384

# plot for R
effects_df = compare_results.effect_df.copy()
effects_df.to_csv("scCODA_effects.csv", index=True)

### 
compare_results.summary_extended(hdi_prob=0.9)
# 
az.plot_trace(
    compare_results,
    divergences=False,
    var_names=["alpha", "beta"],
    coords={"cell_type": compare_results.posterior.coords["cell_type_nb"]},
)
plt.savefig("plot_trace.png", dpi=300, bbox_inches="tight")
plt.show()

# Run scCODA with each cell type as the reference
cell_types = data_compare.var.index
results_cycle = pd.DataFrame(index=cell_types, columns=["times_credible"]).fillna(0)

for ct in cell_types:
    print(f"Reference: {ct}")
    # Run inference
    model_temp = mod.CompositionalAnalysis(data_compare, formula="group", reference_cell_type=ct)
    temp_results = model_temp.sample_hmc(num_results=20000)
    # Select credible effects
    cred_eff = temp_results.credible_effects()
    cred_eff.index = cred_eff.index.droplevel(level=0)
    # add up credible effects
    results_cycle["times_credible"] += cred_eff.astype("int")
    
# save
path = "BPF_ang_vs_ang_1"
compare_results.save(path)
# loading
with open(path, "rb") as f:
    compare_results = pkl.load(f)
compare_results.summary()


