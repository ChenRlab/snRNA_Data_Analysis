# -*- coding: utf-8 -*-
"""
Created on Fri Dec 30 19:25:01 2025

@author: maihuanzhuo
"""

import scanpy as sc
import matplotlib.pyplot as plt
import decoupler as dc
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
# 
cov_df = cov_df.set_index("sample")
# 
collectri = dc.op.collectri(organism="mouse")
collectri
# 
score = dc.pp.get_obsm(adata=adata, key="score_ulm")
score
# score.write_h5ad("decoupler_TF_score.h5ad")
# score = sc.read_h5ad("decoupler_TF_score.h5ad")
score
# 
tf = "Trp53"
sc.pl.umap(score, color=[tf], cmap="RdBu_r", vcenter=0, title=[f"{tf} score"], show=False)
plt.savefig("Trp53_score_umap.png", dpi=300, bbox_inches="tight")
# 
ax = sc.pl.violin(score, keys=[tf], groupby="group", rotation=45, ylabel=f"{tf} score", show=False)
for label in ax.get_xticklabels():
    label.set_horizontalalignment("right")
    label.set_verticalalignment("top")
plt.savefig("Trp53_score_violin.png", dpi=300, bbox_inches="tight")
# 
################## progeny Pathway Scoring
progeny = dc.op.progeny(organism="mouse")
progeny

# 
dc.mt.ulm(data=adata, net=progeny)
progeny_score = dc.pp.get_obsm(adata=adata, key="score_ulm")
progeny_score
# 
adata_cm = adata[adata.obs["cell_type"] == "Cardiomyocytes"].copy()
progeny_score_cm = dc.pp.get_obsm(adata=adata_cm, key="score_ulm")
progeny_score_cm
# 
sc.pl.umap(progeny_score_cm, color=["p53"], cmap="RdBu_r", vcenter=0)
sc.pl.violin(progeny_score_cm, keys=["p53"], groupby="group", rotation=90)
# 
mp = sc.pl.matrixplot(
    adata=progeny_score_cm,
    var_names=progeny_score_cm.var_names,
    groupby="group",
    standard_scale="var", # z-score
    cmap="RdBu_r",
    show=False,
)
ax = mp["mainplot_ax"]
for label in ax.get_xticklabels():
    label.set_rotation(45)
    label.set_rotation_mode("anchor") 
    label.set_ha("right")
    label.set_va("top")
plt.savefig("progeny_heatmap_all_group.pdf", dpi=600, bbox_inches="tight")
### 
X = progeny_score_cm.to_df()
# 
groups = progeny_score_cm.obs["group"]
df = X.copy()
df["group"] = progeny_score_cm.obs["group"].values
df.to_csv("progeny_score_all.csv")
# 
group_score_mean = X.groupby(groups).mean()
group_score_mean
# calculate z-score
from scipy.stats import zscore
group_score_z = group_score_mean.apply(zscore, axis=0, nan_policy="omit")
group_score_z
group_score_z.to_csv("progeny_score_cm_score_z.csv")