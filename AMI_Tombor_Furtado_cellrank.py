#!/usr/bin/env python
# coding: utf-8

# # scVelo Analysis on Myo combined dataset
# # Author: Lukas Tombor
# # Date: 22-03-14
# 
# This is a continued analysis of MNN-corrected Young only samples

# In[1]:


import scvelo as scv
import scanpy
import os
import pandas as pd
import glob
import numpy as np
import scanpy as sc
import scanorama
import re
from collections.abc import Iterable
import cellrank as cr
from matplotlib.pyplot import rc_context


scv.set_figure_params()


# In[2]:


# Import dataset - h5ad

adata = scv.read(filename= "/media/tombor/Helios_scStorage/Lukas/Analysis/EHT/Combined_Seurats/scanpy_mnn_corrected_Young.h5ad")


# In[3]:


# Import GFP.positives
anno=pd.read_csv("/media/tombor/Helios_scStorage/Lukas/Analysis/EHT/Combined_Seurats/GFP_anno.csv", sep = ",", index_col=0)
print(anno)
print(adata.obs.index)

print(adata.obs.index == anno.index)

adata.obs["GFP"] = anno["GFP"]
print(adata.obs)


# In[4]:


# Basic preprocessing
scv.pl.proportions(adata)
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
sc.tl.pca(adata)
sc.pp.neighbors(adata, n_pcs=30, n_neighbors=30)
scv.pp.moments(adata, n_pcs=None, n_neighbors=None)


# In[5]:


# Import Cell Rank

import cellrank as cr
from matplotlib.pyplot import rc_context


# In[6]:


# Build subset

## Identify endothelial cells

sc.tl.leiden(adata, resolution = 0.1)


# In[7]:


# rc_context is used for the figure size, in this case 4x4
with rc_context({'figure.figsize': (4, 4)}):
    sc.pl.umap(adata, color=['Cdh5', 'Pecam1', 'leiden', "GFP"], save = "Innitial_clustering.pdf")
    
ax = sc.pl.dotplot(adata, ["Acta2", "Postn", "Col1a1", "Col1a2", "Col3a1", "Pdgfra", "Col4a1", "Cd74", "Lyz2", "Itgam", "Fcgr4", "Ccr2", "Igf1", "Timd4", "Siglec1", "Cd14", "Cxcr2", "Lgals3", "Napsa", "Cd209a", "Itgb7", "Clec10a", "Flt3", "Vwf", "Eng", "Emcn", "Fabp4", "Pecam1", "Cdh5", "Tie1", "Lyve1", "Pdpn", "Prox1", "Ccl5", "Nkg7", "Ptprc", "Klre1", "Cd4", "Cd8a", "Ctla4", "Icos", "Cd69", "Il2ra", "Itgae", "Cd3e", "Cd44", "Ccr7", "Lat", "Cd79a", "Cd79b", "Pax5", "Plp1", "Kcna1", "Cd59a", "Rgs5", "Pdgfrb", "Tagln", "Des", "Krt19", "Krt8", "Upk3b", "Tnnt2", "Nkain4"], 
                          groupby='leiden', save = "Innitial_clustering_markers.pdf")


# In[8]:


# rough characterisation

new_cluster_names = ["Fibroblasts", "Myeloid cells", "Endothelial cells", "B-cells", "T-cells", "Mural cells", "Epithelial cells"]

adata.rename_categories("leiden", new_cluster_names)

#sc.pl.umap(adata, color='leiden', legend_loc='on data', title='', frameon=False, save='.pdf')

# Subset only 2 populations

adata3 = adata[adata.obs["leiden"].isin(["Myeloid cells", "Endothelial cells"])]


# In[9]:


sc.tl.pca(adata3)
sc.pp.neighbors(adata3, n_pcs=15, n_neighbors=15)
scv.pp.moments(adata3, n_pcs=None, n_neighbors=None)
sc.tl.umap(adata3)


# In[10]:


# run scVelo

scv.tl.recover_dynamics(adata3, n_jobs=30)


# In[11]:


# built velocity graph

scv.tl.velocity(adata3, mode="dynamical")
scv.tl.velocity_graph(adata3)


# In[12]:


scv.pl.velocity_embedding_stream(
    adata3, basis="umap", legend_fontsize=12, title="", smooth=0.8, min_mass=3, color="leiden", save = ".pdf")


# In[14]:


# Annotate Timepoints per integer

adata3.rename_categories("Timepoint", [0,1,3,5,7,14,28])

pd.set_option('display.max_rows', None, 'display.max_columns', None)

adata3.obs.groupby(["Timepoint","Sample"]).size()

# Subsample

#sc.pp.subsample(adata3, fraction=0.5, random_state=0)
adata3


# In[15]:


# Calculate WOTKernel

from cellrank.external.kernels import WOTKernel
from cellrank.tl.kernels import ConnectivityKernel
from cellrank.tl.estimators import GPCCA

wk = WOTKernel(adata3, time_key="Timepoint")


# In[16]:


wk.compute_initial_growth_rates(organism="mouse", key_added="growth_rate_init")
scv.pl.scatter(
    adata3, c="growth_rate_init", legend_loc="right", s=10, save= "Growth_rate.pdf"
)


# In[17]:


wk.compute_transition_matrix(
    growth_iters=3, growth_rate_key="growth_rate_init", last_time_point="connectivities"
)


# In[17]:


wk.plot_random_walks(
    n_sims=300,
    max_iter=200,
    start_ixs={"Timepoint": 0},
    stop_ixs={"Timepoint": 1},
    c="Timepoint",
    legend_loc="none",
    title = "Day 0 to day 1",
    linealpha=0.5,
    dpi=600, save = "Day0to1.pdf"
)

wk = WOTKernel(adata3, time_key="Timepoint")
wk.compute_initial_growth_rates(organism="mouse", key_added="growth_rate_init")
wk.compute_transition_matrix(
    growth_iters=3, growth_rate_key="growth_rate_init", last_time_point="connectivities"
)
wk.plot_random_walks(
    n_sims=300,
    max_iter=200,
    start_ixs={"Timepoint": 1},
    stop_ixs={"Timepoint": 3},
    c="Timepoint",
    legend_loc="none",
    title = "Day 1 to day 3",
    linealpha=0.5,
    dpi=600, save = "Day1to3.pdf"
)
wk = WOTKernel(adata3, time_key="Timepoint")
wk.compute_initial_growth_rates(organism="mouse", key_added="growth_rate_init")
wk.compute_transition_matrix(
    growth_iters=3, growth_rate_key="growth_rate_init", last_time_point="connectivities"
)
wk.plot_random_walks(
    n_sims=300,
    max_iter=200,
    start_ixs={"Timepoint": 3},
    stop_ixs={"Timepoint": 7},
    c="Timepoint",
    legend_loc="none",
    title = "Day 3 to day 7",
    linealpha=0.5,
    dpi=600, save = "Day3to7.pdf"
)
wk = WOTKernel(adata3, time_key="Timepoint")
wk.compute_initial_growth_rates(organism="mouse", key_added="growth_rate_init")
wk.compute_transition_matrix(
    growth_iters=3, growth_rate_key="growth_rate_init", last_time_point="connectivities"
)
wk.plot_random_walks(
    n_sims=300,
    max_iter=200,
    start_ixs={"Timepoint": 7},
    stop_ixs={"Timepoint": 14},
    c="Timepoint",
    legend_loc="none",
    title = "Day 7 to day 14",
    linealpha=0.5,
    dpi=600, save = "Day7to14.pdf"
)
wk = WOTKernel(adata3, time_key="Timepoint")
wk.compute_initial_growth_rates(organism="mouse", key_added="growth_rate_init")
wk.compute_transition_matrix(
    growth_iters=3, growth_rate_key="growth_rate_init", last_time_point="connectivities"
)
wk.plot_random_walks(
    n_sims=300,
    max_iter=200,
    start_ixs={"Timepoint": 14},
    stop_ixs={"Timepoint": 28},
    c="Timepoint",
    legend_loc="none",
    title = "Day 14 to day 28",
    linealpha=0.5,
    dpi=600, save = "Day14to28.pdf"
)


# In[18]:


wk.plot_random_walks(
    n_sims=30,
    max_iter=50,
    c="Timepoint",
    legend_loc="best",
    title = "Day 0 to day 28",
    linealpha=0.75,
    dpi=600, save = "Day0to28.pdf"
)


# In[19]:


sc.tl.leiden(adata3, resolution = 1.2, key_added= "subcluster")
#sc.pl.umap(adata3, color=["leiden", "subcluster"], save="reclustering_EC_MYE.pdf")


# In[20]:


#sc.tl.rank_genes_groups(adata3, 'subcluster', method='logreg')
sc.tl.rank_genes_groups(adata3, 'subcluster', method='wilcoxon')


# In[21]:


result = pd.DataFrame(adata3.uns['rank_genes_groups']["names"]).head(100)

result.to_csv("./subcluster_marker_wilcox.csv")


# In[22]:


marker_genes_dict_EC = { "General Endothelial": ["Pecam1", "Cdh5", "Cd36", "Eng", "Esam", "Tie1", "Kdr"],
                    "Artery Endothelial":["Plk2", "Id1", "Mgll", "Magix"],
                    "Venous Endothelial": ["Cytl1", "Vwf", "Npr3", "Plvap", "Cpe", "Ece1", "Ctla2a", "Heg1"],
                    "Lymphatic Endothelial": ["Nts", "Mmrn1", "Flt4", "Reln", "Lcn2", "Lyve1", "Pard6g", "Clca3a1"],
                    "Interferon Endothelial": ["Fabp4", "Slc28a2", "Gm12002", "Tcf15", "Ly6a", "Rgcc", "Rsad2", "Hspb1"],
                     "Angiogenic Endothelial": ["Sparc", "Pcdh17", "Sox4", "Bst2", "Trp53i11"]
                    }

marker_genes_dict_Myl = {"General Myeloid": ["Ptprc", "Lyz2", "Lgals3", "Itgam", "Runx1"],
                        "Timd4+ resident": ["Folr2", "Lyve1", "Igf1", "Cd163", "Timd4"],
                         "Monocytes": ["Cd14", "Il1b", "Cx3cr1", "Adgre1", "Mmp12", "Ly6c2", "Cxcl3", "Cd86"],
                         "cDCs": ["Itgax", "Itgam", "Cd24a"],
                         "Pro-inflammatory" : ["Cd86", "H2-Aa", "H2-Ab1"],
                         "pDCs": ["Tlr7", "Bst2"],
                         "Ccr2+ circulating": ["Ccr2", "Fcgr1", "Plac8"],
                         "Cd209a cDCs": ["Cd209a", "H2-DMb2"], 
                    "Interferon stimulated": ["Irf7", "Rsad2", "Mndal", "Isg15", "Rtp4", "Ifit1", "Ifi203", "Oasl2", "Ly6a", "Ifit3", "AW112010"]}
sc.pl.dotplot(adata3, marker_genes_dict_EC, 'subcluster', dendrogram=True, save = "EC_marker_bubble.pdf")
sc.pl.dotplot(adata3, marker_genes_dict_Myl, 'subcluster', dendrogram=True, save = "Myeloid_marker_bubble.pdf")


# In[20]:


adata3.rename_categories("subcluster", ["Artery EC (0)", "Class. Monocyte (1)", "Class. Monocyte (2)",
                                       "Inflamm. Monocyte (3)", "Venous EC (4)", "Monocyte (5)", "Cxcl3+ Monocyte (6)",
                                       "Non-Clas. Monocyte (7)", "Artery EC (8)", "Cd209+ cDC (9)", "Monocyte (10)", 
                                       "Timd4+ resident Macro (11)", "Lymphatic EC (12)", "Cxcl3+ Monocyte (13)", 
                                        "Ccr2+ circulating (14)", "Interferon DC (15)", "IMEC (16)", "Artery EC (17)"])

#sc.pl.umap(adata3, color='leiden', legend_loc='on data', title='Celltypes', frameon=False, save='Celltypes_leiden.pdf')
#sc.pl.umap(adata3, color='subcluster', legend_loc='right margin', title='Subclusters', frameon=False, save='Celltypes_subcluster.pdf')


# In[21]:


save_file = './scanpy_ECs_Myeloid.h5ad'
adata3.write_h5ad(save_file)

save_file2 = './scanpy_all_cells.h5ad'
adata.write_h5ad(save_file2)


# In[24]:


sc.tl.dendrogram(adata3, groupby="subcluster")


# In[25]:


marker_genes1 = ["Pecam1", "Cdh5", "Cd36", "Fabp4", "Eng", "Esam", "Tie1", "Kdr",
                "Ptprc", "Lyz2", "Lgals3", "Itgam", "Runx1", 
               "Abhd6", "Ablim1", "Adgrf5", "Amica1", "Arg2", "Cd300lg", "Cdc14a", "Ciita", "Clec14a", "Coro1a", "Galnt15", "Hdc", "Hspa12b", "Il1f9", "Lcp1", "Mal", "Mcam", "Mmp8", "Mtss1", "Mxd1", "Myct1", "Plxnd1", "Ppp1r13b", "Ptgs2", "Rasgrp4", "Rpl3l", "Rps6ka2", "Shroom4", "Tal1", "Tek", "Tmem204", "Tyrobp", "Ushbp1"]

marker_genes = ["Pecam1", "Cdh5", "Cd36", "Fabp4", "Eng", "Esam", "Tie1", "Kdr",
                "Ptprc", "Lyz2", "Lgals3", "Itgam", "Runx1"]
ax = sc.pl.stacked_violin(adata3, marker_genes, groupby='subcluster',
                         var_group_positions=[(0, 7), (8,12)], var_group_labels=['Endothelial', "Myeloid"], dendrogram=True, save = "Violin_general_markers.pdf")
ax = sc.pl.dotplot(adata3, marker_genes1, groupby='subcluster',
                         var_group_positions=[(0, 7), (8,12)], var_group_labels=['Endothelial', "Myeloid"], dendrogram=True, save = "cluster16.pdf")


# In[26]:


ax = wk.plot_single_flow(
    cluster_key="subcluster",
    time_key="Timepoint",
    cluster = "Artery EC (0)",
    min_flow=0.2,
    xticks_step_size=1,
    show=False,
    dpi=300)

import matplotlib.pyplot as plt

# prettify the plot a bit, rotate x-axis tick labels
locs, labels = plt.xticks()
ax.set_xticks(locs)
ax.set_xticklabels(labels, rotation=90)

plt.show()


# In[27]:


with rc_context({'figure.figsize': (4.5, 3)}):
    sc.pl.violin(adata3, ["Cd44", "Nfkbia"], groupby='subcluster', scale = "width", 
                 rotation = 90, save = "Cd44.pdf", order = ["IMEC (16)", "Artery EC (0)", "Artery EC (8)", "Venous EC (4)", "Lymphatic EC (12)", "Artery EC (17)"])


# In[28]:


with rc_context({'figure.figsize': (4.5, 3)}):
    sc.pl.violin(adata3, ["Fn1", "Col1a2", "S100a4", "Postn", "Mmp14"], groupby='subcluster', scale = "width", 
                 rotation = 90, save = "EndMT.pdf", order = ["IMEC (16)", "Artery EC (0)", "Artery EC (8)", "Venous EC (4)", "Lymphatic EC (12)", "Artery EC (17)"])


# In[29]:


adata3.obs.to_csv("Obs.csv")


# In[ ]:





# In[30]:


ck = ConnectivityKernel(adata3)
ck.compute_transition_matrix()

combined_kernel = 0.9 * wk + 0.1 * ck


# In[31]:


g = GPCCA(combined_kernel)


# In[ ]:


g.compute_schur(n_components= 25)
g.plot_spectrum(real_only=True)


# In[ ]:


g.compute_macrostates(n_states=8, cluster_key="subcluster")
g.plot_macrostates(discrete=True, legend_loc="right", save = "Macrostates.pdf")


# In[ ]:


g.plot_macrostate_composition(key="Timepoint", save = "Macrostate_composition.pdf")


# In[ ]:


g.set_terminal_states_from_macrostates(["Artery EC (8)", "Venous EC (4)_1", "Timd4+ resident Macro (11)", "Non-Clas. Monocyte (7)"])
g.compute_absorption_probabilities(solver="gmres", use_petsc=True)


# In[ ]:


g.plot_absorption_probabilities(same_plot=False, perc=[30, 99], save = "Absorption probabilities.pdf")


# In[ ]:


cr.pl.circular_projection(adata3, keys=["cell_sets", "Timepoint"], legend_loc="right", title="")


# In[ ]:


drivers = g.compute_lineage_drivers(lineages= "Endothelial cells_3", return_drivers=True)
drivers.sort_values(by="Endothelial cells_3_corr", ascending=False)


# In[ ]:


g.plot_lineage_drivers("Endothelial cells_3", n_genes=12)


# In[ ]:


wk.plot_random_walks(100, stop_ixs = {"Timepoint"}, max_iter=100,
    show_progress_bar=True,
    ixs_legend_loc="best",
    seed=42)


# In[ ]:


sc.tl.leiden(adata3, resolution = 1.2, key_added= "subcluster")
sc.pl.umap(adata3, color=["leiden", "subcluster"], save="reclustering.pdf")


# In[ ]:


ax = sc.pl.dotplot(adata3, ["Acta2", "Postn", "Col1a1", "Col1a2", "Col3a1", "Pdgfra", "Col4a1", "Cd74", "Lyz2", "Itgam", "Fcgr4", "Ccr2", "Igf1", "Timd4", "Siglec1", "Cd14", "Cxcr2", "Lgals3", "Napsa", "Cd209a", "Itgb7", "Clec10a", "Flt3", "Vwf", "Eng", "Emcn", "Fabp4", "Pecam1", "Cdh5", "Tie1", "Lyve1", "Pdpn", "Prox1", "Ccl5", "Nkg7", "Ptprc", "Klre1", "Cd4", "Cd8a", "Ctla4", "Icos", "Cd69", "Il2ra", "Itgae", "Cd3e", "Cd44", "Ccr7", "Lat", "Cd79a", "Cd79b", "Pax5", "Plp1", "Kcna1", "Cd59a", "Rgs5", "Pdgfrb", "Tagln", "Des", "Krt19", "Krt8", "Upk3b", "Tnnt2", "Nkain4"], 
                          groupby='subcluster')


# In[ ]:


from cellrank.tl.kernels import VelocityKernel

vk = VelocityKernel(adata3)
vk.compute_transition_matrix()


# In[ ]:


from cellrank.tl.kernels import ConnectivityKernel

ck = ConnectivityKernel(adata).compute_transition_matrix()


# In[ ]:


combined_kernel = 0.8 * vk + 0.2 * ck


# In[ ]:


from cellrank.tl.estimators import GPCCA

g = GPCCA(combined_kernel)
g.compute_schur(n_components=20)
g.plot_spectrum()

