#!/usr/bin/env python
# coding: utf-8

# # Recreate the PAGA graph for publication

# In[1]:


import scvelo as scv
import scanpy
import os
import pandas as pd
import glob
import numpy as np
import scanpy as sc
import re
import anndata


scv.set_figure_params()


# In[26]:


def cluster_small_multiples(
    adata, clust_key, size=60, frameon=False, legend_loc=None, **kwargs
):
    tmp = adata.copy()

    for i, clust in enumerate(adata.obs[clust_key].cat.categories):
        tmp.obs[clust] = adata.obs[clust_key].isin([clust]).astype("category")
        tmp.uns[clust + "_colors"] = ["#d3d3d3", adata.uns[clust_key + "_colors"][i]]

    sc.pl.draw_graph(
        tmp,
        groups=tmp.obs[clust].cat.categories[1:].values,
        color=adata.obs[clust_key].cat.categories.tolist(),
        edges=True,
        size=size,
        frameon=frameon,
        legend_loc=legend_loc,
        **kwargs,
    )
    


# In[3]:


adata = sc.read_h5ad('/media/tombor/Helios_scStorage/Lukas/Multiplex.HUVEC/scanpy/objects/HUVEC_d7_merged_PAGA02.h5ad')


# In[4]:


import matplotlib.pyplot as pl
pl.figure(figsize=(8, 2))
for i in range(28):
    pl.scatter(i, 1, c=sc.pl.palettes.zeileis_28[i], s=200)
pl.show()
zeileis_colors = np.array(sc.pl.palettes.zeileis_28)
new_colors = np.array(adata.uns['PAGA-0.2_colors'])
new_colors[[7,6]] = zeileis_colors[[12, 13]]  # EndMT / green
new_colors[[9,12]] = zeileis_colors[[11, 5]]  # IMEC / red
new_colors[[0,2,13]] = zeileis_colors[[9, 4, 10]]  # Transitions / red
new_colors[[4,5,8,10]] = zeileis_colors[[15,14,17,16]]  # IFNG / orange
new_colors[[17,11,15,16]] = zeileis_colors[[0,6,1,25]]  # Baseline / blue
new_colors[[1,3,14]] = zeileis_colors[[23,22,21]]  # Prolif / pink

adata.uns['PAGA-0.2_colors'] = new_colors


# In[5]:


sc.pl.paga(adata, color=['PAGA-0.2', "Treatment"], save = "PAGA_general.pdf")


# In[34]:


sc.pl.paga_compare(
    adata, threshold=0.1, title='', right_margin=0.2, size=10, edge_width_scale=0.8,
    legend_fontsize=12, fontsize=12, frameon=False, edges=True, save=True)


# In[39]:


sc.pl.paga_compare(
    adata, threshold=0.1, color=['RUNX1'], title='', right_margin=0.6, size=10, edge_width_scale=0.8, fontsize=12, frameon=False, edges=True, save='PAGA_Runx1.pdf')


# In[42]:


sc.pl.paga_compare(
    adata, threshold=0.1, color=['HLA-DPA1'], title='', right_margin=0.6, size=10, edge_width_scale=0.8, fontsize=12, frameon=False, edges=True, save='PAGA_HLA_DPA1.pdf')


# In[45]:


sc.pl.paga_compare(
    adata, threshold=0.1, color=['CIITA'], title='', right_margin=0.6, size=10, edge_width_scale=0.8, fontsize=12, frameon=False, edges=True, save='PAGA_CIITA.pdf')


# In[117]:


gene_names = ['VWF', 'AQP1', "SOX18", "RUNX1T1", 'HLA-DPB1','HLA-DPA1', 'CIITA', "VCAM1", 'ICAM1', "SELE", 'RUNX1', "TGFB2", "IL6", "CCL7", "CXCL8", "CXCL1", "CCL5"]


# In[115]:


adata.obs['distance'] = adata.obs['dpt_pseudotime']
paths = [('INFG', [17, 4, 5, 8]),
         ('IMEC', [17, 4, 5, 8, 13, 9, 12]),
         ('IMEC_2', [17, 4, 8, 13, 0]),
         ('IMEC_3', [17, 4, 5, 8, 10, 13, 9, 12]),
         ('EndMT', [3, 2, 7]),
         ('Prolif_EndMT_2', [1, 16, 3, 2, 7])]
adata.obs['clusters'] = adata.obs['PAGA-0.2']  # just a cosmetic change
adata.uns['clusters_colors'] = adata.uns['PAGA-0.2_colors']


# In[118]:


_, axs = pl.subplots(ncols=6, figsize=(14, 5), gridspec_kw={'wspace': 0.05, 'left': 0.12})
pl.subplots_adjust(left=0.05, right=0.98, top=0.82, bottom=0.2)
for ipath, (descr, path) in enumerate(paths):
    _, data = sc.pl.paga_path(
        adata, path, gene_names,
        show_node_names=False,
        ax=axs[ipath],
        ytick_fontsize=12,
        left_margin=0.15,
        n_avg=50,
        annotations=['distance'],
        show_yticks=True if ipath==0 else False,
        show_colorbar=False,
        color_map='Greys',
        groups_key='clusters',
        color_maps_annotations={'distance': 'viridis'},
        title='{} path'.format(descr),
        return_data=True,
        show=False)
    data.to_csv('./write/paga_path_{}.csv'.format(descr))
pl.savefig('./figures/paga_path_paul15.pdf')
pl.show()


# In[36]:


sc.pl.paga_compare(
    adata, threshold=0.1, color=['Treatment'], title='', right_margin=0.2, size=10, edge_width_scale=0.8,
 fontsize=12, frameon=False, edges=True, save='PAGA_Treatment.pdf')

