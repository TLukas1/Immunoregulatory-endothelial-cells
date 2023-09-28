#!/usr/bin/env python
# coding: utf-8

# In[3]:


# Load scanpy object and libraries

import os
import pandas as pd
import glob
import numpy as np
import scanpy as sc
from matplotlib.pyplot import rc_context


# In[ ]:


adata = sc.read(filename= "/media/tombor/Helios_scStorage/Lukas/Analysis/EHT/Combined_Seurats/scanpy_all_cells.h5ad")


# In[28]:


with rc_context({'figure.figsize': (4, 4)}):
    sc.pl.umap(adata, color=['Lyz2', 'Lgals3', 'leiden', "GFP"], save = "_myeloid_clustering.pdf")

    sc.pl.umap(adata, color=["Cdh5", "Pecam1", 'leiden', "GFP"], save = "_EC_clustering.pdf")
    
    sc.pl.umap(adata, color=['leiden'], save = "leiden_clustering.pdf")


# In[10]:


with rc_context({'figure.figsize': (4, 4)}):
    sc.pl.umap(adata, color=['Ccl5', "Nkg7", "Klra8", "Cd3e", "Cd8a", "Cd4","Ifng"], save = "_T-cells_clustering.pdf")


# In[12]:


adata3 = adata[adata.obs["leiden"].isin(["T-cells"])]

sc.tl.pca(adata3)
sc.pp.neighbors(adata3, n_pcs=15, n_neighbors=15)

sc.tl.umap(adata3)
sc.tl.leiden(adata3, resolution = 0.5)


# In[58]:


dicMarker = {
    'Ifng production': ['Ifng', 'Irf8','Klrk1', 'Klre1', 'Sash3'],
    
    'NK-cells': ['Nkg7', 'Klrd1', 'Klrc2', 'Ncr1', 'Eomes', 'Prf1'],
    'T-cells' :['Cd8a', 'Cd4', 'Cd3d', 'Cd69', 'Il7r', 'Ccr7', 'Lat', 'Icos']
}


with rc_context({'figure.figsize': (4.5, 3)}):
    ax = sc.pl.violin(adata3, ['Ifng', 'Irf8','Klrk1', 'Klre1', 'Sash3'],groupby='Timepoint')


# In[97]:


Nk = adata3[adata3.obs["leiden"].isin(["1"])]
sc.tl.rank_genes_groups(Nk, groupby='Timepoint', method='wilcoxon')

sc.pl.rank_genes_groups_dotplot(Nk, groupby = ["Timepoint"], var_names = ['Ifng', 'Irf8','Klrk1', 'Klre1', 'Sash3'] , dendrogram = False, cmap='Oranges', swap_axes = True, save = 'Ifng_bubble_time_cluster1.pdf')

Nk = adata3[adata3.obs["leiden"].isin(["2"])]
sc.tl.rank_genes_groups(Nk, groupby='Timepoint', method='wilcoxon')

sc.pl.rank_genes_groups_dotplot(Nk, groupby = ["Timepoint"], var_names = ['Ifng', 'Irf8','Klrk1', 'Klre1', 'Sash3'] , dendrogram = False, cmap='Greens', swap_axes = True, save = 'Ifng_bubble_time_cluster2.pdf')

Nk = adata3[adata3.obs["leiden"].isin(["5"])]
sc.tl.rank_genes_groups(Nk, groupby='Timepoint', method='wilcoxon')

sc.pl.rank_genes_groups_dotplot(Nk, groupby = ["Timepoint"], var_names = ['Ifng', 'Irf8','Klrk1', 'Klre1', 'Sash3'] , dendrogram = False, cmap='Reds', swap_axes = True, save = 'Ifng_bubble_time_cluster5.pdf')


# In[105]:


Nk = adata3[adata3.obs["leiden"].isin(["2"])]
sc.pl.tracksplot(Nk, ['Ifng', 'Irf8','Klrk1', 'Klre1', 'Sash3'], groupby='Timepoint', dendrogram=False, save = "Tracksplot_Cluster2.pdf")
Nk = adata3[adata3.obs["leiden"].isin(["1"])]
sc.pl.tracksplot(Nk, ['Ifng', 'Irf8','Klrk1', 'Klre1', 'Sash3'], groupby='Timepoint', dendrogram=False, save = "Tracksplot_Cluster1.pdf")
Nk = adata3[adata3.obs["leiden"].isin(["5"])]
sc.pl.tracksplot(Nk, ['Ifng', 'Irf8','Klrk1', 'Klre1', 'Sash3'], groupby='Timepoint', dendrogram=False, save = "Tracksplot_Cluster5.pdf")


# In[ ]:


sc.pl.dotplot(adata3,
              dicMarker,
              'leiden',
              dendrogram=True, var_group_rotation = 0, save = '_Tcells_dot.pdf')


# In[35]:


with rc_context({'figure.figsize': (4, 4)}):
    sc.pl.umap(adata3, color=['Ccl5', "Nkg7", "Klra8", "Cd3e", "Cd8a", "Cd4","Ifng", "leiden"], save = "T-cell.pdf")
    
    ax = sc.pl.dotplot(adata, ['Ccl5', "Nkg7", "Klra8", "Cd3e", "Cd8a", "Cd4","Ifng"], 
                       groupby='leiden', save = "Tcell_markers.pdf")
     
    ax = sc.pl.dotplot(adata3, ['Ccl5', "Nkg7", "Klra8", "Cd3e", "Cd8a", "Cd4","Ifng"], dendrogram = True, 
                       groupby='leiden', save = "ReclusteringTcell_markers.pdf")


# In[24]:


ax = sc.pl.correlation_matrix(adata3, 'leiden', figsize=(5,3.5))

with rc_context({'figure.figsize': (4.5, 3)}):
    sc.pl.violin(adata3,['Ifng'], groupby='leiden' )


# In[25]:


sc.tl.rank_genes_groups(adata3, 'leiden', method='wilcoxon', key_added = "wilcoxon")
sc.pl.rank_genes_groups(adata3, n_genes=25, sharey=False, key="wilcoxon")


# In[22]:


adata3.write("Tcells_22_07_21")
adata3.write_csvs("Tcells")

