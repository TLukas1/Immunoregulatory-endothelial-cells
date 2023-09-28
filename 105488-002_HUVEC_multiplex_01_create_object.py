#!/usr/bin/env python
# coding: utf-8

# # Create Adata of HUVEC treatment single cell sequencing experiments

# # 1- Import

# In[ ]:


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


# # 2-Load datasets
# 
# Use the output .loom files from velocyto which are located here: [scStorage/Lukas/Multiplex.HUVEC/105488-002/combined/velocity/out](files/../velocity/out/)
# 
# Script for running velocyto: [script](/media/tombor/Helios_scStorage/Lukas/Multiplex.HUVEC/105488-002/combined/velocity/run_velocity_seq.sh)
# 
# Commands: [Commands.txt](/media/tombor/Helios_scStorage/Lukas/Multiplex.HUVEC/105488-002/combined/velocity/Commands.txt) [rawcommands.txt](/media/tombor/Helios_scStorage/Lukas/Multiplex.HUVEC/105488-002/combined/velocity/outcommands.txt)
# 
# Velocyto Tutorial: [tutorial](http://velocyto.org/)
# 
# For loop for all .loom files:

# In[ ]:


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
        size=size,
        frameon=frameon,
        legend_loc=legend_loc,
        **kwargs,
    )
    
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization


# In[ ]:


looms=list()
names=list()
mypath="/media/tombor/Helios_scStorage/Lukas/Multiplex.HUVEC/105488-002/combined/velocity/out/"

for f in glob.glob(glob.escape(mypath) + "/*.loom"):
    adata=scv.read(f,cache=True)
    adata.var_names_make_unique()
    looms.append(adata)
    names.append(f)
print(looms)


# In[ ]:


#getsamplenames 
samplenames=list()
for n in names:
    samplename = re.search("/media/tombor/Helios_scStorage/Lukas/Multiplex.HUVEC/105488-002/combined/velocity/out/(.+?)_", n).group(1)
    print(samplename)
    samplenames.append(samplename)


# In[ ]:


anno=pd.read_csv("/media/tombor/Helios_scStorage/Lukas/Multiplex.HUVEC/105488-002/combined/velocity/Anno.csv", sep = ",", index_col=0)
print(anno)
for i in range(len(looms)):
    looms[i].obs["ID"] = samplenames[i]
    looms[i].obs["Timepoint"] = anno["Timepoint"][samplenames[i]]
    looms[i].obs["Treatment"] = anno["Treatment"][samplenames[i]]
    looms[i].obs["Batch"] = anno["Batch"][samplenames[i]]
    looms[i].obs["Sample"] = anno["Sample"][samplenames[i]]


# In[ ]:


processed_adatas=list()
for adata in looms:
    print(adata.shape)
    print("Filter cells and genes")
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    print("MT-filter")
    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')
    print("save to list...")
    processed_adatas.append(adata)


# In[ ]:


adata = anndata.concat(processed_adatas, join = "outer")
print(adata)


# In[ ]:


adata = adata[adata.obs.total_counts < 25000, :]
adata = adata[adata.obs.pct_counts_mt < 35, :]
print(adata)


# In[ ]:


sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.001, max_mean=4, min_disp=0.2)
sc.pl.highly_variable_genes(adata)


# In[ ]:


sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata)


# In[ ]:


sc.external.pp.bbknn(adata, batch_key="Treatment", n_pcs = 10)
sc.tl.leiden(adata, resolution = 0.2)
sc.tl.umap(adata, min_dist = 0.2)


# In[ ]:


sc.pl.umap(adata, color=['leiden'], save= "105488-002_HUVEC_treatment_d7_UMAP.pdf")
sc.pl.umap(adata, color=['Treatment'], save= "105488-002_HUVEC_treatment_d7_UMAP_Treatment.pdf")
sc.pl.umap(adata, color=['Sample'], save= "105488-002_HUVEC_treatment_d7_UMAP_Sample.pdf")


# In[ ]:


# Save this
save_file = '/media/tombor/Helios_scStorage/Lukas/Multiplex.HUVEC/scanpy/objects/HUVEC_d7_merged.h5ad'
adata.write_h5ad(save_file)


# In[ ]:


# Read file
adata = sc.read_h5ad('/media/tombor/Helios_scStorage/Lukas/Multiplex.HUVEC/scanpy/objects/HUVEC_d7_merged.h5ad')


# In[ ]:


# Basic preprocessing velocity
scv.pl.proportions(adata)
#scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=None, n_neighbors=None)


# In[ ]:


# run scVelo
scv.tl.recover_dynamics(adata, n_jobs=30)


# In[ ]:


# built velocity graph
scv.tl.velocity(adata, mode="dynamical")
scv.tl.velocity_graph(adata)


# In[ ]:


scv.pl.velocity_embedding_stream(
    adata, basis="umap", legend_fontsize=12, title="", smooth=0.8, min_mass=3, color="Treatment", save = "velocity.pdf")


# In[ ]:


scv.pl.velocity_embedding(adata, arrow_length=3, arrow_size=0.1, dpi=120, save= ".pdf")


# In[ ]:


scv.pl.velocity(adata, ['IL6', 'TGFBI'])


# In[ ]:


scv.tl.rank_velocity_genes(adata, groupby='Treatment', min_corr=.3)
df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
df.head(n=10)


# In[ ]:


kwargs = dict(frameon=False, size=10, linewidth=1.5)

scv.pl.scatter(adata, df['IFNG'][:5], ylabel='INFg', **kwargs)
scv.pl.scatter(adata, df['TGF-B IL-1B'][:5], ylabel='TGF-B IL-1B', **kwargs)
scv.pl.scatter(adata, df['TGF-B IL-1B IFN-G'][:5], ylabel='TGF-B IL-1B INFg', **kwargs)


# In[ ]:


scv.tl.velocity_pseudotime(adata)
scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot')


# In[ ]:


sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(3, 3), facecolor='white')  # low dpi (dots per inch) yields small inline figures


# In[ ]:


sc.tl.diffmap(adata)
sc.pp.neighbors(adata, n_neighbors=10, use_rep='X_diffmap')
sc.tl.draw_graph(adata)


# In[ ]:


sc.pl.draw_graph(adata, color='Treatment', legend_loc='on data')


# In[ ]:


sc.tl.leiden(adata, resolution = 0.2, key_added="PAGA-0.2")
sc.tl.paga(adata, groups="PAGA-0.2")


# In[ ]:


# Compare DEGs against Homeostasis against Timepoint
adata_copy = adata
sc.tl.rank_genes_groups(adata_copy, 'PAGA-0.2', method='wilcoxon', pts=True, use_raw =False)
sc.pl.rank_genes_groups(adata_copy, n_genes=25, sharey=False)
result = adata_copy.uns['rank_genes_groups']
groups = result['names'].dtype.names
degs_by_cluster = pd.DataFrame({group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'logfoldchanges', 'pvals_adj']})
degs_by_cluster.to_csv("DEG_PAGA-0.2.csv")


# In[ ]:


sc.pl.paga(adata, color=['PAGA-0.2', "Treatment"], save = "PAGA_general.pdf")

sc.pl.paga(adata, color=["RUNX1", "TNC", "IFI44", "CXCL1"], save = "PAGA_0.pdf")

sc.pl.paga(adata, color=["TOP2A", "MKI67", "TYMS", "HIST1H1B"], save = "PAGA_1.pdf")

sc.pl.paga(adata, color=["IFI6", "ICAM1", "CDKN1A", "TAGLN"], save = "PAGA_2.pdf")

sc.pl.paga(adata, color=["PTX3", "KRT18", "ODC1", "AXL"], save = "PAGA_3.pdf")

sc.pl.paga(adata, color=["VWF", "AQP1", "SOX18", "MARCKSL1"], save = "PAGA_4.pdf")

sc.pl.paga(adata, color=["HLA-DPB1", "HLA-DPA1", "HLA-DQA1", "HLA-DQB1"], save = "PAGA_5.pdf")

sc.pl.paga(adata, color=["CCNA2", "CCNB1", "AURKA", "BIRC5"], save = "PAGA_6.pdf")

sc.pl.paga(adata, color=["TGFBI", "MYL9", "TGFBR1", "TGFB2"], save = "PAGA_7.pdf")

sc.pl.paga(adata, color=["HLA-DOA", "BMP4", "BMX", "TCN2"], save = "PAGA_8.pdf")

sc.pl.paga(adata, color=["CXCL8", "CXCL1", "CCL5", "IL6"], save = "PAGA_9.pdf")

sc.pl.paga(adata, color=["MMP1", "ITGA6", "APLN", "LMNA"], save = "PAGA_10.pdf")

sc.pl.paga(adata, color=["MYH10", "IL33", "EMCN", "JUP"], save = "PAGA_11.pdf")

sc.pl.paga(adata, color=["IL7R", "CCL8", "NFKBIZ", "IFI35"], save = "PAGA_12.pdf")

sc.pl.paga(adata, color=["CXCL2", "SELE", "LIPG", "SULF1"], save = "PAGA_13.pdf")

sc.pl.paga(adata, color=["TMSB10", "LGALS1", "PFN1", "RPLP1"], save = "PAGA_14.pdf")

sc.pl.paga(adata, color=["LAMB1", "FAT4", "NID1", "WSB1"], save = "PAGA_15.pdf")

sc.pl.paga(adata, color=["ADGRF5", "RUNX1T1", "ABLIM1", "MAGI1"], save = "PAGA_17.pdf")


# In[ ]:


sc.pl.paga(adata, color=["BATF2", "BATF3", "KDM6B", "KDM7A"])
sc.pl.paga(adata, color=["CIITA", "CNN1", "CD40", "TNFSF4"])
sc.pl.paga(adata, color=["TNFSF9", "CD80", "CD274", "PDCD1LG2"])


# In[ ]:


sc.tl.draw_graph(adata, init_pos='paga')


# In[ ]:


sc.pl.draw_graph(adata, color=["PAGA-0.2", "IL6", "TGFB2", "HLA-DPB1", "TOP2A", "AQP1"], legend_loc='on data', save = "PAGA_maps01.pdf")
sc.pl.draw_graph(adata, color=["CXCL8", "CCL8", "RUNX1", "RUNX1T1"], legend_loc='on data', save = "PAGA_maps02.pdf")


# In[ ]:


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
new_colors[[5,8,10]] = zeileis_colors[[17,17,16]]  # IFNG / orange
new_colors[[17,4,11,15,16]] = zeileis_colors[[0,6,6,1,25]]  # Baseline / blue
new_colors[[1,3,14]] = zeileis_colors[[23,22,22]]  # Prolif / pink

adata.uns['PAGA-0.2_colors'] = new_colors


# In[ ]:


sc.pl.paga(adata, color=['PAGA-0.2', "Treatment"], save = "PAGA_general.pdf")


# In[ ]:


sc.pl.paga_compare(
    adata, threshold=0.1, title='', right_margin=0.2, size=10, edge_width_scale=0.5,
    legend_fontsize=12, fontsize=12, frameon=False, edges=True, save=True)


# In[ ]:


adata.uns['iroot'] = np.flatnonzero(adata.obs['PAGA-0.2']  == '17')[0]
sc.tl.dpt(adata)


# In[ ]:


adata_raw = sc.read_h5ad('/media/tombor/Helios_scStorage/Lukas/Multiplex.HUVEC/scanpy/objects/HUVEC_d7_merged.h5ad')
sc.pp.log1p(adata_raw)
sc.pp.scale(adata_raw)
adata.raw = adata_raw


# In[ ]:


sc.pl.draw_graph(adata, color=['PAGA-0.2'], save = "PAGA_UMAP.pdf")
sc.pl.draw_graph(adata, color=['dpt_pseudotime'], save = "Pseudo_UMAP.pdf")
sc.pl.draw_graph(adata, color=['Treatment'], save = "Treatment_UMAP.pdf")


# In[ ]:


cluster_small_multiples(adata, clust_key = "Treatment", save = "Treatment_highlighted.pdf")


# In[ ]:


# Save this
save_file = '/media/tombor/Helios_scStorage/Lukas/Multiplex.HUVEC/scanpy/objects/HUVEC_d7_merged_PAGA02.h5ad'
adata.write_h5ad(save_file)


# In[ ]:


sc.pl.draw_graph(adata, color=["BATF2", "BATF3", "KDM7A", "KDM6B"], legend_loc='on data', save = "PAGA_maps03.pdf")


# In[ ]:


adata.obs['distance'] = adata.obs['dpt_pseudotime']
paths = [('INFG_1', [17, 4, 5, 8]),
         ('IMEC_1', [17, 4, 5, 8, 13, 9, 12]),
         ('IMEC_2', [17, 4, 8, 13, 0]),
         ('IMEC_3', [17, 4, 5, 8, 10, 13, 9, 12]),
         ('Prolif_EndMT_1', [6, 3, 2, 7]),
         ('Prolif_EndMT_2', [1, 16, 3, 2, 7])]
adata.obs['clusters'] = adata.obs['PAGA-0.2']  # just a cosmetic change
adata.uns['clusters_colors'] = adata.uns['PAGA-0.2_colors']


# In[ ]:


get_ipython().system('mkdir write')


# In[ ]:


gene_names = ['VWF', 'KRT18', 'AQP1', 'HLA-DPB1','HLA-DPA1', 'CIITA', 'BATF2', 'MMP1', 'ICAM1', 'NFKBIZ', 'RUNX1', 'BATF3', 'KDM7A', 'KDM6B', 'TAGLN', 'MYL9', 'TGFB2', 'CXCL3', 'TNC']


# In[ ]:


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


# In[3]:


import pandas as pd
import scanpy as sc

adata = sc.read_h5ad('/media/tombor/Helios_scStorage/Lukas/Multiplex.HUVEC/scanpy/objects/HUVEC_d7_merged_PAGA02.h5ad')


# In[18]:


ax = sc.pl.heatmap(adata, ["VWF", "AQP1", "KRT18", "HLA-DPA1", "HLA-DPB1", "CIITA", "BATF2", "MMP1", "MYL9", "TAGLN", "HEY1", "HES2", "TGFBI", "BATF3", "KDM7A", "KDM6B", "CXCL8", "IL6", "RUNX1"], groupby='Treatment', swap_axes=False, dendrogram=False, standard_scale = "var", cmap = "rocket", save = True)

