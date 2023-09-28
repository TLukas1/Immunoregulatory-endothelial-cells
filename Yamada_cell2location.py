#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
IN_COLAB = "google.colab" in sys.modules
if IN_COLAB:
    get_ipython().system('pip install --quiet scvi-colab')
    from scvi_colab import install
    install()
    get_ipython().system('pip install --quiet git+https://github.com/BayraktarLab/cell2location#egg=cell2location[tutorials]')


# In[5]:


import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os

import cell2location

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text for PDFs


# In[3]:


results_folder = './results'

# create paths and names to results folders for reference regression and cell2location models
ref_run_name = f'{results_folder}/reference_signatures'
run_name = f'{results_folder}/cell2location_map'


# In[165]:


# Load Yamada

# Directory where your .h5ad files are stored
input_directory = '/media/tombor/Helios_scStorage/LukasZ/Spatial/Yamada_mouse_AMI/Analysis/scanpy/objects/'

# Get list of all .h5ad files in the directory
h5ad_files = [f for f in os.listdir(input_directory) if f.endswith('.h5ad')]

# Load each .h5ad file into a dictionary with the filename (minus extension) as the key
adata_samples = {}
for h5ad_file in h5ad_files:
    sample_name = os.path.splitext(h5ad_file)[0]  # remove .h5ad extension to get sample name
    adata_samples[sample_name] = sc.read(os.path.join(input_directory, h5ad_file))

# Now, adata_samples is a dictionary where keys are sample names and values are loaded adata objects


# In[40]:


# Load adata3 reference data
adata_ref = sc.read(filename= "/media/tombor/Helios_scStorage/Lukas/Analysis/EHT/Combined_Seurats/scanpy_all_cells.h5ad")
adata3 = sc.read(filename= "/media/tombor/Helios_scStorage/Lukas/Analysis/EHT/Combined_Seurats/scanpy_ECs_Myeloid.h5ad")


# In[166]:


for sample_name, adata in adata_samples.items():
    adata.var['SYMBOL'] =  adata.var.index
    # Update the adata object in the dictionary
    adata_samples[sample_name] = adata


# In[167]:


for sample_name, adata in adata_samples.items():
    # find mitochondria-encoded (MT) genes
    adata.var['MT_gene'] = [gene.startswith('MT-') for gene in adata.var['SYMBOL']]

    # remove MT genes for spatial mapping (keeping their counts in the object)
    adata.obsm['MT'] = adata[:, adata.var['MT_gene'].values].X.toarray()
    adata = adata[:, ~adata.var['MT_gene'].values]
    # Update the adata object in the dictionary
    adata_samples[sample_name] = adata    


# In[77]:


adata_ref.X = adata_ref.layers["matrix"].toarray()

adata_ref.var['SYMBOL'] = adata_ref.var.index

from cell2location.utils.filtering import filter_genes

selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)

# filter the object
adata_ref = adata_ref[:, selected].copy()


# In[83]:


# Create a new column 'celltypes' with the values from 'leiden' as default
adata_ref.obs['celltypes'] = adata_ref.obs['leiden'].astype(str)

# Identify which cells from adata3 have a 'subcluster' in the specified list
mask = adata3.obs['subcluster'].isin(["IMEC (16)", "Venous EC (4)", "Artery EC (8)", "Artery EC (0)", "Lymphatic EC (12)", "Artery EC (17)"])

# Convert the series to strings to avoid Categorical error
subset_subclusters = adata3[mask].obs['subcluster'].astype(str)

# Update the 'celltypes' column in adata_ref for those cells
adata_ref.obs.loc[subset_subclusters.index, 'celltypes'] = subset_subclusters
unique_celltypes = adata_ref.obs['celltypes'].unique()
print(unique_celltypes)


# In[84]:


# prepare anndata for the regression model
cell2location.models.RegressionModel.setup_anndata(adata=adata_ref,
                        # 10X reaction / sample / batch
                        batch_key='Sample',
                        # cell type, covariate used for constructing signatures
                        labels_key='celltypes',
                        # multiplicative technical effects (platform, 3' vs 5', donor effect)
                        categorical_covariate_keys=['Timepoint']
                       )


# In[85]:


# create the regression model
from cell2location.models import RegressionModel
mod = RegressionModel(adata_ref)

# view anndata_setup as a sanity check
mod.view_anndata_setup()


# In[87]:


mod.train(max_epochs=200, use_gpu=True)


# In[88]:


mod.plot_history(20)


# In[89]:


# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_ref = mod.export_posterior(
    adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True}
)

# Save model
mod.save(f"{ref_run_name}", overwrite=True)

# Save anndata object with results
adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref.write(adata_file)
adata_file


# In[91]:


adata_ref = mod.export_posterior(
    adata_ref, use_quantiles=True,
    # choose quantiles
    add_to_varm=["q05","q50", "q95", "q0001"],
    sample_kwargs={'batch_size': 2500, 'use_gpu': True}
)


# In[142]:


# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']
inf_aver.iloc[0:20, 0:20]


# In[168]:


for sample_name, adata_vis in list(adata_samples.items())[:3]:
    # find shared genes and subset both anndata and reference signatures
    Genes = inf_aver.index.tolist()
    intersect = np.intersect1d(adata_vis.var_names, Genes)
    print(intersect)
    adata_vis = adata_vis[:, intersect].copy()
    inf_aver = inf_aver.loc[intersect, :].copy()
    adata_samples[sample_name] = adata_vis
    


# In[173]:


for sample_name, adata_vis in list(adata_samples.items())[:3]:
    # prepare anndata for cell2location model
    cell2location.models.Cell2location.setup_anndata(adata=adata_vis)
    # create and train the model
    mod = cell2location.models.Cell2location(
    adata_vis, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=10,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=20
    )
    mod.view_anndata_setup()
    mod.train(max_epochs=3000,
          # train using full data (batch_size=None)
          batch_size=None,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1,
          use_gpu=True,
         )

    # plot ELBO loss history during training, removing first 100 epochs from the plot
    mod.plot_history(1000)
    plt.legend(labels=['full data training']);
    # In this section, we export the estimated cell abundance (summary of the posterior distribution).
    adata_vis = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True}
    )

    # Save model
    mod.save(f"{run_name}{sample_name}", overwrite=True)

    # mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

    # Save anndata object with results
    adata_file = f"{sample_name}_sp.h5ad"
    adata_vis.write(adata_file)
    adata_file


# In[175]:


# Load Yamada c2l files

# Directory where your .h5ad files are stored
input_directory = '/media/tombor/Helios_scStorage/LukasZ/Spatial/Yamada_mouse_AMI/Analysis/scanpy/objects/c2l'

# Get list of all .h5ad files in the directory
h5ad_files = [f for f in os.listdir(input_directory) if f.endswith('.h5ad')]

# Load each .h5ad file into a dictionary with the filename (minus extension) as the key
adata_c2l = {}
for h5ad_file in h5ad_files:
    sample_name = os.path.splitext(h5ad_file)[0]  # remove .h5ad extension to get sample name
    adata_c2l[sample_name] = sc.read(os.path.join(input_directory, h5ad_file))

# Now, adata_samples is a dictionary where keys are sample names and values are loaded adata objects


# In[178]:


print(adata_c2l)


# In[189]:


adata = adata_c2l['MI_d1_1_sp']  
# add 5% quantile, representing confident cell abundance, 'at least this amount is present',
# to adata.obs with nice names for plotting
adata.obs[adata.uns['mod']['factor_names']] = adata.obsm['q05_cell_abundance_w_sf']
sc.pl.spatial(adata, cmap='magma',
                  # show first 8 cell types
                  color=["IMEC (16)","T-cells", "Artery EC (0)", "Fibroblasts"],
                  ncols=2, size=1.3,
                  img_key='hires', scale_factor=1,
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2', save = "IMEC_MI_d1_1"
                 )


# In[190]:


adata = adata_c2l['MI_d1_2_sp']  
# add 5% quantile, representing confident cell abundance, 'at least this amount is present',
# to adata.obs with nice names for plotting
adata.obs[adata.uns['mod']['factor_names']] = adata.obsm['q05_cell_abundance_w_sf']
sc.pl.spatial(adata, cmap='magma',
                  # show first 8 cell types
                  color=["IMEC (16)","T-cells", "Artery EC (0)", "Fibroblasts"],
                  ncols=2, size=1.3,
                  img_key='hires', scale_factor=1,
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2', save = "IMEC_MI_d1_2"
                 )


# In[191]:


adata = adata_c2l['MI_d1_3_sp']  
# add 5% quantile, representing confident cell abundance, 'at least this amount is present',
# to adata.obs with nice names for plotting
adata.obs[adata.uns['mod']['factor_names']] = adata.obsm['q05_cell_abundance_w_sf']
sc.pl.spatial(adata, cmap='magma',
                  # show first 8 cell types
                  color=["IMEC (16)","T-cells", "Artery EC (0)", "Fibroblasts"],
                  ncols=2, size=1.3,
                  img_key='hires', scale_factor=1,
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2', save = "IMEC_MI_d1_3"
                 )


# In[184]:


print(adata)

