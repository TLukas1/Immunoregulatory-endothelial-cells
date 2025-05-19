# Import Libraries
import os
import warnings 

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt

warnings.filterwarnings('ignore')

import davidrUtility

input_directory = '/mnt/davidr/scStorage/LukasZ/Spatial/Yamada_mouse_AMI/Analysis/scanpy/objects/c2l'
h5ad_files = [f for f in os.listdir(input_directory) if f.endswith('.h5ad')]

# Load each .h5ad file into a dictionary with the filename (minus extension) as the key
adata_c2l = {}
for h5ad_file in h5ad_files:
    sample_name = os.path.splitext(h5ad_file)[0]  # remove .h5ad extension to get sample name
    adata_c2l[sample_name] = sc.read(os.path.join(input_directory, h5ad_file))
    adata_c2l[sample_name].obs[adata_c2l[sample_name].uns['mod']['factor_names']] = adata_c2l[sample_name].obsm['q05_cell_abundance_w_sf']
    
# Adapt metadata information
adata1 = adata_c2l['MI_d1_1_sp']  
adata1.obs['sample'] = 'MI_d1_1'
adata1.uns['spatial']['MI_d1_1'] = adata1.uns['spatial']['default'] 
del adata1.uns['spatial']['default']

adata2 = adata_c2l['MI_d1_2_sp']
adata2.obs['sample'] = 'MI_d1_2'
adata2.uns['spatial']['MI_d1_2'] = adata2.uns['spatial']['default'] 
del adata2.uns['spatial']['default']

adata3 = adata_c2l['MI_d1_3_sp']
adata3.obs['sample'] = 'MI_d1_3'
adata3.uns['spatial']['MI_d1_3'] = adata3.uns['spatial']['default'] 
del adata3.uns['spatial']['default']

adata = ad.concat([adata1, adata2, adata3], uns_merge='unique')
adata.obs['sample'] = pd.Categorical(adata.obs['sample'])


# Calculate the proportion of IMECs on T cell spots
imec_tcell = []
for batch in adata.obs['sample'].unique():
    adata_subset = davidrUtility.select_slide(adata, batch)

    adata_subset.obs['IMEC_spots'] =  pd.Categorical(adata_subset.obs['IMEC (16)'] > np.percentile(adata_subset.obs['IMEC (16)'], 90)).map({False:'Nothing', True:'IMEC'})
    adata_subset.obs['Tcells_spots'] =  pd.Categorical(adata_subset.obs['T-cells'] > np.percentile(adata_subset.obs['T-cells'], 90)).map({False:'Nothing', True:'T_cells'})
    
    # sc.pl.spatial(adata_subset, color=[ 'IMEC_spots', 'Tcells_spots', 'ECs_spots',], size=1.3, scale_factor=1, palette={'Nothing':'gray', 'IMEC':'firebrick', 'T_cells':'sandybrown','ECs':'royalblue'})

    imec_tcell.append(len(np.where((adata_subset.obs['IMEC_spots'] == 'IMEC') & (adata_subset.obs['Tcells_spots'] == 'T_cells'))[0]) / np.sum(adata_subset.obs['Tcells_spots'] == 'T_cells')* 100) 

imec_tcell = pd.DataFrame(imec_tcell).T
df_main = pd.concat([imec_tcell])
df_main.index = ['IMEC\nT_cells']
df_main.columns = ['d1_1', 'd1_2', 'd1_3']
df_main_melt = df_main.reset_index().melt(id_vars='index')

sns.heatmap(df_main.T, cmap='Reds', annot=True, robust=True)
plt.savefig('/mnt/davidr/scStorage/DavidR/Requests/ForStefanie/IMEC_revision/Heatmap_IMEC_ECs_Tcell_Coloc_PropscRNA.svg', bbox_inches='tight')


