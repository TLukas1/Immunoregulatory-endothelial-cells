#!/usr/bin/env python
# coding: utf-8

# In[32]:


import scanpy as sc
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns


# In[33]:


sc.logging.print_versions()
sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3


# In[34]:


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


# In[35]:


print(adata_samples)


# In[36]:


for sample_name, adata in adata_samples.items():
    # Make variable names unique
    adata.var_names_make_unique()
    
    # Identify mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    
    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    
    #Normalize and Scale
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, key_added="clusters")
    
    # Update the adata object in the dictionary
    adata_samples[sample_name] = adata


# In[40]:


for sample_name, adata in adata_samples.items():
    plt.rcParams["figure.figsize"] = (8, 8)
    #scf= adata.uns['scalefactors']['tissue_lowres_scalef']
    #spotsize= adata.uns['scalefactors']['spot_diameter_fullres']
    sc.pl.spatial(adata, img_key="lowres", color=["Nppb", "Ankrd1","Csrp3", "Nppa", "Myh7"])


# In[51]:


for sample_name, adata in adata_samples.items():
    plt.rcParams["figure.figsize"] = (12, 8)
    #scf= adata.uns['scalefactors']['tissue_lowres_scalef']
    #spotsize= adata.uns['scalefactors']['spot_diameter_fullres']
    sc.pl.spatial(adata, img_key="hires", color=["Nppb", "Ankrd1","Csrp3", "Nppa", "Myh7"], cmap='magma', scale_factor=1, alpha_img=0.2, save=sample_name+"_AMI_markers.pdf")


# In[39]:


output_directory = '/media/tombor/Helios_scStorage/LukasZ/Spatial/Yamada_mouse_AMI/Analysis/scanpy/objects/'

# Create the output directory if it doesn't exist
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

for sample_name, adata in adata_samples.items():
    # define the save path
    save_path = os.path.join(output_directory, f"{sample_name}.h5ad")
    
    # save the AnnData object
    adata.write(save_path)

print(f"All adata samples saved to {output_directory}")


# In[31]:




