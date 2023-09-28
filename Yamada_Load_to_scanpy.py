#!/usr/bin/env python
# coding: utf-8

# In[24]:


from pathlib import Path, PurePath
from typing import Union, Dict, Optional, Tuple, BinaryIO, Literal

import h5py
import json
import numpy as np
import pandas as pd
from matplotlib.image import imread
import anndata
from anndata import (
    AnnData,
    read_csv,
    read_text,
    read_excel,
    read_mtx,
    read_loom,
    read_hdf,
)
from anndata import read as read_h5ad


def process_visium_data(
    adata: AnnData,
    path: Union[str, Path],
    library_id: Optional[str] = None,
    load_images: Optional[bool] = True,
    source_image_path: Optional[Union[str, Path]] = None
) -> AnnData:
    """\
    Process 10x-Genomics-formatted visum dataset already loaded into an AnnData object.
    
    It looks for the `spatial` folder and loads images, coordinates, and scale factors.
    """
    path = Path(path)

    # Check if the spatial information is already present
    if "spatial" not in adata.uns:
        adata.uns["spatial"] = dict()
    
    if library_id is None:
        library_id = "default"

    adata.uns["spatial"][library_id] = dict()

    if load_images:
        tissue_positions_file = (
            path / "spatial/tissue_positions.csv"
            if (path / "spatial/tissue_positions.csv").exists()
            else path / "spatial/tissue_positions_list.csv"
        )
        
        files = dict(
            tissue_positions_file=tissue_positions_file,
            scalefactors_json_file=path / 'spatial/scalefactors_json.json',
            hires_image=path / 'spatial/tissue_hires_image.png',
            lowres_image=path / 'spatial/tissue_lowres_image.png',
        )

        for f in files.values():
            if not f.exists() and not any(x in str(f) for x in ["hires_image", "lowres_image"]):
                raise OSError(f"Could not find '{f}'")
        
        adata.uns["spatial"][library_id]['images'] = dict()
        for res in ['hires', 'lowres']:
            adata.uns["spatial"][library_id]['images'][res] = imread(files[f'{res}_image'])
        
        adata.uns["spatial"][library_id]['scalefactors'] = json.loads(files['scalefactors_json_file'].read_bytes())

        positions = pd.read_csv(
            files['tissue_positions_file'],
            header=1 if tissue_positions_file.name == "tissue_positions.csv" else None,
            index_col=0,
        )
        positions.columns = [
            'in_tissue',
            'array_row',
            'array_col',
            'pxl_col_in_fullres',
            'pxl_row_in_fullres',
        ]

        adata.obs = adata.obs.join(positions, how="left")
        adata.obsm['spatial'] = adata.obs[['pxl_row_in_fullres', 'pxl_col_in_fullres']].to_numpy()
        adata.obs.drop(columns=['pxl_roAdgre5w_in_fullres', 'pxl_col_in_fullres'], inplace=True)

        if source_image_path is not None:
            source_image_path = str(Path(source_image_path).resolve())
            adata.uns["spatial"][library_id]["metadata"]["source_image_path"] = source_image_path

    return adata


# In[25]:


import scanpy as sc
import os
import json
import pandas as pd

base_path = '/media/tombor/Helios_scStorage/LukasZ/Spatial/Yamada_mouse_AMI/spatial_only/'
samples = [d for d in os.listdir(base_path) if d.startswith('MI_')]

adata_samples = {}

for sample in samples:
    sample_path = os.path.join(base_path, sample)
    data_path = os.path.join(sample_path, 'filtered_feature_bc_matrix')
    
    # Find the prefix for this directory
    mtx_files = [f for f in os.listdir(data_path) if f.endswith('_matrix.mtx.gz')]
    if not mtx_files:
        print(f"Matrix file not found for sample {sample}. Skipping...")
        continue
    prefix = mtx_files[0].split('_matrix.mtx.gz')[0]
    
    # Load the gzipped count data with the detected prefix
    adata_1 = sc.read_10x_mtx(data_path, var_names='gene_symbols', cache=True, prefix=prefix+'_')
    # Load spatial data
    adata = process_visium_data(adata_1, path=sample_path)
    
    # Save this adata to the dictionary using sample name as key
    adata_samples[sample] = adata

# Now, adata_samples is a dictionary where keys are sample names and values are corresponding adata objects


# In[26]:


output_directory = '/media/tombor/Helios_scStorage/LukasZ/Spatial/Yamada_mouse_AMI/Analysis/scanpy/objects/'
Adgre5
# Create the output directory if it doesn't exist
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

for sample_name, adata in adata_samples.items():
    # define the save path
    save_path = os.path.join(output_directory, f"{sample_name}.h5ad")
    
    # save the AnnData object
    adata.write(save_path)

print(f"All adata samples saved to {output_directory}")

