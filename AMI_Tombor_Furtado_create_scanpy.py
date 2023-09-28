#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import scvelo as scv
import scanpy
import os
import pandas as pd
import glob
import numpy as np
import scanpy as sc
import scanorama
import re


scv.set_figure_params()


# In[ ]:


looms=list()
names=list()
mypath="/media/tombor/Helios_scStorage/Lukas/Analysis/EHT/Velocity/out/Young"

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
    samplename = re.search("/media/tombor/Helios_scStorage/Lukas/Analysis/EHT/Velocity/out/Young/(.+?)_", n).group(1)
    print(samplename)
    samplenames.append(samplename)


# In[ ]:


anno=pd.read_csv("/media/tombor/Helios_scStorage/Lukas/Analysis/EHT/Velocity/Anno.csv", sep = ",", index_col=0)
print(anno)
for i in range(len(looms)):
    looms[i].obs["Sample"] = samplenames[i]
    looms[i].obs["Timepoint"] = anno["Timepoint"][samplenames[i]]
    looms[i].obs["Age"] = anno["Age"][samplenames[i]]
    looms[i].obs["Batch"] = anno["Batch"][samplenames[i]]


# In[ ]:


processed_adatas=list()

for adata in looms:
    print(adata.shape)
    print("Filter cells")
    sc.pp.filter_cells(adata, min_genes=200)
    print("MT-filter")
    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    print("log1p")
    sc.pp.log1p(adata)
    print("HVG calc")
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    print("Regression")
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'], n_jobs = 12)
    print("Scale")
    sc.pp.scale(adata, max_value=10)
    print("save to list...")
    processed_adatas.append(adata)


# In[ ]:


adata = adata.concatenate(processed_adatas, join = "inner")

save_file = './adata_merge.h5ad'
adata.write_h5ad(save_file)

