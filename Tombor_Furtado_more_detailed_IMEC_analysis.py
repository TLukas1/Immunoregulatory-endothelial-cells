#!/usr/bin/env python
# coding: utf-8

# In[97]:


import os
import pandas as pd
import glob
import numpy as np
import scanpy as sc
from matplotlib.pyplot import rc_context
import seaborn as sns
import matplotlib.pyplot as plt


adata3 = sc.read(filename= "/media/tombor/Helios_scStorage/Lukas/Analysis/EHT/Combined_Seurats/scanpy_ECs_Myeloid.h5ad")


# In[117]:


ECs = adata3[adata3.obs["subcluster"].isin(["IMEC (16)", "Venous EC (4)", "Artery EC (8)", "Artery EC (0)", "Lymphatic EC (12)", "Artery EC (17)"])]
ECs_2 = adata3[adata3.obs["subcluster"].isin(["IMEC (16)", "Venous EC (4)", "Artery EC (8)", "Lymphatic EC (12)", "Artery EC (17)", "Inflamm. Monocyte (3)"])]


# In[119]:


# New Figure 1C
sc.pl.dotplot(ECs, ["Pecam1", "Cdh5", "Cd36", "Fabp4", "Eng", "Esam", "Tie1", "Kdr", "Ptprc", "Runx1", "Irf8", "Cd68", "Cxcl16", "Ccl7", "Adgre5", "Cd83", "H2-Aa", "H2-Ab1", "H2-DMa", "H2-DMb1", "H2-Eb1", "Cnn1", "Fn1", "Serpine1", "Tnc"], groupby='subcluster', cmap = "Reds", dendrogram = False, save = "dotplot_new_Figure1C.pdf")
sc.pl.stacked_violin(ECs, ["Pecam1", "Cdh5", "Cd36", "Fabp4", "Eng", "Esam", "Tie1", "Kdr", "Ptprc", "Itgam", "Runx1", "Irf8", "Cd68", "Cxcl16", "Ccl7", "Adgre5", "Cd83", "H2-Aa", "H2-Ab1", "H2-DMa", "H2-DMb1", "H2-Eb1", "Cnn1", "Fn1", "Serpine1", "Tnc"], groupby='subcluster', save = "vln_new_Figure1C.pdf")


# # Receiving and downstream signals of IMEC

# In[121]:


sc.pl.dotplot(ECs, ["Ifngr1", "Ifngr2", "Irf2", "Irf8", "Ccr1", "Ccr8", "Tgfbr1", "Zeb2", "Mmp9", "Il1r1", "Il1r2","Il1b", "Cxcl2", "Tnfrsf1a", "Tnfrsf1b", "Nfkbia","Tnf", "Cxcl10", "Ccl2"], groupby='subcluster', cmap = "Reds", dendrogram = False, save = "dotplot_Receiving_and_downstream.pdf")
sc.pl.dotplot(ECs, ["Runx1", "Ptprc", "Irf8", "Cd68", "Cxcl16", "Ccl7", "Adgre5", "Cd83", "Cnn1", "Fn1", "Serpine1", "Tnc", "Birc5", "Cdk1", "Tyms"], groupby='subcluster', cmap = "Reds", dendrogram = True, save = "dotplot_Transition_vs_activation.pdf")
sc.pl.dotplot(ECs, ["Runx1", "Ptprc", "Icam1", "Vcam1", "Nos3", "Sele", "Selp", "Vwf"], groupby='subcluster', cmap = "Reds", dendrogram = True, save = "dotplot_Transition_vs_activation_2.pdf")
sc.pl.dotplot(ECs_2, ["H2-Aa", "H2-Ab1", "H2-DMa", "H2-DMb1", "H2-Eb1"], groupby='subcluster', cmap = "Reds", dendrogram = True, save = "Antigen_presentation.pdf")
sc.pl.dotplot(ECs_2, ["Cd80", "Cd86",  "Icosl", "Tnfsf4", "Cd40", "Cd274", "Pdcd1lg2"], groupby='subcluster', cmap = "Reds", dendrogram = True, save = "Co-receptors.pdf")


# In[ ]:





# In[78]:


anno=pd.read_csv("Ifng genes.csv", sep = ",", index_col=0)
sc.pl.matrixplot(ECs, anno["Symbol"], groupby="subcluster")


# In[122]:


def coexpressing_fraction_or_conditions(adata, gene_sets, group_key):
    # Binary matrix where 1 indicates expression and 0 indicates no expression
    binary_matrix = (adata.X > 0).astype(int)

    # Convert each gene set into a boolean matrix where each column represents a gene
    bool_matrices = []
    for gene_set in gene_sets:
        indices = [adata.var_names.get_loc(gene) for gene in gene_set]
        bool_matrices.append(binary_matrix[:, indices].sum(axis=1) > 0)

    # Combine the boolean matrices using AND condition
    combined_matrix = bool_matrices[0]
    for matrix in bool_matrices[1:]:
        combined_matrix = combined_matrix & matrix

    # Get unique groups
    groups = adata.obs[group_key].cat.categories

    results = {}
    for group in groups:
        subset = combined_matrix[adata.obs[group_key] == group]
        fraction = subset.mean()
        results[group] = fraction

    return results



# In[123]:


gene_combinations = [
    [["Irf2", "Irf8"], ["H2-Aa", "H2-DMa", "H2-Ab1", "H2-DMb1", "H2-Eb1"]],
    [["Runx1", "Ptprc"], ["H2-Aa", "H2-DMa", "H2-Ab1", "H2-DMb1", "H2-Eb1"]],
    [["Runx1", "Ptprc"], ["Irf2", "Irf8"]]
]

group_key = "subcluster"

# Dictionary to store the results
results = {}

for combo in gene_combinations:
    combo_name = ' AND '.join([' OR '.join(subcombo) for subcombo in combo])
    results[combo_name] = coexpressing_fraction_or_conditions(ECs, combo, group_key)

# Convert dictionary to DataFrame
df = pd.DataFrame(results)
print(df)


# In[124]:


plt.figure(figsize=(3, 3))

# Create a 2D array of strings for the annotations with percentage format
percentage_annotations = df.applymap(lambda x: f"{x*100:.2f}%").values

ax = sns.heatmap(df, annot=percentage_annotations, fmt='', cmap="viridis", cbar_kws={'label': 'Fraction of Cells'})
# Adjusting the tick labels for better visualization
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')  # `ha` is for horizontal alignment

plt.title("Co-expression fractions")
#plt.tight_layout()

# Save to PDF
plt.savefig("coexpression_heatmap_Ifng_HLA_Runx1.pdf", format='pdf')
plt.show()


# In[61]:


anno=pd.read_csv("Ifng genes.csv", sep = ",", index_col=0)
print(anno["Symbol"])


# In[ ]:


anno=pd.read_csv("Ifng genes.csv", sep = ",", index_col=0)
sc.pl.heatmap(ECs, anno["Symbol"], groupby='subcluster', cmap = "Reds")


# In[87]:


markers = {
"MHC-II": ["H2-Aa", "H2-Ab1", "H2-D1", "H2-DMa", "H2-DMb1", "H2-DMb2", "H2-Eb1", "H2-K1", "H2-Ke6", "H2-M2", "H2-M3", "H2-M5", "H2-M9", "H2-Oa", "H2-Ob", "H2-Q1", "H2-Q10", "H2-Q2", "H2-Q4", "H2-Q6", "H2-Q7", "H2-T22", "H2-T23", "H2-T24"]}

sc.pl.dotplot(ECs, markers, groupby='Timepoint', dendrogram = False,
              save = "H2_all_time.pdf")

IMEC = ECs.obs["subcluster"].isin(["IMEC (16)"])
sc.pl.dotplot(IMEC, markers, groupby='Timepoint', dendrogram = False,
              save = "H2_all_time_IMEC.pdf")

markers = {"Ifng": ["Irf8","Ifngr1"],
           "Cytokines": ["Ccl2", "Ccl3", "Ccl4", "Ccl6", "Ccl9", "Ccl12"],
           "Interleukin 1": ["Il1b", "Il1r1", "Il1r2"], 
           "TGF-b" : ["Zeb2", "Tgfb1", "S100a4", "Mmp9"],
           "Tnf-a" : ["Tnfrsf1a", "Nfkbia","Tnfrsf1b"]}
#sc.pl.dotplot(ECs, markers, groupby='subcluster', dendrogram = True, save = "Newer_2a.pdf")


# In[6]:


sc.tl.rank_genes_groups(IMECs, 'Timepoint', method='wilcoxon', key_added = "wilcoxon")
sc.pl.rank_genes_groups(IMECs, n_genes=25, sharey=False, key="wilcoxon")

