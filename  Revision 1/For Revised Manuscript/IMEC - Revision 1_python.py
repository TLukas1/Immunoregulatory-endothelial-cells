# Import packages
import scanpy as sc
import matplotlib.pyplot as plt
plt.switch_backend('TkAgg')
plt.ion()

# Path where all figures/plots are saved
savepath = '/Users/tilllautenschlager/Library/CloudStorage/OneDrive-JohannWolfgangGoetheUniversität/Uni/PhD/Bioinformatic/24_Lukas IMEC/Revision Lukas/figures'

#Import ECs_Myeloid dataset
adata = sc.read_h5ad('/Users/tilllautenschlager/Library/CloudStorage/OneDrive-JohannWolfgangGoetheUniversität/Uni/PhD/Bioinformatic/24_Lukas IMEC/IMEC paper submission Sept 23/data/scanpy_objects/scanpy_ECs_Myeloid.h5ad')

#Define genes of Ifngr, Tgfb and Il1 receptors
genes = [ "Ifngr1", "Ifngr2", "Tgfbr1", "Il1r1", "Il1r2"]

# Subset the endothelial cell cluster
ECs = adata[adata.obs["leiden"] == "Endothelial cells"].copy()

#Plot and save dotplot (Supplement Figure 8C)
sc.pl.dotplot(ECs, genes, groupby="Timepoint", swap_axes=True,save="ECs_cytokine receptor.svg" )