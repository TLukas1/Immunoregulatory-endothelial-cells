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

#%% Cytokine expression in different myeloid subclusters
adata = sc.read_h5ad('/Users/tilllautenschlager/Library/CloudStorage/OneDrive-JohannWolfgangGoetheUniversität/Uni/PhD/Bioinformatic/24_Lukas IMEC/IMEC paper submission Sept 23/data/scanpy_objects/scanpy_ECs_Myeloid.h5ad')

#Define genes of interest
selected_genes = ["Ifng","Tgfb2", "Il1b"]

# List of cell types to subset
cell_types = adata.obs["subcluster"].unique()

adata_subsets = {}

#Create a new dotplot for each myeloid subcluster in the dataset
for cell_type in cell_types:
    adata_subsets[cell_type] = adata[adata.obs["subcluster"] == cell_type].copy()

for cell_type in cell_types:
    # Generate the dotplot for the current cell type using the 'all_genes'
    sc.pl.dotplot(
        adata_subsets[cell_type],
        var_names=genes["selected_genes"],
        groupby="Timepoint",
        cmap="Reds",
        swap_axes=True,
        title=f"Dotplot of {cell_type}",
        show=False
    )
    print(cell_type)

    # Rotate x-axis labels by 90 degrees
    plt.xticks(rotation=90)

 # Save the plot to the specified path
    plt.savefig(os.path.join(savepath, f"{cell_type}_dotplot.png"), dpi=600, bbox_inches='tight')
    plt.close()