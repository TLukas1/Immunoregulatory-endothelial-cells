# Import .h5a to Seurat

file <- "/media/tombor/Helios_scStorage/Lukas/Analysis/EHT/Scripts/adata_merge.h5ad"
Sys.setenv(RETICULATE_PYTHON="/usr/bin/python3")
library(reticulate)
library(anndata)
library(Seurat)

ad <- import("anndata", convert = FALSE)
adata <- ad$read_h5ad(file)

raw.data.matrix <- tryCatch(
  expr = t(as.matrix(py_to_r(adata$layers$`_data`["matrix"]))),
  error = function(e) {
    stop("No adata.raw.X in provided adata. please make sure adata has adata.raw.X when you tyr to turn it to `seurat object`")
  })

rownames(x = raw.data.matrix) <- as.character(py_to_r(adata$var_names))
colnames(x = raw.data.matrix) <- as.character(py_to_r(adata$obs_names))

meta.data <- as.data.frame(py_to_r(adata$obs))

Myo_concat <- CreateSeuratObject(counts = raw.data.matrix, meta.data = meta.data)

save(Myo_concat, file = "/media/tombor/Helios_scStorage/Lukas/Analysis/EHT/Combined_Seurats/220517_Seurat_of_Anndata_concat.Rds")


###### Anndata 3


file <- "/media/tombor/Helios_scStorage/Lukas/Analysis/EHT/Combined_Seurats/scanpy_ECs_Myeloid.h5ad"
Sys.setenv(RETICULATE_PYTHON="/usr/bin/python3")
library(reticulate)
library(anndata)
library(Seurat)

ad <- import("anndata", convert = FALSE)
adata <- ad$read_h5ad(file)

raw.data.matrix <- tryCatch(
  expr = t(py_to_r(adata$raw$X)),
  error = function(e) {
    stop("No adata.raw.X in provided adata. please make sure adata has adata.raw.X when you tyr to turn it to `seurat object`")
  })
######### Github function #######

suppressMessages(library(Seurat))# the version of Seurat <3.0
suppressMessages(library(reticulate))
#turn category obs in adata into string
change_obs=function(adata,exclude=c("n_genes","n_counts")){
  #@exclude: the columns names
  n_col=as.numeric(as.character(adata$obs$columns$shape[0]))
  tmp=adata$obs$columns
  for (i in seq_len(n_col)){
    cur_colname=as.character(tmp[i-1])
    if (!cur_colname %in% exclude){
      adata$obs[[cur_colname]]=adata$obs[[cur_colname]]$astype("str")
    }
  }
  return(adata)
}


Convert_from_anndata_to_seurat=function(from=adata,X.slot="scale.data",raw.X.slot="logcount.data"){
  if(!py_module_available("anndata")) {
    stop("Please install the anndata python module")
  }
  stopifnot(X.slot%in%c("scale.data","normlizecount.data"))
  stopifnot(raw.X.slot%in%c("count.data","logcount.data"))
  data.matrix=tryCatch(
    expr=t(py_to_r(from$X)),
    error=function(e){
      stop("No adata.X in provided adata. If Both adata.X and adata.raw.X are None")
    }
  )
  rownames(data.matrix)<-rownames(py_to_r(from$var))
  colnames(data.matrix)<-rownames(py_to_r(from$obs))
  if (X.slot=="normalizecount.data") X.slot="data"
  
  raw.data.matrix <- tryCatch(
    expr = t(py_to_r(from$raw$X)),
    error = function(e) {
      stop("No adata.raw.X in provided adata. please make sure adata have adata.raw.X when you tyr to turn it to `seurat object`")
    }
  )
  if(raw.X.slot=="logcount.data"){
    raw.data.matrix=expm1(raw.data.matrix)
  }
  
  rownames(x = raw.data.matrix) <- rownames(x = py_to_r(from$raw$var))
  colnames(x = raw.data.matrix) <- rownames(x =py_to_r(from$obs))
  #get meta.data
  meta.data=py_to_r(from$obs)
  if ("nUMI" %in% colnames(x = meta.data)) {
    colnames(x = meta.data) <- gsub(
      pattern = "nUMI",
      replacement = "nUMI_ori",
      x = colnames(x = meta.data)
    )
  }
  if ("nGene" %in% colnames(x = meta.data)) {
    colnames(x = meta.data) <- gsub(
      pattern = "nGene",
      replacement = "nGene_ori",
      x = colnames(x = meta.data)
    )
  }
  seurat.object <- CreateSeuratObject(counts = raw.data.matrix,meta.data = meta.data)
  seurat.object <- SetAssayData(
    object = seurat.object,
    assay = "RNA",
    slot = X.slot,
    new.data = data.matrix
  )
  #deal with obsm fields that are not dimensional reductions, or have different name structures
  x1=py_to_r(from$obsm$keys())
  drs<-unlist(strsplit(gsub(".{1,50}:|\\s|)",replacement = "",x = x1),split = ","))
  for (dr in drs) {
    em <- py_to_r(from$obsm)
    dr.embed <- em[[eval(dr)]]
    dr.name <- sub(pattern="X_",replacement="",x=dr)
    if (is.na(dr.name)|dr.name=="") {
      dr.name <- dr
    }
    dr.dict <- list(tSNE_ = "tsne", PC = "pca")
    if (dr.name %in% dr.dict) {
      dr.key <- names(x = which(x = dr.dict == dr.name))
    } else {
      dr.key <- toupper(x = dr.name)
    }
    colnames(x = dr.embed) <- paste0(dr.key, 1:ncol(x = dr.embed))
    rownames(x = dr.embed) <- rownames(seurat.object@meta.data)
    seurat.object@reductions[[dr.name]] <- dr.embed
  }
  return(seurat.object)
}

adata=change_obs(adata,exclude=c("n_genes","n_counts"))# columns to exclude
EcHe=Convert_from_anndata_to_seurat(adata,raw.X.slot = "count.data")

save(EcHe, file = "/home/tombor/Desktop/220514_Seurat_of_Anndata3.Rds")


###


