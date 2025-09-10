
library(Seurat)
library(SeuratDisk)
library(reticulate)
# library(anndata)
ad <- reticulate::import("anndata")
adata <- ad$read_h5ad(pbmc10kmono)

adata
# AnnData object with n_obs × n_vars = 3782 × 13483
#     obs: 'orig.ident', 'n_genes', 'celltypeL0'
#     var: 'features', 'n_cells'

adata$T$X[1:5,1:5]
# 5 x 5 sparse Matrix of class "dgRMatrix"
#            AAACAGCCATCCAGGT-1 AAACCAACACAATGCC-1 AAACCAACAGGAACTG-1
# AL627309.5                  .                  .                  .
# LINC01409                   .                  .                  .
# LINC01128                   .                  .                  .
# LINC00115                   1                  .                  .
# NOC2L                       .                  .                  2
#            AAACCAACATAATCCG-1 AAACCAACATTGTGCA-1
# AL627309.5                  1                  .
# LINC01409                   1                  .
# LINC01128                   .                  .
# LINC00115                   .                  .
# NOC2L                       .                  .

# 细胞和基因的元数据
cell_meta <- as.data.frame(adata$obs)
gene_meta <- as.data.frame(adata$var)
# 创建Seurat对象 
seurat_obj <- Seurat::CreateSeuratObject(
  counts = adata$T$X,   ## 这里需要转置，因为seurat中以cell作为列，gene作为行
  meta.data = cell_meta,
  assay = "RNA"
)
# 添加基因注释（如基因类型）
seurat_obj[["RNA"]]@meta.features <- gene_meta