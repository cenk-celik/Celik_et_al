library(Seurat)

rm(list = ls(all.names = T))

# clean R environment
seurat <- readRDS('seurat.rds')

# split SeuratObject by 'condition'
split.list <- SplitObject(seurat, split.by = 'condition')

for (i in 1:length(split.list)) {
  condition <- names(split.list[i])
  seurat <- unlist(split.list[[i]])
  
  # save metadata table:
  seurat$barcode <- colnames(seurat)
  seurat$UMAP_1 <- seurat@reductions$umap@cell.embeddings[, 1]
  seurat$UMAP_2 <- seurat@reductions$umap@cell.embeddings[, 2]
  write.csv(
    seurat@meta.data,
    file = paste0(condition, '_metadata.csv'),
    quote = F,
    row.names = F
  )
  
  # write expression counts matrix
  library(Matrix)
  counts_matrix <-
    GetAssayData(seurat, assay = 'RNA', slot = 'counts')
  writeMM(counts_matrix, file = paste0(condition, '_counts.mtx'))
  
  # write dimesnionality reduction matrix, in this example case pca matrix
  write.csv(
    seurat@reductions$pca@cell.embeddings,
    file = paste0(condition, '_pca.csv'),
    quote = F,
    row.names = F
  )
  
  # write gene names
  write.table(
    data.frame('gene' = rownames(counts_matrix)),
    file = paste0(condition, '_gene_names.csv'),
    quote = F,
    row.names = F,
    col.names = F
  )
}

sessionInfo()
