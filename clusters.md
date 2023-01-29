
Run Principal Component (PC) Analysis and determine the number of PC to include in the downstream analysis:

```R
seurat <- RunPCA(seurat, nfeatures.print = 5)

# Calculate where the principal components start to elbow by taking the larger value of
# - The point where the principal components only contribute 5% of standard deviation and the principal components 
#   cumulatively contribute 90% of the standard deviation
# - The point where the percent change in variation between the consecutive PCs is less than 0.1%.

pct <- seurat[["pca"]]@stdev / sum(seurat[["pca"]]@stdev) * 100
cumulative <- cumsum(pct)
co1 <- which(cumulative > 90 & pct < 5)[1]
co2 <-
  sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
nPCs <- min(co1, co2)
plot_df <-
  data.frame(pct = pct,
             cumu = cumulative,
             rank = 1:length(pct))

# Visualise the PCs to use
ggplot(plot_df, aes(cumulative, pct, label = rank, color = rank > nPCs)) + geom_text() + theme_bw()
```

Find neighbours and estimate clusters, then run UMAP for dimensinal reduction:

```R
seurat <- FindNeighbors(seurat, dims = 1:nPCs, reduction = 'pca')
gc(full = T)

# Determine the resolution for finding clusters using cutree package
library(clustree)


# Find clusters for a range of resolutions 
seuratObj <- seurat
resolutionRange <- seq(from = 0, to = 1.2, by = 0.1)

# Find clusters using a range of resolutions
seuratObj <- FindClusters(object = seuratObj, resolution = resolutionRange)

# resolution prefix in Seurat Object
head(seuratObj[[]])

# Build clustree
clustree(seuratObj, prefix = "integrated_snn_res.")
rm(seuratObj); gc(full = T)

sessionInfo()

# Find clusters for the integrated Seurat Object
seurat <-
  FindClusters(seurat, resolution = 0.7) # determined from the chunk above

# UMAP projection
seurat <- RunUMAP(seurat, dims = 1:nPCs, reduction = 'pca')
```

Reorder the conditions for the desired order and save the processed Seurat Object:

```R
seurat$condition <-
  factor(seurat$condition, levels = c('uninfected', 'infected'))
seurat$orig.ident <-
  factor(
    seurat$orig.ident,
    levels = c('uninfected_1', 'uninfected_2', 'infected_1', 'infected_2')
  )

saveRDS(seurat, file = paste0(date, '_integrated_seurat_clusters.rds'))
gc(full = T)
```
