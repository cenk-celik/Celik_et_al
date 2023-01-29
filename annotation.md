Load the required packages for cluster annotation:
```R
library(HGNChelper)
library(openxlsx)

# the following script was adapted from https://github.com/IanevskiAleksandr/sc-type
source("gene_sets_prepare.r")
source("sctype_score_.R")

db_ <- "cell_annotation_1stLevel.xlsx" #for main cell annotation
# db_ <- "cell_annotation_2ndLevel.xlsx" #for 2nd level cell annotation

tissue <- "Skin"

gs_list <-  gene_sets_prepare(db_, tissue)
es.max <-
  sctype_score(
    scRNAseqData = seurat[["integrated"]]@scale.data,
    scaled = TRUE,
    gs = gs_list$gs_positive,
    gs2 = gs_list$gs_negative
  )

cL_resutls <-
  do.call("rbind", lapply(unique(seurat@meta.data$seurat_clusters), function(cl) {
    es.max.cl <-
      sort(rowSums(es.max[, rownames(seurat@meta.data[seurat@meta.data$seurat_clusters ==
                                                        cl,])]), decreasing = !0)
    head(data.frame(
      cluster = cl,
      type = names(es.max.cl),
      scores = es.max.cl,
      ncells = sum(seurat@meta.data$seurat_clusters == cl)
    ),
    10)
  }))

sctype_scores <-  cL_resutls %>%
  group_by(cluster) %>%
  top_n(n = 1, wt = scores)

sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells /
                     4] = "Unknown"
print(sctype_scores[, 1:4])

seurat@meta.data$celltype = ""
for (j in unique(sctype_scores$cluster)) {
  cl_type = sctype_scores[sctype_scores$cluster == j, ]
  
  seurat@meta.data$celltype[seurat@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}
```

Visualise annotations:
```R
# visualise the main cell populations:
celltype <-
  DimPlot(
    seurat,
    reduction = "umap",
    label = F,
    repel = TRUE,
    group.by = 'celltype',
    pt.size = 0.5
  ) +
  ggtitle('Cell types') + NoAxes() + NoLegend()
pdf('./plots/integrated_seurat_celltype.pdf',
    height = 5,
    width = 5)
plot(celltype)
dev.off()

# visualise 2nd Level cell populations:
celltype_2nd <-
  DimPlot(
    seurat,
    reduction = "umap",
    label = F,
    repel = TRUE,
    group.by = 'celltype_2nd',
    pt.size = 0.5
  ) +
  ggtitle('Cell types - 2nd Level') + NoAxes() + NoLegend()
pdf(
  './plots/integrated_seurat_celltype_2nd.pdf',
  height = 5,
  width = 5
)
plot(celltype_2nd)
dev.off()

# save integrated Seurat Object with main and 2nd level cell annotations:
saveRDS(seurat,
        file = paste0(date, '_integrated_seurat_main_cell_types.rds'))
```

Rename ```Idents``` of the integrated Seurat Object based on the output from the script above:

```R
cluster.ids <- c(
  'MAC1',
  'FIB1',
  'FIB2',
  'BAS1',
  'SUP1',
  'HFSUP1',
  'SG',
  'NEUT1',
  'SUP2',
  'MAC2',
  'SUP3',
  'EC1',
  'EC2',
  'SUP4',
  'OB',
  'FIB3',
  'HFSUP2',
  'DC',
  'VSM',
  'NEU',
  'NEUT2',
  'BAS3',
  'TMEM',
  'BAS2'
)

names(cluster.ids) <- levels(seurat)
seurat <- RenameIdents(seurat, cluster.ids)

# save the Seurat Object with renamed Idents:
saveRDS(seurat, file = paste0(date, '_integrated_seurat_idents_annotated.rds'))
```
