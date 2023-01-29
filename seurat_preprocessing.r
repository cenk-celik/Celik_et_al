home <- system("echo ${HOME}/", intern = TRUE)
main <- paste0(home, "nmrc/")

#----Load libraries & tweak processors----
library(Seurat)
library(celldex)
library(patchwork)
library(dplyr)
library(ggplot2)
library(future)

plan(multisession, workers = 4)
options(future.globals.maxSize = 10 * 1024 ^ 3,
        future.seed = NULL)

##----Read Mtx files----
data.list <- list()

data.list[[1]] <-
  ReadMtx(mtx = './data/uninfected_1/outs/filtered_feature_bc_matrix/matrix.mtx.gz',
          cells = './data/uninfected_1/outs/filtered_feature_bc_matrix/barcodes.tsv.gz',
          features = './data/uninfected_1/outs/filtered_feature_bc_matrix/features.tsv.gz')

data.list[[2]] <-
  ReadMtx(mtx = './data/uninfected_2/outs/filtered_feature_bc_matrix/matrix.mtx.gz',
          cells = './data/uninfected_2/outs/filtered_feature_bc_matrix/barcodes.tsv.gz',
          features = './data/uninfected_2/outs/filtered_feature_bc_matrix/features.tsv.gz')

data.list[[3]] <-
  ReadMtx(mtx = './data/infected_1/outs/filtered_feature_bc_matrix/matrix.mtx.gz',
          cells = './data/infected_1/outs/filtered_feature_bc_matrix/barcodes.tsv.gz',
          features = './data/infected_1/outs/filtered_feature_bc_matrix/features.tsv.gz')

data.list[[4]] <-
  ReadMtx(mtx = './data/infected_2/outs/filtered_feature_bc_matrix/matrix.mtx.gz',
          cells = './data/infected_2/outs/filtered_feature_bc_matrix/barcodes.tsv.gz',
          features = './data/infected_2/outs/filtered_feature_bc_matrix/features.tsv.gz')

seurat.list <- list()
samples <-
  c('uninfected_1', 'uninfected_2', 'infected_1', 'infected_2')
conditions <- c('uninfected', 'uninfected', 'infected', 'infected')

##----Create Seurat Objects----
for (i in 1:length(data.list)) {
  seurat.list[[i]] <-
    CreateSeuratObject(
      counts = data.list[[i]],
      min.cells = 3,
      min.features = 200,
      project = samples[i]
    )
  seurat.list[[i]][["samples"]] <-  samples[i]
  seurat.list[[i]][['percent.mt']] <-
    PercentageFeatureSet(seurat.list[[i]], pattern = '^mt-')
  seurat.list[[i]]@meta.data$condition <- conditions[i]
}

rm(data.list)

i = 4
p1 <-
  VlnPlot(
    seurat.list[[i]],
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3,
    pt.size = 0.1
  ) +
  theme(legend.position = 'none', axis.title.x = element_blank())
p2 <-
  FeatureScatter(seurat.list[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt") +
  theme(legend.position = 'none', axis.title.x = element_blank())
p3 <-
  FeatureScatter(seurat.list[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  theme(legend.position = 'none', axis.title.x = element_blank())
cowplot::plot_grid(p1, p2 + p3, nrow = 2)

#determine percent.mt to remove in mitochondria_quality.r
seurat.list[[1]] <-
  subset(seurat.list[[1]],
         subset = nFeature_RNA > 250 &
           nFeature_RNA < 7500 & percent.mt < mito[[1]] & nCount_RNA < 75000)
seurat.list[[2]] <-
  subset(seurat.list[[2]],
         subset = nFeature_RNA > 200 &
           nFeature_RNA < 7500 & percent.mt < mito[[2]] & nCount_RNA < 60000)
seurat.list[[3]] <-
  subset(seurat.list[[3]],
         subset = nFeature_RNA > 1000 &
           nFeature_RNA < 7000 & percent.mt < mito[[3]] & nCount_RNA < 50000)
seurat.list[[4]] <-
  subset(seurat.list[[4]],
         subset = nFeature_RNA > 200 &
           nFeature_RNA < 7500 & percent.mt < mito[[4]] & nCount_RNA < 60000)

##----Remove doublets----
library(scDblFinder)
library(SingleCellExperiment)

seurat.list.sce <- lapply(
  X = seurat.list,
  FUN = function(x) {
    x <- as.SingleCellExperiment(x)
    x <- scDblFinder(x, samples = 'orig.ident', clusters = F)
  }
)

for (i in 1:length(seurat.list)) {
  seurat.list[[i]] <-
    as.Seurat(seurat.list.sce[[i]], counts = "counts", data = "logcounts")
  print(table(seurat.list.sce[[i]]$scDblFinder.class))
}

seurat.list <- lapply(seurat.list, function(x) {
  x <- subset(x, subset = scDblFinder.class == 'singlet')
})

rm(seurat.list.sce)
gc(full = T)

##----Cell cycles----
library(nichenetr)
geneList <- list()
geneList <- lapply(
  cc.genes.updated.2019,
  FUN = function(x) {
    x <- convert_human_to_mouse_symbols(x)
  }
)

s.genes <- geneList$s.genes
g2m.genes <- geneList$g2m.genes

for (i in 1:length(seurat.list)) {
  seurat.list[[i]] <-
    CellCycleScoring(
      seurat.list[[i]],
      s.features = s.genes,
      g2m.features = g2m.genes,
      set.ident = FALSE
    )
  seurat.list[[i]]$CC.Difference <-
    seurat.list[[i]]$S.Score - seurat.list[[i]]$G2M.Score
  print(VlnPlot(
    seurat.list[[i]],
    features = c("S.Score", "G2M.Score"),
    group.by = 'orig.ident',
    ncol = 2,
    pt.size = 0
  ))
}

gc(full = T)

seurat.list <- lapply(
  X = seurat.list,
  FUN = function(x) {
    x <-
      SCTransform(x,
                  vars.to.regress = c('percent.mt', 'CC.Difference'),
                  method = 'glmGamPoi')
  }
)

for (i in 1:length(seurat.list)) {
  seurat.list[[i]]$scDblFinder.sample <- NULL
  seurat.list[[i]]$scDblFinder.class <- NULL
  seurat.list[[i]]$scDblFinder.weighted <- NULL
  seurat.list[[i]]$scDblFinder.score <- NULL
  seurat.list[[i]]$scDblFinder.cxds_score <- NULL
  seurat.list[[i]]$CC.Difference <- NULL
  seurat.list[[i]]$S.Score <- NULL
  seurat.list[[i]]$G2M.Score <- NULL
}

gc(full = T)

##----Integrate datasets----
features <-
  SelectIntegrationFeatures(object.list = seurat.list, nfeatures = 3000)
seurat.list <-
  PrepSCTIntegration(object.list = seurat.list, anchor.features = features)
gc(full = T)
anchors <-
  FindIntegrationAnchors(
    object.list = seurat.list,
    normalization.method = "SCT",
    anchor.features = features
  )
gc(full = T)
seurat <-
  IntegrateData(anchorset = anchors, normalization.method = "SCT")

dim(seurat)

rm(seurat.list)
rm(anchors)
date <- format(Sys.time(), '%Y%m%d')
saveRDS(seurat, file = paste0(date, '_integrated_seurat_raw.rds'))

gc(full = T)

##----PCA----
seurat <- RunPCA(seurat, nfeatures.print = 5)

#Determine the number of PCs
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

ggplot(plot_df, aes(cumulative, pct, label = rank, color = rank > nPCs)) + geom_text() + theme_bw()

seurat <- FindNeighbors(seurat, dims = 1:nPCs, reduction = 'pca')
gc(full = T)
seurat <-
  FindClusters(seurat, resolution = 0.7) #determine in cluster_determination.R
seurat <- RunUMAP(seurat, dims = 1:nPCs, reduction = 'pca')

seurat$condition <-
  factor(seurat$condition, levels = c('uninfected', 'infected'))
seurat$orig.ident <-
  factor(
    seurat$orig.ident,
    levels = c('uninfected_1', 'uninfected_2', 'infected_1', 'infected_2')
  )

sample_colours <- c('#FFCDB2', '#FFB4A2', '#6D6875', '#B5838D')
condition_colours <- c('#E5989B', '#404040')
samples <-
  DimPlot(
    seurat,
    reduction = 'umap',
    group.by = 'orig.ident',
    pt.size = 0.5,
    cols = sample_colours
  ) + NoAxes() + ggtitle('Samples') + NoLegend()
pdf('./plots/integrated_seurat_orig.ident.pdf',
    height = 5,
    width = 5)
plot(samples)
dev.off()

conditions <-
  DimPlot(
    seurat,
    reduction = 'umap',
    group.by = 'condition',
    pt.size = 0.5,
    cols = condition_colours
  ) + NoAxes() + ggtitle('Conditions') + NoLegend()
pdf('./plots/integrated_seurat_condition.pdf',
    height = 5,
    width = 5)
plot(conditions)
dev.off()

clusters <-
  DimPlot(
    seurat,
    reduction = 'umap',
    label = T,
    group.by = 'seurat_clusters',
    pt.size = 0.5
  ) +
  NoAxes() + ggtitle('') + NoLegend()
pdf('./plots/integrated_seurat_clusters.pdf',
    height = 5,
    width = 5)
plot(clusters)
dev.off()

saveRDS(seurat, file = paste0(date, '_integrated_seurat_clusters.rds'))
gc(full = T)

##----Identify Markers----
logfc <- log2(1.5)

seurat <- PrepSCTFindMarkers(seurat)

seurat.markers <-
  FindAllMarkers(
    seurat,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = logfc
  )
seurat.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) %>% 
  write.csv(paste0(date, '_seurat.markers.csv'))

saveRDS(seurat.markers, file = paste0(date, '_seurat.markers.rds'))

##----Annotation----
library(HGNChelper)
library(openxlsx)
source("gene_sets_prepare.r")
source("sctype_score_.R")

db_ <- "Joost_2ndLevel.xlsx"
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

celltype_2nd_separated <-
  DimPlot(
    seurat,
    reduction = "umap",
    split.by = 'condition',
    label = F,
    repel = TRUE,
    group.by = 'celltype_2nd',
    pt.size = 0.5
  ) +
  ggtitle('') + NoAxes() + NoLegend()
pdf(
  './plots/integrated_seurat_celltype_2nd_separated.pdf',
  height = 5,
  width = 9
)
plot(celltype_2nd_separated)
dev.off()

saveRDS(seurat,
        file = paste0(date, '_integrated_seurat_main_cell_types.rds'))

##----Add idents manually----
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
clusters_annotation <-
  DimPlot(
    seurat,
    reduction = "umap",
    split.by = 'condition',
    label = F,
    repel = TRUE,
    pt.size = 0.5
  ) +
  ggtitle('') + ggtitle('') + NoAxes()
pdf(
  './plots/integrated_seurat_clusters_annotation.pdf',
  height = 5,
  width = 11
)
plot(clusters_annotation)
dev.off()

#saveRDS(seurat, file = paste0(date, '_integrated_seurat_idents_annotated.rds'))

##----Split by condition----
split.list <- SplitObject(seurat, split.by = "condition")

uninfected <- split.list[[1]]
infected <- split.list[[2]]

saveRDS(uninfected, file = paste0(date, '_integrated_uninfected.rds'))
saveRDS(infected, file = paste0(date, '_integrated_infected.rds'))

rm(split.list)
gc()

sessionInfo()