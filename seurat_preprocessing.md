Set the project folder:

```R
home <- system("echo ${HOME}/", intern = TRUE)
main <- paste0(home, "nmrc/")
```

Load the libraries:

```R
library(Seurat)
library(celldex)
library(patchwork)
library(dplyr)
library(ggplot2)
library(future)

# set parallellisation parameters
plan(multisession, workers = 4)
options(future.globals.maxSize = 10 * 1024 ^ 3,
        future.seed = NULL)
```

Read data in ```filtered_feature_bc_matrix``` folder for all samples:

```R
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
```

Create sample and condition vectors for annotation for ```SeuratObjects```:

```R
seurat.list <- list()

samples <-
  c('uninfected_1', 'uninfected_2', 'infected_1', 'infected_2')
conditions <- c('uninfected', 'uninfected', 'infected', 'infected')
```

Create Seurat Objects:

```R
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

rm(data.list) # remove unnecessary objects
```

Visualise mitochondrial percentage in violin plots:

```R
i = 1 #2, 3, 4
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
```

Determine low quality cells:

```R
# Adapted from miQC package
library(miQC)
library(scater)
library(flexmix)
library(splines)
library(nichenetr)

seurat.list.sce <- lapply(X = seurat.list, FUN = function(x) {
  x <- as.SingleCellExperiment(x)
})

mito <- list() # record percentage of cells to remove from the data in a list

for (i in 1:length(seurat.list)) {
  mt_genes <- grepl("^mt-",  rownames(seurat.list.sce[[i]]))
  feature_ctrls <- list(mito = rownames(seurat.list.sce[[i]])[mt_genes])
  
  seurat.list.sce[[i]] <- addPerCellQC(seurat.list.sce[[i]], subsets = feature_ctrls)
  head(colData(seurat.list.sce[[i]]))
  
  plotMetrics(seurat.list.sce[[i]])
  model <- mixtureModel(seurat.list.sce[[i]])
  plotModel(seurat.list.sce[[i]], model)
  sce <- seurat.list.sce[[i]]
  metrics <- as.data.frame(colData(sce))
  
  intercept1 <- parameters(model, component = 1)[1]
  intercept2 <- parameters(model, component = 2)[1]
  if (intercept1 > intercept2) {
    compromised_dist <- 1
    intact_dist <- 2
  } else {
    compromised_dist <- 2
    intact_dist <- 1
  }
  
  post <- posterior(model); posterior_cutoff = 0.75
  prob_compromised <- post[, compromised_dist]
  keep <- prob_compromised <= posterior_cutoff
  
  metrics <- cbind(metrics, prob_compromised = prob_compromised, keep = keep)
  
  predictions <- fitted(model)[, intact_dist]
  metrics$intact_prediction <- predictions
  metrics[metrics$subsets_mito_percent <
            metrics$intact_prediction, ]$keep <- TRUE
  
  print(min_discard <- min(metrics[!metrics$keep, ]$subsets_mito_percent))
  mito[[i]] <- min_discard
  print(plotFiltering(seurat.list.sce[[i]], model, detected = 'detected', subsets_mito_percent = 'subset_mito_percent')+
          geom_hline(yintercept=min_discard, linetype="dashed", color = "red"))

}

rm(seurat.list.sce) # remove unnecessary objects
```

Remove low quality cells from the samples:
```R
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
```

Remove doublet GEMs:

```R
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

rm(seurat.list.sce) # remove unnecessary objects
gc(full = T)
```

Estimate cell cycle genes to regress out during the ```SCTransform```:

```R
##----Cell cycles----
library(nichenetr)
geneList <- list()

# convert Human gene symbols embedded in Seurat to mouse gene symbols:
geneList <- lapply(
  cc.genes.updated.2019,
  FUN = function(x) {
    x <- convert_human_to_mouse_symbols(x)
  }
)

s.genes <- geneList$s.genes # S phase
g2m.genes <- geneList$g2m.genes # G2/M phase

# Cell cycle scores
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
```

Apply ```SCTransform``` normalisation using ```glmGamPoi``` method:
```R
seurat.list <- lapply(
  X = seurat.list,
  FUN = function(x) {
    x <-
      SCTransform(x,
                  vars.to.regress = c('percent.mt', 'CC.Difference'),
                  method = 'glmGamPoi')
  }
)
```

Clean up unnecessary slots from the S4 Object:

```R
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
```

Integrate Seurat objects using 3000 features:
```R
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

# save integrated Seurat Object:
saveRDS(seurat, file = paste0(date, '_integrated_seurat_raw.rds'))
```
