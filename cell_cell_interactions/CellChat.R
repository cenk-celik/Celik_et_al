main <- paste0(system("echo ${HOME}/", intern = TRUE),"nmrc/")

condition_oi <- 'infected' # 'uninfected' for the healing dataset.
setwd(paste0(main, '/plots/CellChat/', condition_oi,'/'))

# devtools::install_github("sqjin/CellChat")
library(Seurat)
library(CellChat)
library(patchwork)
options(stringsAsFactors = F)

# Load Seurat Object----
seurat <- readRDS(paste0(main, 'seurat.rds'))
DefaultAssay(seurat) <- 'SCT'
seurat[['integrated']] <- NULL

seurat <- RenameIdents(
  seurat,
  c(
    'BAS1' = 'Keratinocyte',
    'BAS2' = 'Keratinocyte',
    'BAS3' = 'Keratinocyte',
    'SUP1' = 'Keratinocyte',
    'SUP2' = 'Keratinocyte',
    'SUP3' = 'Keratinocyte',
    'SUP4' = 'Keratinocyte',
    'TRM' = 'TRM',
    'M2' = 'M2-like',
    'EC1' = 'Endothelial',
    'EC2' = 'Endothelial',
    'HFSUP1' = 'Keratinocyte',
    'HFSUP2' = 'Keratinocyte',
    'OB' = 'Keratinocyte',
    'SG' = 'Keratinocyte',
    'NEUT1' = 'Neutrophil',
    'NEUT2' = 'Neutrophil',
    'TMEM' = 'Memory T cell',
    'DC.LC' = 'DC/LC',
    'FIB1' = 'Fibroblast',
    'FIB2' = 'Fibroblast',
    'FIB3' = 'Fibroblast',
    'NEU' = 'Neuron',
    'VSM' = 'VSM'
  )
)

seurat$seurat_clusters <- Idents(seurat)

split.list <- SplitObject(seurat, split.by = 'condition')

# Seurat to CellChat----
cellchat <- createCellChat(object = split.list$uninfected,
                           group.by = 'seurat_clusters',
                           assay = 'SCT')

# CellChat mouse DB----
CellChatDB <- subsetDB(CellChatDB.mouse, search = 'Secreted Signaling')
cellchat@DB <- CellChatDB

# Preprocess exprs data----
cellchat <- subsetData(cellchat)
future::plan('multisession', workers = 4)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- projectData(cellchat, PPI.mouse)

# Infer cell-cell communication network----
cellchat <- computeCommunProb(cellchat, population.size = T)
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Signalling pathway----
cellchat <- computeCommunProbPathway(cellchat)

# Calculate aggregated cell-cell network----
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))

pdf('counts_weights.pdf', width = 10, height = 10)
par(mfrow = c(1,2), xpd = T)
netVisual_circle(cellchat@net$count,
                 vertex.weight = groupSize,
                 weight.scale = T,
                 label.edge = F,
                 title.name = 'Number of interactions')
netVisual_circle(cellchat@net$weight,
                 vertex.weight = groupSize,
                 weight.scale = T,
                 label.edge = F,
                 title.name = 'Interaction weights/strength')
dev.off()

mat <- cellchat@net$weight
pdf('celltype_weights.pdf', width = 12, height = 12)
par(mfrow = c(2,5), xpd = T)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

# Visualise cell-cell network----
pathways.list <- cellchat@netP$pathways
vertex.receiver <- seq(1,4)

for (i in 1:length(pathways.list)) {
  par(mfrow = c(1,1))
  pdf(paste0('cell_cell_network/', pathways.list[i],'.pdf'), width = 5, height = 5)
  netVisual_aggregate(cellchat,
                      signaling = pathways.list[i],
                      vertex.receiver = vertex.receiver,
                      layout = 'circle')
  dev.off()
}

# Compute contribution of LR pairs----
netAnalysis_contribution(cellchat, signaling = pathways.list)

# Loop to save all inferred network for exploration----
pathways.show.all <- cellchat@netP$pathways
levels(cellchat@idents)

vertex.receiver <- seq(1,4)
for (i in 1:length(pathways.show.all)) {
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = 'hierarchy')
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename = paste0('./all_pathways/', pathways.show.all[i], '_LR_contribution.pdf'),
         plot = gg, device = 'pdf', width = 3, height = 5, units = 'in', dpi = 300)
}

# Visualise communication mediated by multiple LR or pathways----
library(ggplot2)
sources.use <- c('Keratinocyte','M2-like','TRM','Fibroblast','Neutrophil')
targets.use <- c('Endothelial')
signalling <- head(pathways.list)

i = 6
pdf(paste0('multiple_LR/', signalling[i], '_bubble_LR.pdf'), width = 5, height = 5)
netVisual_bubble(cellchat, sources.use = sources.use, targets.use = targets.use,
                remove.isolate = F, signaling = signalling[i]) +
  scale_color_gradient(high = "indianred2", low = "whitesmoke")
dev.off()

pairLR.use <- extractEnrichedLR(cellchat, signaling = signalling)
pdf(paste0('bubble_enrichedLR.pdf'), width = 3.5, height = 5.5)
netVisual_bubble(cellchat,
                 sources.use = sources.use,
                 targets.use = targets.use,
                 remove.isolate = T,
                 pairLR.use = pairLR.use)+
  scale_color_gradient(high = "indianred2", low = "whitesmoke")
dev.off()

# Plot signalling gene exprs----
i = 6
par(mfrow = c(1,1), xpd = T)
pdf(paste0('expressions/', signalling[i], '_violin_exprs.pdf'), width = 5, height = 5)
plotGeneExpression(cellchat, signaling = signalling[i])
dev.off()

# System analysis of cell-cell network----
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = 'netP')

# Visualise dominant senders and receivers in 2D
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = signalling)
gg1 + gg2

setwd(main)
saveRDS(cellchat, file = paste0(condition_oi,'_cellchat.rds'))
cellchat <- readRDS(paste0(main, condition_oi,'_cellchat.rds'))

sessionInfo()
