main <- paste0(system("echo ${HOME}/", intern = TRUE),"nmrc/")

library(CellChat)
library(patchwork)
options(stringsAsFactors = F)

# Load cellchat objects----
cellchat.Inf <- readRDS('infected_cellchat.rds')
cellchat.Uninf <- readRDS('uninfected_cellchat.rds')

# Merge objects----
object.list <- list(Uninfected = cellchat.Uninf, Infected = cellchat.Inf)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

setwd(paste0(main, 'plots/CellChat/comparison/'))

colours <- c('Keratinocyte' = '#d0a2c9','TRM' = '#00c23c','M2-like' = '#bb00ff',
             'Endothelial' = '#954321','Neutrophil' = '#f4cb9f','Memory T cell' = '#005d31',
             'DC/LC' = '#3473b6', 'Fibroblast' = '#4b1c5b', 'Neuron' = '#52b549',
             'VSM' = '#bd1b86')

# Compare number of interactions-strengths----
int1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), color.use = c('#6E385D', '#5c761f'))
int2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = 'weight', color.use = c('#6E385D', '#5c761f'))

pdf(paste0('interaction_comparison.pdf'), width = 3, height = 1.5)
int1 + int2
dev.off()

# Differential number of interactions----
pdf('Differential_number_of_interactions_strength.pdf', width = 10, height = 10)
netVisual_diffInteraction(cellchat, weight.scale = T, color.use = colours)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = 'weight', color.use = colours)
dev.off()

weight.max <- getMaxWeight(object.list, attribute = c('idents', 'count'))
pdf('Number_of_interactions.pdf', width = 10, height = 10)
par(mfrow = c(1,2), xpd = T)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge = T, 
                   edge.weight.max = weight.max[2], color.use = colours,
                   edge.width.max = 12, title.name = paste0('Number of interactions - ', names(object.list)[i]))
}
dev.off()

# Differential number of interactions----
group.cellType <-
  c(
    'Keratinocyte',
    'TRM',
    'M2-like',
    'Endothelial',
    'Neutrophil',
    'Fibroblast'
  )

object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

weight.max <- getMaxWeight(object.list, slot.name = c('idents','net','net'), attribute = c('idents','count','count.merged'))

pdf('Differential_number_of_interactions_strength_selected.pdf', width = 10, height = 10)
par(mfrow = c(1,2), xpd = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = 'count.merged', label.edge = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = 'weight.merged', label.edge = T)
dev.off()

# Compare major sources in 2D----
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link))
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}

patchwork::wrap_plots(plots = gg)

gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = 'Fibroblast', color.use = c('darkgrey', '#6E385D', '#5c761f')) + theme(aspect.ratio = 1)
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = 'M2-like', color.use = c('darkgrey', '#6E385D', '#5c761f')) + theme(aspect.ratio = 1)
gg3 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = 'Keratinocyte', color.use = c('darkgrey', '#6E385D', '#5c761f')) + theme(aspect.ratio = 1)
gg4 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = 'Neutrophil', color.use = c('darkgrey', '#6E385D', '#5c761f')) + theme(aspect.ratio = 1)

pdf('Signalling_changes.pdf', width = 10, height = 10)
patchwork::wrap_plots(plots = list(gg1, gg2, gg3, gg4))
dev.off()

# Compare overall information flow----
gg1 <- rankNet(cellchat, mode = 'comparison', stacked = T, do.stat = T , color.use = c('darkgrey', '#6E385D', '#5c761f'))
gg2 <- rankNet(cellchat, mode = 'comparison', stacked = F, do.stat = T, color.use = c('darkgrey', '#6E385D', '#5c761f'))
pdf('Information_flow.pdf', width = 4, height = 6)
rankNet(cellchat, mode = 'comparison', stacked = T, do.stat = T, color.use = c('darkgrey', '#6E385D', '#5c761f'))
dev.off()

# Compare outgoing signalling----
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
signalling <- c('SPP1','TNF','CXCL','EGF','CSF3','ANNEXIN','IL1','CHEMERIN', # INFECTED
                'IGF','WNT','IL10','ncWNT','TWEAK','CD137','IL6','GAS','EDA','TGFb','PDGF','VISFATIN','VEGF','MIF' # UNINFECTED
                )
# Identify LR pairs----
sources.use <- c('Fibroblast','Keratinocyte','M2-like')
targets.use <- c('Endothelial')
library(scales)

pdf('./signalling_bubble_pairs.pdf', width = 4, height = 6)
netVisual_bubble(
  cellchat,
  sources.use = sources.use,
  targets.use = targets.use,
  comparison = c(1, 2),
  angle.x = 45, n.colors = 3,
  signaling = signalling,
  color.text = c('#6E385D', '#5c761f')
) + 
  scale_color_gradient2(high = "darkred", low = "whitesmoke", midpoint = 0.09)
dev.off()

# LR up&down in infected
gg1 <-
  netVisual_bubble(
    cellchat,
    signaling = signalling,
    sources.use = sources.use,
    targets.use = targets.use,
    comparison = c(1, 2),
    max.dataset = 2,
    title.name = 'Increased signalling in Infected',
    angle.x = 45,
    remove.isolate = T,
    color.text = c('#6E385D', '#5c761f')
  ) + 
  scale_color_gradient2(high = "darkred", low = "whitesmoke", midpoint = 0.09)

gg2 <-
  netVisual_bubble(
    cellchat,
    signaling = signalling,
    sources.use = sources.use,
    targets.use = targets.use,
    comparison = c(1, 2),
    max.dataset = 1,
    title.name = 'Decreased signalling in Infected',
    angle.x = 45,
    remove.isolate = T,
    color.text = c('#6E385D', '#5c761f')
  ) + 
  scale_color_gradient2(high = "darkred", low = "whitesmoke", midpoint = 0.09)

pdf(paste0('./signalling_bubble.pdf'), width = 8, height = 5)
gg1 + gg2
dev.off()

# Identify dysfunctional signalling by DEGs----
pos.dataset <- 'Infected'
features.name <- pos.dataset

cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = 'datasets',
                                       pos.dataset = pos.dataset,
                                       features.name = features.name,
                                       only.pos = F, thresh.pc = 0.1,
                                       thresh.fc = 0.1, thresh.p = 1)

net <- netMappingDEG(cellchat, features.name = features.name)

net.up <- subsetCommunication(cellchat, net = net, datasets = 'Infected',
                              ligand.logFC = 0.2, receptor.logFC = NULL)

net.down <- subsetCommunication(cellchat, net = net, datasets = 'Uninfected',
                                ligand.logFC = -0.1, receptor.logFC = -0.1)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

pairLR.use.up <- net.up[ , 'interaction_name', drop = F]
gg1 <- netVisual_bubble(cellchat, direction = 1, #color.heatmap = 'viridis',
                        pairLR.use = pairLR.use.up, sources.use = sources.use,
                        targets.use = targets.use, comparison = c(1, 2),
                        angle.x = 45, remove.isolate = T,
                        title.name = paste0('Upregulated signalling in ', names(object.list)[2]), 
                        color.text = c('#6E385D', '#5c761f')) + 
  scale_color_gradient2(high = "indianred2", low = "white")

pairLR.use.down <- net.down[ , 'interaction_name', drop = F]
gg2 <- netVisual_bubble(cellchat, direction = 1, #color.heatmap = 'viridis',
                        pairLR.use = pairLR.use.down, sources.use = sources.use,
                        targets.use = targets.use, comparison = c(1, 2),
                        angle.x = 45, remove.isolate = T,
                        title.name = paste0('Downregulated signalling in ', names(object.list)[2]), 
                        color.text = c('#6E385D', '#5c761f')) + 
  scale_color_gradient2(high = "indianred2", low = "white")

pdf('signalling_comparison_all.pdf', width = 12, height = 8)
gg1 + gg2
dev.off()

# Chord
par(mfrow = c(1,2), xpd = T)
pdf('./chord_gene_infected_up_down.pdf', width = 15, height = 15)
netVisual_chord_gene(object.list[[2]], sources.use = sources.use,
                     targets.use = targets.use, slot.name = 'net',
                     net = net.up, lab.cex = 0.5, small.gap = 3.5,
                     title.name = paste0('Upregulated signalling in ', names(object.list)[2]))

netVisual_chord_gene(object.list[[2]], sources.use = sources.use,
                     targets.use = targets.use, slot.name = 'net',
                     net = net.down, lab.cex = 0.5, small.gap = 3.5,
                     title.name = paste0('Downregulated signalling in ', names(object.list)[2]))
dev.off()

# Compare communication using hierarchy plot----
infected.signalling <- c('SPP1','TNF','CXCL','EGF','CSF3','ANNEXIN')
uninfected.signalling <- c('IL6', 'GAS', 'EDA', 'TGFb', 'PDGF', 'VEGF', 'MIF')

lapply(infected.signalling, FUN = function(x) {
  weight.max <- getMaxWeight(object.list, slot.name = 'netP', attribute = x)
  pdf(paste0('networks/', x, '_network_infected.pdf'), width = 10, height = 10)
  par(mfrow = c(1,2), xpd = T)
  for (i in 2:(length(object.list)) ) {
    netVisual_aggregate(object.list[[i]], signaling = x, layout = 'circle',
                        edge.weight.max = weight.max[1], edge.width.max = 10,
                        signaling.name = paste(x, names(object.list)[i]),
                        color.use = colours)
    
    netVisual_aggregate(object.list[[i]], signaling = x, layout = 'chord',
                        signaling.name = paste(x, names(object.list)[i]), 
                        color.use = colours)
  }
  dev.off()
})

# Compare gene expression----
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c('Uninfected','Infected'))

i = 1 # to 7
pdf(paste0('violin/', uninfected.signalling[i], '_violin.pdf'), width = 10, height = 5)
plotGeneExpression(cellchat,
                   signaling = uninfected.signalling[i],
                   split.by = 'datasets',
                   colors.ggplot = T,
                   enriched.only = T)
dev.off()

# Compare cell-cell communictation----
lapply(infected.signalling, FUN = function(x){
  pdf(paste0('aggregate/', x, '_aggregate.pdf'), width = 12, height = 10)
  weight.max <- getMaxWeight(object.list, slot.name = c('netP'), attribute = x)
  par(mfrow = c(1,2), xpd = T)
  for (i in 1:length(object.list)) {
    netVisual_aggregate(object.list[[i]],
                        signaling = x,
                        layout = 'chord',
                        edge.weight.max = weight.max[1],
                        edge.width.max = 10,
                        signaling.name = paste(x, names(object.list)[i]),
                        color.use = colours)
  }
  dev.off()
})

# Signalling from----
sources.use <-  c('TRM', 'M2-like')

pdf(paste0('signallingFrom/', sources.use, '_interactions.pdf'), width = 20, height = 8)
par(mfrow = c(1,2), xpd = T)
#layout(matrix(1:2, 1, 2))
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = sources.use, targets.use = c(1:10),
                       signaling = signalling, color.use = colours, scale = T,
                       title.name = paste0('Signalling from ', sources.use, '-', names(object.list)[i]))
}
dev.off()

sessionInfo()
