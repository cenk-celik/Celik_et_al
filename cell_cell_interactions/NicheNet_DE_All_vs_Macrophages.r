library(nichenetr)
library(RColorBrewer)
library(tidyverse)
library(Seurat)

rm(list = ls(all.names = T)); gc(full = T)

##----Prepare SeuratObject----
seurat_obj <- readRDS('seurat.rds')
DefaultAssay(seurat_obj) <- 'RNA'
seurat_obj[['SCT']] <- NULL
seurat_obj[['integrated']] <- NULL

seurat_obj <- RenameIdents(
  seurat_obj,
  c(
    'BAS1' = 'Keratinocyte',
    'BAS2' = 'Keratinocyte',
    'BAS3' = 'Keratinocyte',
    'SUP1' = 'Keratinocyte',
    'SUP2' = 'Keratinocyte',
    'SUP3' = 'Keratinocyte',
    'SUP4' = 'Keratinocyte',
    'M2' = 'M2-like',
    'TRM' = 'Tissue resident macrophage',
    'EC1' = 'Endothelial',
    'EC2' = 'Endothelial',
    'NEUT1' = 'Neutrophil',
    'NEUT2' = 'Neutrophil',
    'TMEM' = 'Memory T cell',
    'DC' = 'Dendritic cell',
    'FIB1' = 'Fibroblast',
    'FIB2' = 'Fibroblast',
    'FIB3' = 'Fibroblast',
    'NEU' = 'Neuron',
    'VSM' = 'Vascular smooth muscle'
  )
)

seurat_obj <-
  subset(
    seurat_obj,
    idents = c(
      'Keratinocyte',
      'Fibroblast',
      'Tissue resident macrophage',
      'M2-like',
      'Neutrophil',
      'Memory T cell',
      'Dendritic cell',
      'Neuron',
      'Vascular smooth muscle',
      'Endothelial'
    )
  )

seurat_obj <- NormalizeData(seurat_obj)

seurat_obj$seurat_clusters <- Idents(seurat_obj)

seurat_obj@meta.data$celltype_condition = paste(seurat_obj@meta.data$seurat_clusters,
                                                seurat_obj@meta.data$condition,
                                                sep = '_')

celltype_id = 'celltype_condition'
seurat_obj = SetIdent(seurat_obj, value = seurat_obj[[celltype_id]])

ligand_target_matrix = readRDS('nichenetr/ligand_target_matrix.rds')
ligand_target_matrix[1:5, 1:5]

lr_network = readRDS('nichenetr/lr_network.rds')
lr_network = lr_network %>% mutate(bonafide = !database %in% c('ppi_prediction', 'ppi_prediction_go'))
lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor, bonafide)

head(lr_network)

organism = 'mouse'

if (organism == 'mouse') {
  lr_network = lr_network %>% mutate(ligand = convert_human_to_mouse_symbols(ligand), receptor = convert_human_to_mouse_symbols(receptor)) %>% drop_na()
  
  colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
  rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols()
  ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]
}

##----Tissue resident macrophages----
# Define the niches/microenvironments of interest
niches = list(
  'infected_niche' = list(
    'sender' = c(
      'Keratinocyte_infected',
      'Fibroblast_infected',
      'Endothelial_infected',
      'Dendritic cell_infected',
      'Neutrophil_infected',
      'Neuron_infected',
      'Vascular smooth muscle_infected',
      'Tissue resident macrophage_infected',
      'M2-like_infected',
      'Memory T cell_infected'
    ),
    'receiver' = c('Tissue resident macrophage_infected')
  ),
  'uninfected_niche' = list(
    'sender' = c(
      'Keratinocyte_uninfected',
      'Fibroblast_uninfected',
      'Endothelial_uninfected',
      'Dendritic cell_uninfected',
      'Neutrophil_uninfected',
      'Neuron_uninfected',
      'Vascular smooth muscle_uninfected',
      'Tissue resident macrophage_uninfected',
      'M2-like_uninfected',
      'Memory T cell_uninfected'
    ),
    'receiver' = c('Tissue resident macrophage_uninfected')
  )
)

# Calculate differential expression between the niches
assay_oi = 'RNA'

DE_sender = calculate_niche_de(
  seurat_obj = seurat_obj %>% subset(features = lr_network$ligand %>% unique()),
  niches = niches,
  type = 'sender',
  assay_oi = assay_oi
)

DE_receiver = calculate_niche_de(
  seurat_obj = seurat_obj %>% subset(features = lr_network$receptor %>% unique()),
  niches = niches,
  type = 'receiver',
  assay_oi = assay_oi
)

DE_sender = DE_sender %>% mutate(avg_log2FC = ifelse(
  avg_log2FC == Inf,
  max(avg_log2FC[is.finite(avg_log2FC)]),
  ifelse(avg_log2FC == -Inf, min(avg_log2FC[is.finite(avg_log2FC)]), avg_log2FC)
))

DE_receiver = DE_receiver %>% mutate(avg_log2FC = ifelse(
  avg_log2FC == Inf,
  max(avg_log2FC[is.finite(avg_log2FC)]),
  ifelse(avg_log2FC == -Inf, min(avg_log2FC[is.finite(avg_log2FC)]), avg_log2FC)
))

# process DE results
expression_pct = 0.10

DE_sender_processed = process_niche_de(
  DE_table = DE_sender,
  niches = niches,
  expression_pct = expression_pct,
  type = 'sender'
)

DE_receiver_processed = process_niche_de(
  DE_table = DE_receiver,
  niches = niches,
  expression_pct = expression_pct,
  type = 'receiver'
)

# combine sender-receiver DR based on L-R pairs
specificity_score_LR_pairs = 'min_lfc'

DE_sender_receiver = combine_sender_receiver_de(DE_sender_processed,
                                                DE_receiver_processed,
                                                lr_network,
                                                specificity_score = specificity_score_LR_pairs)

# Skip spatial info
include_spatial_info_sender = F
include_spatial_info_receiver = F

spatial_info = tibble(
  celltype_region_oi = '',
  celltype_other_region = '',
  niche =  'infected_niche',
  celltype_type = 'sender'
)

specificity_score_spatial = 'lfc'

if (include_spatial_info_sender == F &
    include_spatial_info_receiver == F) {
  spatial_info = tibble(celltype_region_oi = NA,
                        celltype_other_region = NA) %>% mutate(niche =  niches %>% names() %>% head(1),
                                                               celltype_type = 'sender')
}

if (include_spatial_info_sender == T) {
  sender_spatial_DE = calculate_spatial_DE(
    seurat_obj = seurat_obj %>% subset(features = lr_network$ligand %>% unique()),
    spatial_info = spatial_info %>% filter(celltype_type == 'sender')
  )
  sender_spatial_DE_processed = process_spatial_de(
    DE_table = sender_spatial_DE,
    type = 'sender',
    lr_network = lr_network,
    expression_pct = expression_pct,
    specificity_score = specificity_score_spatial
  )
  
  sender_spatial_DE_others = get_non_spatial_de(
    niches = niches,
    spatial_info = spatial_info,
    type = 'sender',
    lr_network = lr_network
  )
  sender_spatial_DE_processed = sender_spatial_DE_processed %>% bind_rows(sender_spatial_DE_others)
  
  sender_spatial_DE_processed = sender_spatial_DE_processed %>% mutate(scaled_ligand_score_spatial = scale_quantile_adapted(ligand_score_spatial))
  
} else {
  sender_spatial_DE_processed = get_non_spatial_de(
    niches = niches,
    spatial_info = spatial_info,
    type = 'sender',
    lr_network = lr_network
  )
  sender_spatial_DE_processed = sender_spatial_DE_processed %>% mutate(scaled_ligand_score_spatial = scale_quantile_adapted(ligand_score_spatial))
  
}

if (include_spatial_info_receiver == T) {
  receiver_spatial_DE = calculate_spatial_DE(
    seurat_obj = seurat_obj %>% subset(features = lr_network$receptor %>% unique()),
    spatial_info = spatial_info %>% filter(celltype_type == 'receiver')
  )
  receiver_spatial_DE_processed = process_spatial_de(
    DE_table = receiver_spatial_DE,
    type = 'receiver',
    lr_network = lr_network,
    expression_pct = expression_pct,
    specificity_score = specificity_score_spatial
  )
  
  receiver_spatial_DE_others = get_non_spatial_de(
    niches = niches,
    spatial_info = spatial_info,
    type = 'receiver',
    lr_network = lr_network
  )
  receiver_spatial_DE_processed = receiver_spatial_DE_processed %>% bind_rows(receiver_spatial_DE_others)
  
  receiver_spatial_DE_processed = receiver_spatial_DE_processed %>% mutate(scaled_receptor_score_spatial = scale_quantile_adapted(receptor_score_spatial))
  
} else {
  receiver_spatial_DE_processed = get_non_spatial_de(
    niches = niches,
    spatial_info = spatial_info,
    type = 'receiver',
    lr_network = lr_network
  )
  receiver_spatial_DE_processed = receiver_spatial_DE_processed %>% mutate(scaled_receptor_score_spatial = scale_quantile_adapted(receptor_score_spatial))
}

# Calculate ligand activities and infer active ligand-target links
lfc_cutoff = 0.5
specificity_score_targets = 'min_lfc'

DE_receiver_targets = calculate_niche_de_targets(
  seurat_obj = seurat_obj,
  niches = niches,
  lfc_cutoff = lfc_cutoff,
  expression_pct = expression_pct,
  assay_oi = assay_oi
)

DE_receiver_processed_targets = process_receiver_target_de(
  DE_receiver_targets = DE_receiver_targets,
  niches = niches,
  expression_pct = expression_pct,
  specificity_score = specificity_score_targets
)

background = DE_receiver_processed_targets  %>% pull(target) %>% unique()

geneset_niche1 = DE_receiver_processed_targets %>% filter(
  receiver == niches[[1]]$receiver &
    target_score >= lfc_cutoff &
    target_significant == 1 &
    target_present == 1
) %>% pull(target) %>% unique()

geneset_niche2 = DE_receiver_processed_targets %>% filter(
  receiver == niches[[2]]$receiver &
    target_score >= lfc_cutoff &
    target_significant == 1 &
    target_present == 1
) %>% pull(target) %>% unique()

geneset_niche1 %>% setdiff(rownames(ligand_target_matrix))
geneset_niche2 %>% setdiff(rownames(ligand_target_matrix))

length(geneset_niche1)
length(geneset_niche2)

top_n_target = 250

niche_geneset_list = list(
  'infected_niche' = list(
    'receiver' = niches[[1]]$receiver,
    'geneset' = geneset_niche1,
    'background' = background
  ),
  'uninfected_niche' = list(
    'receiver' = niches[[2]]$receiver,
    'geneset' = geneset_niche2 ,
    'background' = background
  )
)

ligand_activities_targets = get_ligand_activities_targets(
  niche_geneset_list = niche_geneset_list,
  ligand_target_matrix = ligand_target_matrix,
  top_n_target = top_n_target
)

# Calculate scaled expression of ligands, receptors and targets across cell types of interest
features_oi = union(lr_network$ligand, lr_network$receptor) %>% 
  union(ligand_activities_targets$target) %>% 
  setdiff(NA)

dotplot = suppressWarnings(Seurat::DotPlot(
  seurat_obj %>% 
    subset(idents = niches %>% unlist() %>% unique()),
  features = features_oi,
  assay = assay_oi
))

exprs_tbl = dotplot$data %>% as_tibble()

exprs_tbl = exprs_tbl %>% rename(
  celltype = id,
  gene = features.plot,
  expression = avg.exp,
  expression_scaled = avg.exp.scaled,
  fraction = pct.exp) %>%
  mutate(fraction = fraction / 100) %>%
  as_tibble() %>%
  select(celltype, gene, expression, expression_scaled, fraction) %>%
  distinct() %>%
  arrange(gene) %>% 
  mutate(gene = as.character(gene))

exprs_tbl_ligand = exprs_tbl %>%
  filter(gene %in% lr_network$ligand) %>%
  rename(
    sender = celltype,
    ligand = gene,
    ligand_expression = expression,
    ligand_expression_scaled = expression_scaled,
    ligand_fraction = fraction
  )

exprs_tbl_receptor = exprs_tbl %>%
  filter(gene %in% lr_network$receptor) %>%
  rename(
    receiver = celltype,
    receptor = gene,
    receptor_expression = expression,
    receptor_expression_scaled = expression_scaled,
    receptor_fraction = fraction
  )

exprs_tbl_target = exprs_tbl %>% 
  filter(gene %in% ligand_activities_targets$target) %>%
  rename(
    receiver = celltype,
    target = gene,
    target_expression = expression,
    target_expression_scaled = expression_scaled,
    target_fraction = fraction
  )

exprs_tbl_ligand = exprs_tbl_ligand %>% 
  mutate(scaled_ligand_expression_scaled = scale_quantile_adapted(ligand_expression_scaled)) %>% 
  mutate(ligand_fraction_adapted = ligand_fraction) %>% 
  mutate_cond(ligand_fraction >= expression_pct, ligand_fraction_adapted = expression_pct) %>% 
  mutate(scaled_ligand_fraction_adapted = scale_quantile_adapted(ligand_fraction_adapted))

exprs_tbl_receptor = exprs_tbl_receptor %>% 
  mutate(scaled_receptor_expression_scaled = scale_quantile_adapted(receptor_expression_scaled)) %>% 
  mutate(receptor_fraction_adapted = receptor_fraction) %>% 
  mutate_cond(receptor_fraction >= expression_pct, receptor_fraction_adapted = expression_pct) %>% 
  mutate(scaled_receptor_fraction_adapted = scale_quantile_adapted(receptor_fraction_adapted))

# Expression fraction and receptor
exprs_sender_receiver = lr_network %>%
  inner_join(exprs_tbl_ligand, by = c('ligand')) %>%
  inner_join(exprs_tbl_receptor, by = c('receptor')) %>% 
  inner_join(DE_sender_receiver %>% distinct(niche, sender, receiver))

ligand_scaled_receptor_expression_fraction_df = exprs_sender_receiver %>%
  group_by(ligand, receiver) %>%
  mutate(
    rank_receptor_expression = dense_rank(receptor_expression),
    rank_receptor_fraction  = dense_rank(receptor_fraction)
  ) %>%
  mutate(ligand_scaled_receptor_expression_fraction = 0.5 * ((
    rank_receptor_fraction / max(rank_receptor_fraction)
  ) +
    ((
      rank_receptor_expression / max(rank_receptor_expression)
    )))) %>%
  distinct(ligand,
           receptor,
           receiver,
           ligand_scaled_receptor_expression_fraction,
           bonafide) %>% distinct() %>% ungroup()

# Prioritization of ligand-receptor and ligand-target links
prioritizing_weights = c(
  'scaled_ligand_score' = 5,
  'scaled_ligand_expression_scaled' = 1,
  'ligand_fraction' = 1,
  'scaled_ligand_score_spatial' = 0,
  'scaled_receptor_score' = 0.5,
  'scaled_receptor_expression_scaled' = 0.5,
  'receptor_fraction' = 1,
  'ligand_scaled_receptor_expression_fraction' = 1,
  'scaled_receptor_score_spatial' = 0,
  'scaled_activity' = 0,
  'scaled_activity_normalized' = 1,
  'bona_fide' = 1
)

output = list(
  DE_sender_receiver = DE_sender_receiver,
  ligand_scaled_receptor_expression_fraction_df = ligand_scaled_receptor_expression_fraction_df,
  sender_spatial_DE_processed = sender_spatial_DE_processed,
  receiver_spatial_DE_processed = receiver_spatial_DE_processed,
  ligand_activities_targets = ligand_activities_targets,
  DE_receiver_processed_targets = DE_receiver_processed_targets,
  exprs_tbl_ligand = exprs_tbl_ligand,
  exprs_tbl_receptor = exprs_tbl_receptor,
  exprs_tbl_target = exprs_tbl_target
)

prioritization_tables = get_prioritization_tables(output, prioritizing_weights)

prioritization_tables$prioritization_tbl_ligand_receptor %>% 
  filter(receiver == niches[[1]]$receiver) %>% head(10)
prioritization_tables$prioritization_tbl_ligand_target %>% 
  filter(receiver == niches[[1]]$receiver) %>% head(10)
prioritization_tables$prioritization_tbl_ligand_receptor %>% 
  filter(receiver == niches[[2]]$receiver) %>% head(10)
prioritization_tables$prioritization_tbl_ligand_target %>% 
  filter(receiver == niches[[2]]$receiver) %>% head(10)

# Visualization of the Differential NicheNet output
top_ligand_niche_df = prioritization_tables$prioritization_tbl_ligand_receptor %>%
  select(niche, sender, receiver, ligand, receptor, prioritization_score) %>%
  group_by(ligand) %>%
  top_n(1, prioritization_score) %>%
  ungroup() %>%
  select(ligand, receptor, niche) %>%
  rename(top_niche = niche)

top_ligand_receptor_niche_df = prioritization_tables$prioritization_tbl_ligand_receptor %>%
  select(niche, sender, receiver, ligand, receptor, prioritization_score) %>%
  group_by(ligand, receptor) %>%
  top_n(1, prioritization_score) %>%
  ungroup() %>%
  select(ligand, receptor, niche) %>%
  rename(top_niche = niche)

ligand_prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>%
  select(niche, sender, receiver, ligand, prioritization_score) %>%
  group_by(ligand, niche) %>%
  top_n(2, prioritization_score) %>%
  ungroup() %>%
  distinct() %>%
  inner_join(top_ligand_niche_df) %>%
  filter(niche == top_niche) %>%
  group_by(niche) %>%
  top_n(20, prioritization_score) %>%
  ungroup()

receiver_oi = 'Tissue resident macrophage_infected'

filtered_ligands = ligand_prioritized_tbl_oi %>%
  filter(receiver == receiver_oi) %>%
  pull(ligand) %>%
  unique()

prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>%
  filter(ligand %in% filtered_ligands) %>%
  select(niche,
         sender,
         receiver,
         ligand,
         receptor,
         ligand_receptor,
         prioritization_score) %>%
  distinct() %>%
  inner_join(top_ligand_receptor_niche_df) %>%
  group_by(ligand) %>%
  filter(receiver == receiver_oi) %>%
  top_n(2, prioritization_score) %>%
  ungroup()

lfc_plot = make_ligand_receptor_lfc_plot(
  receiver_oi,
  prioritized_tbl_oi,
  prioritization_tables$prioritization_tbl_ligand_receptor,
  plot_legend = F,
  heights = NULL,
  widths = NULL
) 

lfc_plot

pdf(paste0('./plots/nichenet/nichenet_DE_lfc_plot_all_', receiver_oi,'.pdf'),
    width = 8,
    height = 9)
plot(lfc_plot)
dev.off()

exprs_activity_target_plot = make_ligand_activity_target_exprs_plot(
  receiver_oi,
  prioritized_tbl_oi,
  prioritization_tables$prioritization_tbl_ligand_receptor,
  prioritization_tables$prioritization_tbl_ligand_target,
  output$exprs_tbl_ligand,
  output$exprs_tbl_target,
  lfc_cutoff,
  ligand_target_matrix,
  plot_legend = F,
  heights = NULL,
  widths = NULL
)

exprs_activity_target_plot$combined_plot

pdf(paste0(
  './plots/nichenet/nichenet_DE_exprs_activity_target_plot_all_',receiver_oi,'.pdf'),
  width = 18,
  height = 10
); plot(exprs_activity_target_plot$combined_plot); dev.off()

filtered_ligands = ligand_prioritized_tbl_oi %>%
  filter(receiver == receiver_oi) %>%
  top_n(20, prioritization_score) %>%
  pull(ligand) %>%
  unique()

prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>%
  filter(ligand %in% filtered_ligands) %>%
  select(niche,
         sender,
         receiver,
         ligand,
         receptor,
         ligand_receptor,
         prioritization_score) %>%
  distinct() %>%
  inner_join(top_ligand_receptor_niche_df) %>%
  group_by(ligand) %>%
  filter(receiver == receiver_oi) %>%
  top_n(3, prioritization_score) %>%
  ungroup()

colors_sender <- c(
  'Keratinocyte_infected' = '#8D7E30',
  'Dendritic cell_infected' = '#3473B8',
  'Endothelial_infected' = '#974221',
  'Fibroblast_infected' = '#4B1C5B',
  'Memory T cell_infected' = '#005D31',
  'Neuron_infected' = '#51B749',
  'Neutrophil_infected' = '#FCCB99',
  'Vascular smooth muscle_infected' = '#C01988',
  'Tissue resident macrophage_infected' = '#59BA5D',
  'M2-like_infected' = '#9b32f5'
)

colors_receiver = c('#59BA5D') %>% #M2 A12BE0
  magrittr::set_names(prioritized_tbl_oi$receiver %>%
                        unique() %>%
                        sort())

transparency <- prioritized_tbl_oi %>% 
  mutate(weight = (prioritization_score - min(prioritization_score)) / (max(prioritization_score) - min(prioritization_score))) %>% 
  mutate(transparency = 1 - prioritization_score) %>% .$transparency

pdf(paste0(
  './plots/nichenet/nichenet_DE_exprs_activity_circos_all_',receiver_oi,'.pdf'),
  width = 10,
  height = 10
) 
circos_output = make_circos_lr(prioritized_tbl_oi, 
                               colors_sender, 
                               colors_receiver,
                               scale = T, border = F,
                               transparency = transparency)
dev.off()

##----M2-like----
# Define the niches/microenvironments of interest
niches = list(
  'infected_niche' = list(
    'sender' = c(
      'Keratinocyte_infected',
      'Fibroblast_infected',
      'Endothelial_infected',
      'Dendritic cell_infected',
      'Neutrophil_infected',
      'Neuron_infected',
      'Vascular smooth muscle_infected',
      'Tissue resident macrophage_infected',
      'M2-like_infected',
      'Memory T cell_infected'
    ),
    'receiver' = c('M2-like_infected')
  ),
  'uninfected_niche' = list(
    'sender' = c(
      'Keratinocyte_uninfected',
      'Fibroblast_uninfected',
      'Endothelial_uninfected',
      'Dendritic cell_uninfected',
      'Neutrophil_uninfected',
      'Neuron_uninfected',
      'Vascular smooth muscle_uninfected',
      'Tissue resident macrophage_uninfected',
      'M2-like_uninfected',
      'Memory T cell_uninfected'
    ),
    'receiver' = c('M2-like_uninfected')
  )
)

# Calculate differential expression between the niches
assay_oi = 'RNA'

DE_sender = calculate_niche_de(
  seurat_obj = seurat_obj %>% subset(features = lr_network$ligand %>% unique()),
  niches = niches,
  type = 'sender',
  assay_oi = assay_oi
)

DE_receiver = calculate_niche_de(
  seurat_obj = seurat_obj %>% subset(features = lr_network$receptor %>% unique()),
  niches = niches,
  type = 'receiver',
  assay_oi = assay_oi
)

DE_sender = DE_sender %>% mutate(avg_log2FC = ifelse(
  avg_log2FC == Inf,
  max(avg_log2FC[is.finite(avg_log2FC)]),
  ifelse(avg_log2FC == -Inf, min(avg_log2FC[is.finite(avg_log2FC)]), avg_log2FC)
))

DE_receiver = DE_receiver %>% mutate(avg_log2FC = ifelse(
  avg_log2FC == Inf,
  max(avg_log2FC[is.finite(avg_log2FC)]),
  ifelse(avg_log2FC == -Inf, min(avg_log2FC[is.finite(avg_log2FC)]), avg_log2FC)
))

# process DE results
expression_pct = 0.10

DE_sender_processed = process_niche_de(
  DE_table = DE_sender,
  niches = niches,
  expression_pct = expression_pct,
  type = 'sender'
)

DE_receiver_processed = process_niche_de(
  DE_table = DE_receiver,
  niches = niches,
  expression_pct = expression_pct,
  type = 'receiver'
)

# combine sender-receiver DR based on L-R pairs
specificity_score_LR_pairs = 'min_lfc'

DE_sender_receiver = combine_sender_receiver_de(DE_sender_processed,
                                                DE_receiver_processed,
                                                lr_network,
                                                specificity_score = specificity_score_LR_pairs)

# Skip spatial info
include_spatial_info_sender = F
include_spatial_info_receiver = F

spatial_info = tibble(
  celltype_region_oi = 'Basal_infected',
  celltype_other_region = 'Suprabasal_infected',
  niche =  'infected_niche',
  celltype_type = 'sender'
)

specificity_score_spatial = 'lfc'

if (include_spatial_info_sender == F &
    include_spatial_info_receiver == F) {
  spatial_info = tibble(celltype_region_oi = NA,
                        celltype_other_region = NA) %>% mutate(niche =  niches %>% names() %>% head(1),
                                                               celltype_type = 'sender')
}

if (include_spatial_info_sender == T) {
  sender_spatial_DE = calculate_spatial_DE(
    seurat_obj = seurat_obj %>% subset(features = lr_network$ligand %>% unique()),
    spatial_info = spatial_info %>% filter(celltype_type == 'sender')
  )
  sender_spatial_DE_processed = process_spatial_de(
    DE_table = sender_spatial_DE,
    type = 'sender',
    lr_network = lr_network,
    expression_pct = expression_pct,
    specificity_score = specificity_score_spatial
  )
  
  sender_spatial_DE_others = get_non_spatial_de(
    niches = niches,
    spatial_info = spatial_info,
    type = 'sender',
    lr_network = lr_network
  )
  sender_spatial_DE_processed = sender_spatial_DE_processed %>% bind_rows(sender_spatial_DE_others)
  
  sender_spatial_DE_processed = sender_spatial_DE_processed %>% mutate(scaled_ligand_score_spatial = scale_quantile_adapted(ligand_score_spatial))
  
} else {
  sender_spatial_DE_processed = get_non_spatial_de(
    niches = niches,
    spatial_info = spatial_info,
    type = 'sender',
    lr_network = lr_network
  )
  sender_spatial_DE_processed = sender_spatial_DE_processed %>% mutate(scaled_ligand_score_spatial = scale_quantile_adapted(ligand_score_spatial))
  
}

if (include_spatial_info_receiver == T) {
  receiver_spatial_DE = calculate_spatial_DE(
    seurat_obj = seurat_obj %>% subset(features = lr_network$receptor %>% unique()),
    spatial_info = spatial_info %>% filter(celltype_type == 'receiver')
  )
  receiver_spatial_DE_processed = process_spatial_de(
    DE_table = receiver_spatial_DE,
    type = 'receiver',
    lr_network = lr_network,
    expression_pct = expression_pct,
    specificity_score = specificity_score_spatial
  )
  
  receiver_spatial_DE_others = get_non_spatial_de(
    niches = niches,
    spatial_info = spatial_info,
    type = 'receiver',
    lr_network = lr_network
  )
  receiver_spatial_DE_processed = receiver_spatial_DE_processed %>% bind_rows(receiver_spatial_DE_others)
  
  receiver_spatial_DE_processed = receiver_spatial_DE_processed %>% mutate(scaled_receptor_score_spatial = scale_quantile_adapted(receptor_score_spatial))
  
} else {
  receiver_spatial_DE_processed = get_non_spatial_de(
    niches = niches,
    spatial_info = spatial_info,
    type = 'receiver',
    lr_network = lr_network
  )
  receiver_spatial_DE_processed = receiver_spatial_DE_processed %>% mutate(scaled_receptor_score_spatial = scale_quantile_adapted(receptor_score_spatial))
}

# Calculate ligand activities and infer active ligand-target links
lfc_cutoff = 0.5
specificity_score_targets = 'min_lfc'

DE_receiver_targets = calculate_niche_de_targets(
  seurat_obj = seurat_obj,
  niches = niches,
  lfc_cutoff = lfc_cutoff,
  expression_pct = expression_pct,
  assay_oi = assay_oi
)

DE_receiver_processed_targets = process_receiver_target_de(
  DE_receiver_targets = DE_receiver_targets,
  niches = niches,
  expression_pct = expression_pct,
  specificity_score = specificity_score_targets
)

background = DE_receiver_processed_targets  %>% pull(target) %>% unique()

geneset_niche1 = DE_receiver_processed_targets %>% filter(
  receiver == niches[[1]]$receiver &
    target_score >= lfc_cutoff &
    target_significant == 1 &
    target_present == 1
) %>% pull(target) %>% unique()

geneset_niche2 = DE_receiver_processed_targets %>% filter(
  receiver == niches[[2]]$receiver &
    target_score >= lfc_cutoff &
    target_significant == 1 &
    target_present == 1
) %>% pull(target) %>% unique()

geneset_niche1 %>% setdiff(rownames(ligand_target_matrix))
geneset_niche2 %>% setdiff(rownames(ligand_target_matrix))

length(geneset_niche1)
length(geneset_niche2)

top_n_target = 250

niche_geneset_list = list(
  'infected_niche' = list(
    'receiver' = niches[[1]]$receiver,
    'geneset' = geneset_niche1,
    'background' = background
  ),
  'uninfected_niche' = list(
    'receiver' = niches[[2]]$receiver,
    'geneset' = geneset_niche2 ,
    'background' = background
  )
)

ligand_activities_targets = get_ligand_activities_targets(
  niche_geneset_list = niche_geneset_list,
  ligand_target_matrix = ligand_target_matrix,
  top_n_target = top_n_target
)

# Calculate scaled expression of ligands, receptors and targets across cell types of interest
features_oi = union(lr_network$ligand, lr_network$receptor) %>% 
  union(ligand_activities_targets$target) %>% 
  setdiff(NA)

dotplot = suppressWarnings(Seurat::DotPlot(
  seurat_obj %>% 
    subset(idents = niches %>% unlist() %>% unique()),
  features = features_oi,
  assay = assay_oi
))

exprs_tbl = dotplot$data %>% as_tibble()

exprs_tbl = exprs_tbl %>% rename(
  celltype = id,
  gene = features.plot,
  expression = avg.exp,
  expression_scaled = avg.exp.scaled,
  fraction = pct.exp) %>%
  mutate(fraction = fraction / 100) %>%
  as_tibble() %>%
  select(celltype, gene, expression, expression_scaled, fraction) %>%
  distinct() %>%
  arrange(gene) %>% 
  mutate(gene = as.character(gene))

exprs_tbl_ligand = exprs_tbl %>%
  filter(gene %in% lr_network$ligand) %>%
  rename(
    sender = celltype,
    ligand = gene,
    ligand_expression = expression,
    ligand_expression_scaled = expression_scaled,
    ligand_fraction = fraction
  )

exprs_tbl_receptor = exprs_tbl %>%
  filter(gene %in% lr_network$receptor) %>%
  rename(
    receiver = celltype,
    receptor = gene,
    receptor_expression = expression,
    receptor_expression_scaled = expression_scaled,
    receptor_fraction = fraction
  )

exprs_tbl_target = exprs_tbl %>% 
  filter(gene %in% ligand_activities_targets$target) %>%
  rename(
    receiver = celltype,
    target = gene,
    target_expression = expression,
    target_expression_scaled = expression_scaled,
    target_fraction = fraction
  )

exprs_tbl_ligand = exprs_tbl_ligand %>% 
  mutate(scaled_ligand_expression_scaled = scale_quantile_adapted(ligand_expression_scaled)) %>% 
  mutate(ligand_fraction_adapted = ligand_fraction) %>% 
  mutate_cond(ligand_fraction >= expression_pct, ligand_fraction_adapted = expression_pct) %>% 
  mutate(scaled_ligand_fraction_adapted = scale_quantile_adapted(ligand_fraction_adapted))

exprs_tbl_receptor = exprs_tbl_receptor %>% 
  mutate(scaled_receptor_expression_scaled = scale_quantile_adapted(receptor_expression_scaled)) %>% 
  mutate(receptor_fraction_adapted = receptor_fraction) %>% 
  mutate_cond(receptor_fraction >= expression_pct, receptor_fraction_adapted = expression_pct) %>% 
  mutate(scaled_receptor_fraction_adapted = scale_quantile_adapted(receptor_fraction_adapted))

# Expression fraction and receptor
exprs_sender_receiver = lr_network %>%
  inner_join(exprs_tbl_ligand, by = c('ligand')) %>%
  inner_join(exprs_tbl_receptor, by = c('receptor')) %>% 
  inner_join(DE_sender_receiver %>% distinct(niche, sender, receiver))

ligand_scaled_receptor_expression_fraction_df = exprs_sender_receiver %>%
  group_by(ligand, receiver) %>%
  mutate(
    rank_receptor_expression = dense_rank(receptor_expression),
    rank_receptor_fraction  = dense_rank(receptor_fraction)
  ) %>%
  mutate(ligand_scaled_receptor_expression_fraction = 0.5 * ((
    rank_receptor_fraction / max(rank_receptor_fraction)
  ) +
    ((
      rank_receptor_expression / max(rank_receptor_expression)
    )))) %>%
  distinct(ligand,
           receptor,
           receiver,
           ligand_scaled_receptor_expression_fraction,
           bonafide) %>% distinct() %>% ungroup()

# Prioritization of ligand-receptor and ligand-target links
prioritizing_weights = c(
  'scaled_ligand_score' = 5,
  'scaled_ligand_expression_scaled' = 1,
  'ligand_fraction' = 1,
  'scaled_ligand_score_spatial' = 0,
  'scaled_receptor_score' = 0.5,
  'scaled_receptor_expression_scaled' = 0.5,
  'receptor_fraction' = 1,
  'ligand_scaled_receptor_expression_fraction' = 1,
  'scaled_receptor_score_spatial' = 0,
  'scaled_activity' = 0,
  'scaled_activity_normalized' = 1,
  'bona_fide' = 1
)

output = list(
  DE_sender_receiver = DE_sender_receiver,
  ligand_scaled_receptor_expression_fraction_df = ligand_scaled_receptor_expression_fraction_df,
  sender_spatial_DE_processed = sender_spatial_DE_processed,
  receiver_spatial_DE_processed = receiver_spatial_DE_processed,
  ligand_activities_targets = ligand_activities_targets,
  DE_receiver_processed_targets = DE_receiver_processed_targets,
  exprs_tbl_ligand = exprs_tbl_ligand,
  exprs_tbl_receptor = exprs_tbl_receptor,
  exprs_tbl_target = exprs_tbl_target
)

prioritization_tables = get_prioritization_tables(output, prioritizing_weights)

prioritization_tables$prioritization_tbl_ligand_receptor %>% 
  filter(receiver == niches[[1]]$receiver) %>% head(10)
prioritization_tables$prioritization_tbl_ligand_target %>% 
  filter(receiver == niches[[1]]$receiver) %>% head(10)
prioritization_tables$prioritization_tbl_ligand_receptor %>% 
  filter(receiver == niches[[2]]$receiver) %>% head(10)
prioritization_tables$prioritization_tbl_ligand_target %>% 
  filter(receiver == niches[[2]]$receiver) %>% head(10)

# Visualization of the Differential NicheNet output
top_ligand_niche_df = prioritization_tables$prioritization_tbl_ligand_receptor %>%
  select(niche, sender, receiver, ligand, receptor, prioritization_score) %>%
  group_by(ligand) %>%
  top_n(1, prioritization_score) %>%
  ungroup() %>%
  select(ligand, receptor, niche) %>%
  rename(top_niche = niche)

top_ligand_receptor_niche_df = prioritization_tables$prioritization_tbl_ligand_receptor %>%
  select(niche, sender, receiver, ligand, receptor, prioritization_score) %>%
  group_by(ligand, receptor) %>%
  top_n(1, prioritization_score) %>%
  ungroup() %>%
  select(ligand, receptor, niche) %>%
  rename(top_niche = niche)

ligand_prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>%
  select(niche, sender, receiver, ligand, prioritization_score) %>%
  group_by(ligand, niche) %>%
  top_n(1, prioritization_score) %>%
  ungroup() %>%
  distinct() %>%
  inner_join(top_ligand_niche_df) %>%
  filter(niche == top_niche) %>%
  group_by(niche) %>%
  top_n(20, prioritization_score) %>%
  ungroup()

receiver_oi = 'M2-like_infected'

filtered_ligands = ligand_prioritized_tbl_oi %>%
  filter(receiver == receiver_oi) %>%
  pull(ligand) %>%
  unique()

prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>%
  filter(ligand %in% filtered_ligands) %>%
  select(niche,
         sender,
         receiver,
         ligand,
         receptor,
         ligand_receptor,
         prioritization_score) %>%
  distinct() %>%
  inner_join(top_ligand_receptor_niche_df) %>%
  group_by(ligand) %>%
  filter(receiver == receiver_oi) %>%
  top_n(2, prioritization_score) %>%
  ungroup()

lfc_plot = make_ligand_receptor_lfc_plot(
  receiver_oi,
  prioritized_tbl_oi,
  prioritization_tables$prioritization_tbl_ligand_receptor,
  plot_legend = F,
  heights = NULL,
  widths = NULL
) 

lfc_plot

pdf(paste0('./plots/nichenet/nichenet_DE_lfc_plot_all_', receiver_oi,'.pdf'),
    width = 8,
    height = 9)
plot(lfc_plot)
dev.off()

exprs_activity_target_plot = make_ligand_activity_target_exprs_plot(
  receiver_oi,
  prioritized_tbl_oi,
  prioritization_tables$prioritization_tbl_ligand_receptor,
  prioritization_tables$prioritization_tbl_ligand_target,
  output$exprs_tbl_ligand,
  output$exprs_tbl_target,
  lfc_cutoff,
  ligand_target_matrix,
  plot_legend = T,
  heights = NULL,
  widths = NULL
)

exprs_activity_target_plot$combined_plot

pdf(paste0(
  './plots/nichenet/nichenet_DE_exprs_activity_target_plot_all_',receiver_oi,'.pdf'),
  width = 18,
  height = 10
); plot(exprs_activity_target_plot$combined_plot); dev.off()

filtered_ligands = ligand_prioritized_tbl_oi %>%
  filter(receiver == receiver_oi) %>%
  top_n(15, prioritization_score) %>%
  pull(ligand) %>%
  unique()

prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>%
  filter(ligand %in% filtered_ligands) %>%
  select(niche,
         sender,
         receiver,
         ligand,
         receptor,
         ligand_receptor,
         prioritization_score) %>%
  distinct() %>%
  inner_join(top_ligand_receptor_niche_df) %>%
  group_by(ligand) %>%
  filter(receiver == receiver_oi) %>%
  top_n(3, prioritization_score) %>%
  ungroup()

colors_sender <- c(
  'Keratinocyte_infected' = '#8D7E30',
  'Dendritic cell_infected' = '#3473B8',
  'Endothelial_infected' = '#974221',
  'Fibroblast_infected' = '#4B1C5B',
  'Memory T cell_infected' = '#005D31',
  'Neuron_infected' = '#51B749',
  'Neutrophil_infected' = '#FCCB99',
  'Vascular smooth muscle_infected' = '#C01988',
  'Tissue resident macrophage_infected' = '#59BA5D',
  'M2-like_infected' = '#9b32f5'
)

colors_receiver = c('#9b32f5') %>% #M2 A12BE0
  magrittr::set_names(prioritized_tbl_oi$receiver %>%
                        unique() %>%
                        sort())

transparency <- prioritized_tbl_oi %>% 
  mutate(weight = (prioritization_score - min(prioritization_score)) / (max(prioritization_score) - min(prioritization_score))) %>% 
  mutate(transparency = 1 - prioritization_score) %>% .$transparency

pdf(paste0(
  './plots/nichenet/nichenet_DE_exprs_activity_circos_all_',receiver_oi,'.pdf'),
  width = 10,
  height = 10
) 
circos_output = make_circos_lr(prioritized_tbl_oi, 
                               colors_sender, 
                               colors_receiver,
                               scale = T, border = F,
                               transparency = transparency)
dev.off()

sessionInfo()
