# Compute lineage drivers for M2-like macrophages
cr.tl.lineage_drivers(adata,lineages='MAC1')

# plot top 12 driver genes
cr.pl.lineage_drivers(adata,
                      lineage = 'MAC1',
                      n_genes = 12,
                      ncols = 3,
                      figsize = (6,24),
                      title_fmt='{gene} qval={qval:.2e}',
                      save = 'infected/infected_immune_lineage_drivers_rescaled.pdf',
                      color_map = 'Reds',
                      smooth = True,
                      wspace = 5,
                     )

# Compute gene expression trends
## compute DPT, starting from CellRank defined root cell
root_idx = np.where(adata.obs["initial_states"] == "MAC2")[0][0]
adata.uns["iroot"] = root_idx
sc.tl.diffmap(adata)
sc.tl.dpt(adata)

scv.pl.scatter(adata,
               color = ["seurat_clusters", root_idx, "latent_time", "dpt_pseudotime"],
               fontsize = 16,
               cmap = "viridis",
               perc = [2, 98],
               colorbar = True,
               rescale_color = [0, 1],
               figsize = (2,6),
               title = ["clusters", "root cell", "latent time", "dpt pseudotime"],
               save = 'immune_root_latent_dpt.pdf'
              )

# Plot gene trends
model = cr.ul.models.GAM(adata)

genes = adata.varm['terminal_lineage_drivers']["MAC1_corr"].sort_values(ascending=False).index[:12]

cr.pl.gene_trends(adata,
                  model = model,
                  data_key = "Ms",
                  genes = genes,
                  lineages = 'MAC1',
                  time_key = "latent_time",
                  same_plot = False,
                  hide_cells = False,
                  show_progress_bar = False,
                  figsize = (3, 15),
                  gene_as_title = True,
                  cell_color = 'latent_time',
                  n_jobs = 10,
                  time_range=[(0,0.6)],
                  save = 'infected/infected_immune_gene_trends.pdf'
                 )

# plot latent time driver genes' heat map
cr.pl.heatmap(adata,
              model,
              adata.var_names[:30],
              time_key = 'latent_time',
              lineages = 'MAC1',
              show_progress_bar = False,
              show_absorption_probabilities = True,
              figsize = (10, 10),
              save ='infected/infected_immune_latent_time.pdf',
              return_genes = False)

# plot cluster fates
cr.pl.cluster_fates(adata, mode = "violin",
                    cluster_key = "seurat_clusters",
                    lineages = 'MAC1',
                    figsize = (3,2),
                    save = 'infected/infected_immune_aggregated_cluster_fates.pdf'
                   )

import session_info
session_info.show()
