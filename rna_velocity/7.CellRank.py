# Compute terminal states
cr.tl.terminal_states(adata,
                      cluster_key="seurat_clusters",
                      weight_connectivities=0.2)

cr.pl.terminal_states(adata,
                      figsize = (2,6),
                      save = 'figures/infected/infected_immune_cellrank_terminal_states.pdf'
                      )

# Compute initial states
cr.tl.initial_states(adata,
                     cluster_key = 'seurat_clusters')

cr.pl.initial_states(adata,
                     discrete = True,
                     figsize = (2,6),
                     save = 'figures/infected/infected_immune_cellrank_initial_states.pdf'
                    )

# Compute fate maps
cr.tl.lineages(adata)

cr.pl.lineages(adata,
               same_plot = False,
               figsize = (2,6),
               save = 'figures/infected/infected_immune_cellrank_lineages.pdf'
              )

# Recover latent time for directed PAGA
scv.tl.recover_latent_time(adata,
                           root_key="initial_states_probs",
                           end_key="terminal_states_probs")

scv.tl.paga(adata,
            groups="seurat_clusters",
            root_key="initial_states_probs",
            end_key="terminal_states_probs",
            use_time_prior="velocity_pseudotime")

cr.pl.cluster_fates(adata,
                    mode="paga_pie",
                    cluster_key="seurat_clusters",
                    basis="umap",
                    legend_kwargs={"loc": "top right out"},
                    legend_loc="top left out",
                    node_size_scale=5,
                    edge_width_scale=1,
                    max_edge_width=4,
                    figsize = (2.5,6),
                    title="directed PAGA",
                    save = 'infected/infected_immune_cluster_fates.svg'
                   )

import session_info
session_info.show()
