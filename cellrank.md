
Using the adata object saved from ```scVelo``` analysis, run ```Cell Rank``` to estimate cluster fates and and genes that drive the differentiation.

```python
adata = scv.read('infected.h5ad')

# run CellRank
cr.tl.terminal_states(adata, cluster_key = 'seurat_clusters', weight_connectivities = 0.1)
cr.pl.terminal_states(adata, save='cellrank_terminal_states.pdf', figsize = (5,5))
```

Iidentify initial states:

```python
cr.tl.initial_states(adata, cluster_key = 'seurat_clusters')
cr.pl.initial_states(adata, discrete = True, save = 'cellrank_initial_states.pdf', figsize = (5,5))
```

Compute fate maps:

```python
cr.tl.lineages(adata, tol=1e-6)
cr.pl.lineages(adata, same_plot = False, ncols = 4, figsize = (5,5), save = 'cellrank_fatemaps.pdf')

cr.pl.cluster_fates(adata, mode = 'paga_pie', cluster_key = 'seurat_clusters', basis = 'umap',
                   legen_kwargs = {'loc': 'top right out'}, legend_loc = 'top left out', node_size_scale = 5,
                   edge_width_scale = 1, max_edge_width = 4, title = 'Directed PAGA', save = 'cellrank_directed_paga.pdf',
                   figsize = (5,5))
```

Compute lineage drivers for M1-like macrophages:

```python
cr.tl.lineage_drivers(adata)
#adata.varm['terminal_lineage_drivers']['MAC2_pval'].isnull().sum()

cr.pl.lineage_drivers(adata, lineage = 'MAC2', n_genes = 20, ncols = 2, figsize = (10,50),
                     save = 'cellrank_mac2_lineage_drivers.pdf')
```

Gene expression trends in M1-like macrophages:

```python
root_idx = np.where(adata.obs['initial_states'] == 'HFSUP1')[0][0]

adata.uns['iroot'] = root_idx
sc.tl.dpt(adata)

scv.pl.scatter(adata, color = ['seurat_clusters', root_idx, 'latent_time', 'dpt_pseudotime'], fontsize = 16,
              cmap = 'viridis', perc = [2, 98], colorbar = True, rescale_color = [0, 1],
              title = ['seurat clusters', 'root cell', 'latent time', 'dpt pseudotime'], ncols = 2,
              save = 'cellrank_gene_expression_trends.pdf', figsize = (5, 5))

model = cr.ul.models.GAM(adata)

genes = ['Ctss', 'Lyz2', 'Psap', 'Pf4', 'Ctsb', 'Arg1', 'Ftl1', 'Npc2', 'Laptm5', 'Ctsl', 'Fcgr2b', 'Lgmn', 'Abca1',
        'F13a1', 'Ccl9', 'Csf1r', 'Grn', 'Thbs1', 'Apoe', 'Ms4a6d'] # selected baesd on gene expression trends

cr.pl.gene_trends(adata, model = model, data_key = 'X', genes = genes,
                 ncols = 2, time_key = 'latent_time', same_plot = True, hide_cells = True, n_test_points = 200,
                 save = 'cellrank_mac2_gene_trends.pdf', show_progress_bar = False)
```

Latent time heat map in M1-like macrophages:

```python
cr.pl.heatmap(adata, model, 
              genes = adata.varm['terminal_lineage_drivers']['MAC2_corr'].sort_values(ascending = False).index[:100],
              show_absorption_probabilities = True, lineages = 'MAC2', n_jobs = 1, backend = 'loky',
              save = 'cellrank_mac2_terminal_lineage_drivers_latent_time_heatmap.pdf', figsize = (15,6))
```
