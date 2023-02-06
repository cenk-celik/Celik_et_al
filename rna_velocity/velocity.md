## Preprocess data for scVelo

First, extract meta data using the integrated Seurat Object in ```R``` for each sample:

```R
split.list <- SplitObject(seurat, split.by = "samples")

uninfected.1 <- split.list[[1]]
uninfected.2 <- split.list[[2]]
infected.1 <- split.list[[3]]
infected.2 <- split.list[[4]]

rm(split.list); gc()

uninfected.1.barcodes <- gsub('_1', '', colnames(uninfected.1))
write.csv(uninfected.1.barcodes, file = 'uninfected_1_barcodes.tsv')

uninfected.2.barcodes <- gsub('_2', '', colnames(uninfected.2))
write.csv(uninfected.2.barcodes, file = 'uninfected_2_barcodes.tsv')

infected.1.barcodes <- gsub('_3', '', colnames(infected.1))
write.csv(infected.1.barcodes, file = 'infected_1_barcodes.tsv')

infected.2.barcodes <- gsub('_4', '', colnames(infected.2))
write.csv(infected.2.barcodes, file = 'infected_2_barcodes.tsv')
```

Two datasets are needed for each sample for creating ```.loom``` files:

- barcodes.tsv extracted from your processed Seurat object
- Raw .BAM file.

In shell, run the following command to sort the data based on barcodes:

```
samtools sort -t CB -O BAM -o cellsorted_possorted_genome_bam.bam possorted_genome_bam.bam
```

After sorting is completed, move the sorted .BAM file into the same output folder for that sample with ```cellsorted_[original_bam_name].bam```. 
Both unsorted and sorted .BAM files must be in the same folder. Velocyto will skip sorting ```[original_bam_name].bam```.

To generate ```.loom``` files with spliced and unspliced genes, run ```velocyto``` in shell with the following ```repeats mask``` (mm10_rmsk.gtf)
and mouse reference annotation file used during the alignment in Cell Ranger (refdata-gex-mm10-2020-A).

- Repeats mask file (mouse: [link](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=611454127_NtvlaW6xBSIRYJEBI0iRDEWisITa&clade=mammal&org=Mouse&db=mm10&hgta_group=allTracks&hgta_track=rmsk&hgta_table=0&hgta_regionType=genome&position=chr12%3A56694976-56714605&hgta_outputType=primaryTable&hgta_outputType=gff&hgta_outFileName=mm10_rmsk.gtf); human: [link](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=611454127_NtvlaW6xBSIRYJEBI0iRDEWisITa&clade=mammal&org=Human&db=0&hgta_group=allTracks&hgta_track=rmsk&hgta_table=rmsk&hgta_regionType=genome&position=&hgta_outputType=gff&hgta_outFileName=mm10_rmsk.gtf))
- Reference .GTF annotation file (refdata-gex-mm10-2020-A)

```Bash
velocyto run -b path_to/barcodes_1.tsv \ # filtered barcodes extracted from Seurat object
						 -o path_to/sample/outs \ # output dir for .loom file
						 -m path_to/mm10_rmsk.gtf \ # repeats mask annotation file
						 path_to/possorted_genome_bam.bam # input .BAM
						 path_to/refdata-gex-mm10-2020-A/genes/genes.gtf \ # reference .GTF annotation file
```

In ```R```, extract more metadata from the integrated Seurat Object:

```R
# save metadata table:
seurat$barcode <- colnames(seurat)
seurat$UMAP_1 <- seurat@reductions$umap@cell.embeddings[,1]
seurat$UMAP_2 <- seurat@reductions$umap@cell.embeddings[,2]
write.csv(seurat@meta.data, file='integrated_metadata.csv', quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(seurat, assay='RNA', slot='counts')
writeMM(counts_matrix, file='integrated_counts.mtx')

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(seurat@reductions$pca@cell.embeddings, file='_integrated_pca.csv', 
          quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='integrated_gene_names.csv',
  quote=F,row.names=F,col.names=F
)
```

Merge all data generated above to create ```adata``` for each condition in ```Phyton```.

```Python
import scanpy as sc
import anndata as ad
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import os
import pandas as pd

X = io.mmread("integrated_counts.mtx")

# create anndata object
adata = ad.AnnData(
    X=X.transpose().tocsr()
)

# load cell metadata pre-calculated in Seurat:
cell_meta = pd.read_csv("integrated_metadata.csv")

# load gene names pre-calculated in Seurat:
with open("integrated_gene_names.csv", 'r') as f:
    gene_names = f.read().splitlines()

# load barcodes pre-calculated in Seurat:
adata.obs = cell_meta
adata.obs.index = adata.obs['barcode']
adata.var.index = gene_names

# load dimensional reduction pre-calculated in Seurat:
pca = pd.read_csv("integrated_pca.csv")
pca.index = adata.obs.index

# set pca and umap pre-calculated in Seurat:
adata.obsm['X_pca'] = pca.to_numpy()
adata.obsm['X_umap'] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T

# plot a UMAP colored by celltypes to test:
sc.pl.umap(adata, color=['celltype'], frameon=False, save='umap_uninfected_celltype.svg') #celltype is the cell annotation slot in the Seurat object

# save dataset as anndata format
adata.write('seurat.h5ad')
```

Create colour palettes for matching colour schemes for the plots generated in ```R```:

```python
main_celltype_palette=['#F6766D', '#A2A402', '#00BD7C', '#00AEF4', '#E46AF1']
colour_palette=["#F0A0FF", "#0075DC", "#993F00", "#4C005C", "#191919", 
                "#005C31", "#2BCE48", "#FFCC99", "#808080", "#94FFB5",
                "#8F7C00", "#9DCC00", "#C20088" ]
cluster_palette=["#AA0DFE", "#3283FE", "#85660D", "#782AB6", "#565656",
                 "#1C8356", "#16FF32", "#F7E1A0", "#E2E2E2", "#1CBE4F",
                 "#C4451C", "#DEA0FD", "#FE00FA", "#325A9B", "#FEAF16",
                 "#F8A19F", "#90AD1C", "#F6222E", "#1CFFCE", "#2ED9FF",
                 "#B10DA1", "#C075A6", "#FC1CBF", "#B00068"]
```

Set parameters for scVelo plots:

```python
import scvelo as scv
import cellrank as cr
```

```python
scv.settings.verbosity = 2
scv.settings.presenter_view = True
scv.settings.set_figure_params(facecolor='white', dpi_save=300, frameon=False, vector_friendly=True,
                              figsize=(5,5), format='svg', transparent=True)
scv.settings.plot_prefix = 'infected_'

cr.settings.verbosity = 2
```

Load ```adata```:
```python
adata = sc.read_h5ad('infected_seurat.h5ad')
```

Read ```.loom``` files created using velocyto for each sample in the integrated dataset:

```python
ldata1 = scv.read('/Users/cenk/nmrc/data/infected_1/outs/possorted_genome_bam.loom', cache=False, validate=False)
ldata2 = scv.read('/Users/cenk/nmrc/data/infected_2/outs/possorted_genome_bam.loom', cache=False, validate=False)
```

Rename barcodes in order to merge the samples for a condition

```python
barcodes = [bc.split(':')[1] for bc in ldata1.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_1' for bc in barcodes]
ldata1.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata2.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_2' for bc in barcodes]
ldata2.obs.index = barcodes

ldata1.var_names_make_unique()
ldata2.var_names_make_unique()

# concatenate loom data
ldata = ldata1.concatenate(ldata2)
```

Merge matrices into the original adata object:

```python
adata = scv.utils.merge(adata, ldata)
```

Check if UMAP is identical to that of created in Seurat:

```python
sc.pl.umap(adata, color='seurat_clusters', frameon=False, legend_loc='on data', 
           title='', palette=cluster_palette)
```

Visualise spliced/unspliced proportions for the cell 2nd level cell populations:

```python
scv.pl.proportions(adata, groupby='celltype_2nd', save = 'spliced_unspliced.svg')
```

Filter and normalise the data:

```python
scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)
```

Compute velocity using ```dynamical``` mode in scVelo:

```python
scv.tl.recover_dynamics(adata)
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)
```

Estimate latent time and top genes:

```python
df = adata.var
df = df[(df['fit_likelihood'] > .1) & df['velocity_genes'] == True]

kwargs = dict(xscale='log', fontsize=16)
with scv.GridSpec(ncols=3) as pl:
    pl.hist(df['fit_alpha'], xlabel='transcription rate', **kwargs)
    pl.hist(df['fit_beta'] * df['fit_scaling'], xlabel='splicing rate', xticks=[.1, .4, 1], **kwargs)
    pl.hist(df['fit_gamma'], xlabel='degradation rate', xticks=[.1, .4, 1], **kwargs)

scv.get_df(adata, 'fit*', dropna=True).head()

scv.tl.latent_time(adata)
scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80, save='scatter_latent_time.svg')

top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
scv.pl.heatmap(adata, var_names=top_genes, sortby='latent_time', col_color='seurat_clusters', n_convolve=100,
              save='overall_latent_time.pdf')
```


Visualise top likelihood genes

```python
kwargs = dict(frameon=False, size=10, linewidth=1.5)

top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index
scv.pl.scatter(adata, basis=top_genes[:15], ncols=5, **kwargs, save='top_likelihood_genes.svg',
              color='seurat_clusters')
```

Cluster specific top likelihood genes:

```python
scv.tl.rank_dynamical_genes(adata, groupby='seurat_clusters')
df = scv.get_df(adata, 'rank_dynamical_genes/names')
df.head(5)
```

Identify important genes:
```python
scv.tl.rank_velocity_genes(adata, groupby='seurat_clusters', min_corr=.3)
df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
df.head()
```

## Cell Rank


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
