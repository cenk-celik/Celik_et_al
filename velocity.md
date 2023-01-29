First, extract meta data using the integrated Seurat Object in ```R``` for each sample:

```R
split.list <- SplitObject(seurat, split.by = "samples")

uninfected.1 <- split.list[[1]]
uninfected.2 <- split.list[[2]]
infected.1 <- split.list[[3]]
infected.2 <- split.list[[4]]

rm(split.list); gc()

uninfected.1.barcodes <- gsub('_1', '', colnames(uninfected.1))
write.csv(uninfected.1.barcodes, file = '20220929_uninfected_1_barcodes.tsv')

uninfected.2.barcodes <- gsub('_2', '', colnames(uninfected.2))
write.csv(uninfected.2.barcodes, file = '20220929_uninfected_2_barcodes.tsv')

infected.1.barcodes <- gsub('_3', '', colnames(infected.1))
write.csv(infected.1.barcodes, file = '20220929_infected_1_barcodes.tsv')

infected.2.barcodes <- gsub('_4', '', colnames(infected.2))
write.csv(infected.2.barcodes, file = '20220929_infected_2_barcodes.tsv')
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

import scvelo as scv
import cellrank as cr

# read .loom files created using velocyto for each sample in the integrated dataset
ldata1 = scv.read('path_to/possorted_genome_bam.loom', cache=True, validate=False) # sample_1
ldata2 = scv.read('path_to/possorted_genome_bam.loom', cache=True, validate=False) # sample_2, create as many as in your Seurat object

# rename barcodes in order to merge
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

# merge matrices into the original adata object
adata = scv.utils.merge(adata, ldata)

# check if UMAP is identical to that of created in Seurat
sc.pl.umap(adata, color='celltype', frameon=False, legend_loc='on data', 
           title='', save='umap.svg')
```
