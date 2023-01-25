# Celik_2023
Scripts used in:

Celik C, XXX

# Pre-processing single-cell rawdata

This [script](cellranger_scripts.md) contains preprocessing and obtaining counts, .BAM files and .cloupe files from the fastq files using Cell Ranger.

```cellranger count``` will generate three matrices (barcodes.tsv, features.csv and matrix.mtx). These will be fed into [Seurat](https://satijalab.org/seurat/) in R to create a ```SeuratObject```.
```.BAM``` files will later be used for creating ```.loom``` files required for RNA velocity analysis using [Velocyto](http://velocyto.org/velocyto.py/tutorial/analysis.html#analysis) and [scVelo](http://scvelo.readthedocs.io).

# Integration of data



## Identifying clusters

## Cell annotation

## Gene Ontology

# RNA velocity

# Terminal cell fates

# Cell-cell interactions

