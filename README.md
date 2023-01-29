# Celik_et_al_2023
Scripts used in:

Celik C, XXX

# Pre-processing single-cell rawdata

This [script](cellranger_scripts.md) contains preprocessing and obtaining counts, .BAM files and .cloupe files from the fastq files using Cell Ranger.

```cellranger count``` will generate three matrices (barcodes.tsv, features.csv and matrix.mtx). These will be fed into [Seurat](https://satijalab.org/seurat/) in R to create a ```SeuratObject```.
```.BAM``` files will later be used for creating ```.loom``` files required for RNA velocity analysis using [Velocyto](http://velocyto.org/velocyto.py/tutorial/analysis.html#analysis) and [scVelo](http://scvelo.readthedocs.io).

# Integration of data

This ```R``` [script](seurat_preprocessing.r) will take the Cell Ranger output ```filtered_feature_bc_matrix``` data as input, filter out doublet GEMs and low quality cells, and integrate the two conditions in the Seurat package (v4.3.0) using ```SCTransform``` function.

## Identifying clusters

After ```SCTransform```, the number of principal components and resolution for UMAP reduction will be estimated using probabilistic approaches as described [here](clusters.md).

## Cell annotation

Using the integrated object, annotate the cell populations using the extensive dataset for mouse skin cell population markers published in [Joost _et al._ (2015)](https://doi.org/10.1016/j.cels.2016.08.010) and [Joost _et al._ (2020)](https://doi.org/10.1016/j.stem.2020.01.012) using this [script](annotation.md).

# RNA velocity

To determine the terminal cell populations, first sort the ```.BAM``` files created in Cell Ranger, then run [velocyto](http://velocyto.org/velocyto.py/tutorial/analysis.html) to create ```.loom``` files with _spliced_ and _unspliced_ transcript matrix. Merge the loom files in _Python_ for [scVelo](https://scvelo.readthedocs.io/en/stable/). Finally, utilise [cellrank](http://cellrank.readthedocs.io) for estimating terminal cell states. The entire pipeline can be found [here].

# Terminal cell fates

# Cell-cell interactions

