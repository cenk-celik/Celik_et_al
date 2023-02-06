## Preprocess data for scVelo

First, extract meta data using the integrated ```SeuratObject``` in ```R``` for each sample:

[1. Extract data](1.extract_barcodes.r)

Two datasets are needed for each sample for creating ```.loom``` files:

- barcodes.tsv extracted from your processed Seurat object
- Raw .BAM file.

In shell, run the following script to sort the data based on barcodes for each sample:

[2. Sort barcodes in ```.BAM``` files](2.sort_barcodes.sh)

After sorting is completed, move the sorted .BAM file into the same output folder for that sample with ```cellsorted_[original_bam_name].bam```. 
Both unsorted and sorted .BAM files must be in the same folder. Velocyto will skip sorting ```[original_bam_name].bam```.

To generate ```.loom``` files with spliced and unspliced genes, run ```velocyto``` in shell with the following ```repeats mask``` (mm10_rmsk.gtf)
and mouse reference annotation file used during the alignment in Cell Ranger (refdata-gex-mm10-2020-A).

- Repeats mask file (mouse: [link](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=611454127_NtvlaW6xBSIRYJEBI0iRDEWisITa&clade=mammal&org=Mouse&db=mm10&hgta_group=allTracks&hgta_track=rmsk&hgta_table=0&hgta_regionType=genome&position=chr12%3A56694976-56714605&hgta_outputType=primaryTable&hgta_outputType=gff&hgta_outFileName=mm10_rmsk.gtf); human: [link](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=611454127_NtvlaW6xBSIRYJEBI0iRDEWisITa&clade=mammal&org=Human&db=0&hgta_group=allTracks&hgta_track=rmsk&hgta_table=rmsk&hgta_regionType=genome&position=&hgta_outputType=gff&hgta_outFileName=mm10_rmsk.gtf))
- Reference .GTF annotation file (refdata-gex-mm10-2020-A)

[3. Generate ```.loom``` files for each sample](3.velocyto_loom.sh)

In ```R```, extract more information from the integrated ```SeuratObject``` required for ```scVelo``` for each condition:

[4. Extract PCA and UMAP reduction informations, count matrices and feature names](4.scVelo_preparation.r)

Merge all data generated above to create ```adata``` for each condition in ```Phyton```:

Count matrices can be downloaded from [here](https://doi.org/10.5281/zenodo.7608772). PCA and UMAP reductions and feature name matrices can be found [here](data/).

[5. Create adata for each condition](5.create_adata.py)

## Run scVelo

All the objects are ready required for ```scVelo``` and ```Cellrank```. First, we run ```scVelo``` to estimate velocity using dynamical model:

[6. Run scVelo](6.scVelo_dynamical.py)

The data is now ready for estimating terminal cell states using ```CellRank```:

[7. Identify terminal and initial states](7.CellRank.py)

Compute lineage drivers in ```CellRank```:

[8. Compute lineage drivers and visualise top lineage driving genes for a cluster of interest](8.lineage_drivers.py)

More can be found in [this Jupyter notebook](infected_rna_velocity_complete.ipynb).
