velocyto run -b path_to/barcodes_1.tsv \ # filtered barcodes extracted from Seurat object
						 -o path_to/sample/outs \ # output dir for .loom file
						 -m path_to/mm10_rmsk.gtf \ # repeats mask annotation file
						 path_to/possorted_genome_bam.bam # input .BAM
						 path_to/refdata-gex-mm10-2020-A/genes/genes.gtf \ # reference .GTF annotation file
