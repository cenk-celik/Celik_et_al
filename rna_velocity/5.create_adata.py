import scanpy as sc
import anndata as ad
from scipy import io
import numpy as np
import pandas as pd
import scvelo as scv
import cellrank as cr

# set scVelo settings
scv.settings.verbosity = 2
scv.settings.presenter_view = True
scv.set_figure_params(facecolor='white',
                               dpi_save=600,
                               frameon=False,
                               vector_friendly=True,
                               figsize=(4,6),
                               transparent=True)
scv.settings.plot_prefix = 'infected/infected_'
cr.settings.verbosity = 2

# To match ggplot2 colours in Seurat:
main_celltype_palette={'Epithelial':'#F6766D','Fibroblast':'#A2A402',
                       'Immune':'#00BD7C','Neurons':'#00AEF4','Vascular':'#E46AF1'}

colour_palette={'Basal':"#F0A0FF", 'Dendritic cells':"#0075DC",
                'Endothelial cells':"#993F00", 'Fibroblast':"#4C005C",
                'Macrophages':"#191919", 'Memory T cells':"#005C31",
                'Neurons':"#2BCE48", 'Neutrophils':"#FFCC99",
                'Outer bulge': "#808080", 'Sebacious gland':"#94FFB5",
                'Suprabasal':"#8F7C00",
                'Upper hair follicle suprabasal':"#9DCC00",
                'Vascular smooth muscle':"#C20088"}

cluster_palette={'MAC1':"#aa0dfe", 'FIB1':"#3283fe",
                 'FIB2':"#85660d", 'BAS1':"#782ab6",
                 'SUP1':"#565656",'HFSUP1':"#1c8356",
                 'SG':"#16FF32", 'NEUT1':"#F7E1A0",
                 'SUP2':"#E2E2E2", 'TRM':"#1CBE4F",
                 'SUP3':"#C4451C", 'EC1':"#DEA0FD",
                 'EC2':"#FE00FA", 'SUP4':"#325A9B",
                 'OB':"#FEAF16", 'FIB3':"#F8A19F",
                 'HFSUP2':"#90AD1C", 'DC':"#F6222E",
                 'VSM':"#1CFFCE", 'NEU':"#2ED9FF",
                 'NEUT2':"#B10DA1", 'BAS3':"#C075A6",
                 'TMEM':"#FC1CBF", 'BAS2':"#B00068"}

# read count matrix for the condition of interest
## here I create an Anndata for intected samples first.
## the same script can be used for the other condition
X = io.mmread("data/infected_counts.mtx")

adata = ad.AnnData(
    X=X.transpose().tocsr(), dtype=X.dtype
)

# load cell metadata pre-calculated in Seurat:
cell_meta = pd.read_csv("data/infected_metadata.csv")

# load gene names pre-calculated in Seurat:
with open("data/infected_gene_names.csv", 'r') as f:
    gene_names = f.read().splitlines()

# load barcodes pre-calculated in Seurat:
adata.obs = cell_meta
adata.obs.index = adata.obs['barcode']
adata.var.index = gene_names

# load dimensional reduction pre-calculated in Seurat:
pca = pd.read_csv("data/infected_pca.csv")
pca.index = adata.obs.index

# set pca and umap pre-calculated in Seurat:
adata.obsm['X_pca'] = pca.to_numpy()
adata.obsm['X_umap'] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T

# plot a UMAP colored by celltypes to test:
sc.pl.umap(adata,
           color='seurat_clusters',
           frameon=False,
           palette=cluster_palette)

# read .loom files created in velocyto for each sample in the integrated dataset
ldata1 = scv.read('loom/infected_1.loom', cache=False, validate=False)
ldata2 = scv.read('loom/infected_2.loom', cache=False, validate=False)

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

# clean up barcodes
scv.utils.clean_obs_names(adata)
scv.utils.clean_obs_names(ldata)

# merge matrices into the original adata object
adata = scv.utils.merge(adata, ldata)

# save dataset as anndata format
adata.write('adata/infected.h5ad')

import session_info
session_info.show()
