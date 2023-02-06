# Read adata created in 5.create_adata.py
adata = scv.read('adata/infected.h5ad') # use the same script for other condition

# look at the ratio of unspliced/spliced transcripts for each cluster 
scv.pl.proportions(adata, groupby='seurat_clusters', save = 'spliced_unspliced.pdf')

# Subset immune cells for velocity estimation
immune = ['MAC1', 'MAC2', 'DC', 'TMEM']
adata = adata[adata.obs['seurat_clusters'].isin(immune)]

# Use a minimum of 20 shared counts for top 2000 features
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
sc.tl.pca(adata)
sc.pp.neighbors(adata, n_pcs=14, n_neighbors=30)
scv.pp.moments(adata, n_pcs=None, n_neighbors=None)

# recover dynamics for dynamical velocity estimation
scv.tl.recover_dynamics(adata, n_jobs = 10)

scv.tl.velocity(adata, mode = 'dynamical', n_jobs = 10)
scv.tl.velocity_graph(adata, n_jobs = 10)

# save for later
adata.write('processed_adata/infected.h5ad', compression='gzip')
#adata = scv.read('processed_adata/infected.h5ad')

scv.pl.velocity_embedding_stream(adata,
                                 basis = 'umap',
                                 smooth = 0.8,
                                 min_mass = 2.5,
                                 color = 'seurat_clusters',
                                 figsize = (2,6),
                                 palette = cluster_palette,
                                 save = 'immune_velocity_stream.svg'
                                )

import session_info
session_info.show()
