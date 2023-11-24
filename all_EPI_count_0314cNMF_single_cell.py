import os
import pandas as pd
import numpy as np
import scanpy as sc
from cnmf import cNMF
import anndata as ad


adata = sc.read_h5ad('All_EPI_0313.h5ad')
adata.var['mt'] = adata.var_names.str.startswith('MT-')
notMT= adata.var['mt'][adata.var['mt']==False]
adata = adata[:, notMT.index]

var_names=adata.var.features.tolist()
n_vars=len(var_names)
var=pd.DataFrame(data=var_names,index=var_names)
var=var.rename(columns={0: 'index'})
adata1=ad.AnnData(adata.X.astype(np.float32),obs=adata.obs,var=var,dtype='float32')
count_adat_fn='./EPI_counts_without_MT.h5ad'
sc.write(count_adat_fn, adata1)


numiter=20 # Number of NMF replicates. Set this to a larger value ~200 for real data. We set this to a relatively low value here for illustration at a faster speed
numhvgenes=2000 ## Number of over-dispersed genes to use for running the actual factorizations
## Results will be saved to [output_directory]/[run_name] which in this example is example_PBMC/cNMF/pbmc_cNMF
output_directory = './'
if not os.path.exists(output_directory):
    os.mkdir(output_directory)
run_name = 'EPI_cNMF'

## Specify the Ks to use as a space separated list in this case "5 6 7 8 9 10"
K = ' '.join([str(i) for i in range(3,26)])

## To speed this up, you can run it for only K=7-8 with the option below
#K = ' '.join([str(i) for i in range(7,9)])


seed = 14 ## Specify a seed pseudorandom number generation for reproducibility
countfn='./EPI_counts_without_MT.h5ad'

cnmf_obj = cNMF(output_dir=output_directory, name=run_name)
cnmf_obj.prepare(counts_fn=countfn, components=np.arange(3,26), n_iter=20, seed=14, num_highvar_genes=2000)
cnmf_obj.factorize(worker_i=0, total_workers=1)
cnmf_obj.combine()
cnmf_obj.k_selection_plot(close_fig=False)

selected_K = 8
density_threshold = 2.00
cnmf_obj.consensus(k=selected_K, density_threshold=density_threshold, show_clustering=True, close_clustergram_fig=False)
density_threshold = 0.10
cnmf_obj.consensus(k=selected_K, density_threshold=density_threshold, show_clustering=True, close_clustergram_fig=False)

adata = sc.read(countfn)
hvgs = open('./EPI_cNMF/EPI_cNMF.overdispersed_genes.txt').read().split('\n')
sc.pp.normalize_per_cell(adata, counts_per_cell_after=10**4) 
adata.raw = sc.pp.log1p(adata.copy(), copy=True)
adata = adata[:,hvgs]
sc.pp.scale(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata, n_neighbors=50, n_pcs=15)
sc.tl.umap(adata)

usage_norm, gep_scores, gep_tpm, topgenes = cnmf_obj.load_results(K=selected_K, density_threshold=density_threshold)
usage_norm.columns = ['Usage_%d' % i for i in usage_norm.columns]
usage_file = cnmf_obj.paths['consensus_usages__txt'] % (selected_K, '0_1')
print(usage_file)

gene_scores_file = cnmf_obj.paths['gene_spectra_score__txt'] % (selected_K, '0_1')
print(gene_scores_file)

gene_tpm_file = cnmf_obj.paths['gene_spectra_tpm__txt'] % (selected_K, '0_1')
print(gene_tpm_file)
adata.obs = pd.merge(left=adata.obs, right=usage_norm, how='left', left_index=True, right_index=True)






