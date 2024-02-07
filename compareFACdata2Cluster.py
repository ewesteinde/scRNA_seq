# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 16:23:37 2024

@author: ewest
"""

# download Smart-seq2 data from ageing brain dataset
# clean reads using fastp - already done??
# Map to BDGP6 D. melanogaster and quantified using STAR -already done??
# Create pseudo-bulk data from cluster marker genes. 
# Create some kind of averaged expression from dfb cells? Look into FAC sorted cell analysis standards
# use a non-negative leaseter squares regression model to match dFB Smart seq2 cells to their best corresponding cluster
# "Cluster-to-bulk mapping was performed through a non-negative least-squares regression model (NNLS)"
# "Before fitting, all datasets were scaled to CPM using the common genes, independent of previously applied normalization techniques"
# to read tsv files: pd.read_csv('file,sep='\t') 
# reads of each gene for each cell 

import pandas as pd
import scanpy as sc

import os
import platform 

# replace paths to code directory if needed
if platform.system() == 'Windows':
    os.chdir('Z:\Dropbox (HMS)\Wilson_Lab_Data\Omics_datasets\processedData')
    saveDir = 'Z:\Dropbox (HMS)\Wilson_Lab_Data\Omics_datasets\processedData'
    file = 'Z:/Dropbox (HMS)/Wilson_Lab_Data/Omics_datasets/datasets/AgingBrain/R23E10_gal4/GSE107451_R23E10_SMART-seq2_DGEM.tsv'
    adata = sc.read_h5ad('Z:/Dropbox (HMS)/Wilson_Lab_Data/Omics_datasets/processedData/allconds_allgenes_clusterRes_1.5.h5ad')
    os.chdir('C:\Code\scRNA_seq')
else:
    os.chdir('/Users/elenawesteinde/Documents/Omics_proj/scRNA_seq')
    saveDir = '/Users/elenawesteinde/Dropbox (HMS)/Wilson_Lab_Data/Omics_datasets/processedData'
    file = '/Users/elenawesteinde/Dropbox (HMS)/Wilson_Lab_Data/Omics_datasets/datasets/AgingBrain/R23E10_gal4/GSE107451_R23E10_SMART-seq2_DGEM.tsv'
    adata = sc.read_h5ad('/Users/elenawesteinde/Dropbox (HMS)/Wilson_Lab_Data/Omics_datasets/processedData/allconds_allgenes_postCluster_normBeforeCocat.h5ad')
    markers = pd.read_csv('/Users/elenawesteinde/Dropbox (HMS)/Wilson_Lab_Data/Omics_datasets/processedData/Markers_allconds_allgenes_clusterRes_1.5.csv')

#dfb = pd.read_csv(file,sep = '\t') 
dfb = sc.read_csv(file,  delimiter='\t')
dfb = dfb.transpose()


#%% look at dfb data

# keep genes shared b/w datasets 
# lognormalize dfb cells
# coefficients are each reference transcriptome --> each cluster 

#%% Strategy 1 (ageing & sleep dataset)

# non linear least squares regression
# scipy.optimize.nnls(A, b, maxiter=None, *, atol=None)

# n_samples: cells
# n_features: transcripts/genes
# X: dfb(n_samples, n_features)

# Compare against pseudobulk data of each cluster (1, n_features)
# find which cluster the dfb cells have best match with 

#%% Strategy 2 map query dataset (dfb) onto ref dataset (adata - whole brain)

# only keep shared features
var_names = adata.var_names.intersection(dfb.var_names)
adata = adata[:,var_names]
dfb = dfb[:, var_names]

#%% recluster reference dataset 
import clusterFunctions as cf
adata, model, cluster_dic, markers_dict = cf.clusterData(adata, 'all', 1)

sc.tl.ingest(dfb, adata, obs='leiden', embedding_method = 'umap')

dfb.uns['leiden_colors'] = adata.uns['leiden_colors']  # fix colors
sc.pl.umap(dfb, color=['leiden'], wspace=0.5)
adata_concat = adata.concatenate(dfb, batch_categories = ['ref','new'])
adata_concat.obs.leiden = adata_concat.obs.leiden.astype('category')














