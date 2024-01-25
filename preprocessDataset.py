# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 11:56:21 2024

@author: ewest
"""

import os
import platform 
# replace paths to code directory if needed
if platform.system() == 'Windows':
    os.chdir('Z:\Dropbox (HMS)\Wilson_Lab_Data\Code\OmicsCode')
else:
    os.chdir('/Users/elenawesteinde/Dropbox (HMS)/Wilson_Lab_Data/Code/OmicsCode')

import preprocessFunctions as ppf 
import scanpy as sc 
from scipy.sparse import csr_matrix
sc.set_figure_params(dpi=100, dpi_save=100)

#%% gather data
dataFiles, conds = ppf.gatherFiles(rootDir = 'Z:\Dropbox (HMS)\Wilson_Lab_Data\Omics_datasets\GSE207799_RAW')
uniqueConds = set(conds)
print(uniqueConds)
#%% preprocess data

# preprocesses each individual dataset seperately and returns a list with all of them
allData = [] 
for file in dataFiles:
    allData.append(ppf.preprocess(file))
    
#%% Cocatinate datasets

cond = 'all'
if cond == 'all':
    datasetIndices = list(range(0,len(conds)))
elif cond == 'test':
    datasetIndices = [0,1]
else:
    datasetIndices = [i for i in range(len(conds)) if conds[i] == cond]

# join = outer causes union of genes, inner only keeps genes shared between the dataset 
# can play around with inner/outer joins, concat vs merging 
adata = sc.concat([allData[i] for i in datasetIndices], join = "outer") 
# makes sure each obs has unique label 
adata.obs_names_make_unique()
# convert X to sparse matrix to reduce file size, can also do before cocatination,
adata.X = csr_matrix(adata.X)
# save complete dataset prior to further processing
adata.write_h5ad('allconds_all_genes_postPreProcess_normbeforeCocat.h5ad')

#%% Setup & train model for clustering & further analysis 
# for scvi want number of cells to be at least half the number of genes you have, otherwise might need a diff model
# or only x number most variable genes
adata = sc.read_h5ad('Z:/Dropbox (HMS)/Wilson_Lab_Data/Code/OmicsCode/all_postPreProcess_normbeforeCocat.h5ad')
num_genes = 'all'
adata, model, cluster_dic, markers_scvi = ppf.clusterData(adata, num_genes)

# save model & data
adata.write_h5ad('allconds_allgenes_postCluster_normBeforeCocat.h5ad')
model.save('model_allconds_allgenes_postCluster_normBeforeCocat.model')