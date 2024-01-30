# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 11:56:21 2024

@author: ewest

Run this file to preprocess a 10xgenomics scRNA seq dataset with one or more batches/conditions.
Performs doublet prediction & removal, removal of cells with high mitochondrial & ribosomal RNA expression, 
and normalization of data.
Following preprocessing will perfrom clustering on the data and save the result. 
Doublet prediction and clustering performed by scVI model: https://www.nature.com/articles/s41592-018-0229-2 

"""
#%% set code, data directories & figure settings
import os
import platform 

# replace paths to code directory if needed
if platform.system() == 'Windows':
    os.chdir('C:\Code\scRNA_seq')
    saveDir = 'Z:\Dropbox (HMS)\Wilson_Lab_Data\Omics_datasets\processedData'
    datasetDir = 'Z:\Dropbox (HMS)\Wilson_Lab_Data\Omics_datasets\datasets\GSE207799_RAW'
else:
    os.chdir('/Users/elenawesteinde/Documents/Omics_proj/scRNA_seq')
    saveDir = '/Users/elenawesteinde/Dropbox (HMS)/Wilson_Lab_Data/Omics_datasets/processedData'
    datasetDir = '/Users/elenawesteinde/Dropbox (HMS)/Wilson_Lab_Data/Omics_datasets/datasets/GSE207799_RAW'
    
#%% load modules

import preprocessFunctions as ppf 
import clusterFunctions as cf
import scanpy as sc 
from scipy.sparse import csr_matrix
    
sc.set_figure_params(dpi=100, dpi_save=100)

#%% gather data
dataFiles, conds = ppf.gatherFiles(datasetDir)
uniqueConds = set(conds)
print(uniqueConds)
#%% preprocess data

# preprocesses each individual dataset (predicts doublets, cleans & normalizes data) seperately and returns a list with all of them
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
adata.write_h5ad(os.path.join(saveDir, 'allconds_all_genes_postPreProcess_normbeforeCocat.h5ad'))

#%% Setup & train model for clustering & further analysis 
# for scvi want number of cells to be at least half the number of genes you have, otherwise might need a diff model
# or only x number most variable genes
adata = sc.read_h5ad('Z:/Dropbox (HMS)/Wilson_Lab_Data/Code/OmicsCode/all_postPreProcess_normbeforeCocat.h5ad')
# set number of genes (ranked by variability) to keep & define clusters by
num_genes = 'all'
adata, model, _, _ = cf.clusterData(adata, num_genes)

# save model & data
adata.write_h5ad(os.path.join(saveDir, 'allconds_allgenes_postCluster_normBeforeCocat.h5ad'))
model.save(os.path.join(saveDir, 'model_allconds_allgenes_postCluster_normBeforeCocat.model'))