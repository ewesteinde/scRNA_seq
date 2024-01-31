# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 11:47:26 2024

@author: ewest

Code to compare cluster similarity across datasets
Iterative clustering to identify top markers per genes: https://github.com/felixhorns/FlyPN

"""

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
    
    
import scanpy as sc
import clusterFunctions as cf
import preprocessFunctions as ppf
import scvi
import pandas as pd


#%% find marker genes for each cluster if not already done

# have to change to dropbox folder before loading in
if platform.system() == 'Windows':
    os.chdir('Z:\Dropbox (HMS)\Wilson_Lab_Data\Omics_datasets\processedData')
    adata = sc.read_h5ad('allconds_all_genes_postCluster_normbeforeCocat.h5ad')
    model = scvi.model.SCVI.load('model_allconds_allgenes_postCluster_normBeforeCocat.model', adata)
    os.chdir('C:\Code\scRNA_seq')
else:
    os.chdir('/Users/elenawesteinde/Documents/Omics_proj/scRNA_seq')

adata, cluster_dic = cf.defineClusters(adata, 1.5)
adata, marker_dic = cf.identifyClusterMarkers(adata, model, 1)
markers_scvi = marker_dic['markers_scvi']
markers = marker_dic['markers']

adata.write_h5ad(os.path.join(saveDir, 'allconds_allgenes_clusterRes_1.5.h5ad'))
markers.to_csv(os.path.join(saveDir, 'Markers_allconds_allgenes_clusterRes_1.5.csv'))
markers_scvi.Dataframe.to_csv(os.path.join(saveDir, 'Markers_scvi_allconds_allgenes_clusterRes_1.5.csv'))

#%%  load in marker dataframes
ageBrain_markersPerCluster = pd.read_excel('Z:\Dropbox (HMS)\Wilson_Lab_Data\Omics_datasets\datasets\AgingBrain\MarkerGenesperCluster.xlsx')
markers = pd.read_csv(os.path.join(saveDir, 'Markers_allconds_allgenes_clusterRes_1.5.csv'))
# ageBrain_markersPerCluster.obs: gene names
# ageBrain_markersPerCluster.var: p-Value, Average Difference, Cluster ID

#%% clean up markers to contain only the same genes as in ageBrain_markersPerCluster
markers2 = markers.copy()
markers2 = markers2.drop(labels = ['pvals_adj','scores'], axis = 1)
# rename so column labels match (CHECK AVERAGE DIFF IS ACTUALLY LFC)
ageBrain = ageBrain_markersPerCluster.rename(columns ={'Gene Symbol':'names','p-Value':'pvals','Average Difference':'logfoldchanges','Cluster ID':'group'})

# only keep rows with genes shared between the two 
markers2 = markers2.loc[markers2.names.isin(ageBrain.names)]
ageBrain = ageBrain.loc[ageBrain.names.isin(markers2.names)]
# diff num of clusters b/w the two, and diff num & identity of marker genes for each cluster

#%% start by looking for best match to cluster 61 (dFB)
ageBrain61 = ageBrain.loc[ageBrain.group == 61]


#%% pivot to compare entire datasets to eachother, aging brain = reference dataset
rootDir = 'Z:\Dropbox (HMS)\Wilson_Lab_Data\Omics_datasets\datasets\AgingBrain\wholeDataset'
allData, cond = ppf.gatherFiles(rootDir)

## SCOPE loom files are non-standard can't just be loaded into python
import tarfile 
file = tarfile.open('Z:/Dropbox (HMS)/Wilson_Lab_Data/Omics_datasets/datasets/AgingBrain/wholeDataset/GSE107451_DGRP-551_w1118_WholeBrain_57k_0d_1d_3d_6d_9d_15d_30d_50d_10X_DGEM_MEX.mtx.tsv.tar.gz')
file.extractall('Z:\Dropbox (HMS)\Wilson_Lab_Data\Omics_datasets\datasets\AgingBrain\wholeDataset')
file.close()

test = sc.read_10x_mtx('Z:\Dropbox (HMS)\Wilson_Lab_Data\Omics_datasets\datasets\AgingBrain\wholeDataset')



