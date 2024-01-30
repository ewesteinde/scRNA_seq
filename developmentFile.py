#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 10:18:06 2024

@author: elenawesteinde

This file is used as a testing ground for code development & playing around with different datasets as this repo is still in
development. 

Currently working with the Parker et al thirst dataset & following code @ https://github.com/sims-lab/FlyThirst/blob/main/
"""

#%% load modules
import os
import platform 
if platform.system() == 'Windows':
    os.chdir('Z:\Dropbox (HMS)\Wilson_Lab_Data\Code\OmicsCode')
else:
    os.chdir('/Users/elenawesteinde/Dropbox (HMS)/Wilson_Lab_Data/Code/OmicsCode')

import preprocessFunctions as ppf 
import analysisFunctions as af
import scanpy as sc 
from scipy.sparse import csr_matrix
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sc.set_figure_params(dpi=100, dpi_save=100)
#%% load saved adata & model 

# adata = sc.read_h5ad('/Users/elenawesteinde/Dropbox (HMS)/Wilson_Lab_Data/Code/OmicsCode/all_normBeforeCocat_4000.h5ad')
# adata_preCluster = sc.read_h5ad('Z:/Dropbox (HMS)/Wilson_Lab_Data/Code/OmicsCode/all_postPreProcess_normbeforeCocat.h5ad')
#%% gather data
dataFiles, conds = ppf.gatherFiles(rootDir = 'Z:\Dropbox (HMS)\Wilson_Lab_Data\Omics_datasets\GSE207799_RAW')
uniqueConds = set(conds)
print(uniqueConds)
#%% preprocess data

# play around w/ norm before or after concatenation & ignore or not highly expressed genes

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
# play around with inner/outer joins, concat vs merging 
adata = sc.concat([allData[i] for i in datasetIndices], join = "outer") 
# make sure each obs has unique label 
adata.obs_names_make_unique()
# convert X to sparse matrix to reduce file size, can also do before cocatination,
# really should do to each sample when loading them in to save memory
adata.X = csr_matrix(adata.X)

adata.write_h5ad('allconds_all_genes_postPreProcess_normbeforeCocat.h5ad')

#%% Setup & train model 
# for scvi want number of cells to be at least half the number of genes you have, otherwise might need a diff model
# or only x number most variable genes
adata = sc.read_h5ad('Z:/Dropbox (HMS)/Wilson_Lab_Data/Code/OmicsCode/all_postPreProcess_normbeforeCocat.h5ad')
num_genes = 'all'
adata, model, cluster_dic, markers_scvi = ppf.clusterData(adata, num_genes)

# save model & data
adata.write_h5ad('allconds_allgenes_postCluster_normBeforeCocat.h5ad')
model.save('model_allconds_allgenes_postCluster_normBeforeCocat.model')

# plot with cluster number showing
sc.pl.umap(adata, color = ['leiden'], frameon = False, legend_loc = 'on data', legend_fontsize = 'x-small')


#%% Label broad/known cell types

# neuron markers: nSyb, elav
# cholinergic: VAChT
# glutamatergic: VGlut
# gabaergic: Gad1
# kenyon: ey, Dop1R2, Pka-C1
# monoaminergic: Vmat
# IPCs: Ilp2, Ilp3
# glia: CG10433 
# Ekar found in only photoreceptor neurons neurons

# Looks to have well clustered glia & IPCs
sc.pl.umap(adata, color = ['CG10433','Ilp2','Ilp3'], frameon = False, layer = 'scvi_normalized')
# chol, glut, gaba cells are in 3 diff major clusters but broadly scattered within those
sc.pl.umap(adata, color = ['VAChT','VGlut','Gad1'], frameon = False, layer = 'scvi_normalized')
# kenyon looks to largely be in 3 nearby clusters
sc.pl.umap(adata, color = ['ey','Dop1R2', 'Pka-C1'], frameon = False, layer = 'scvi_normalized')
# monoaminergic neurons also fairly well isolated
sc.pl.umap(adata, color = ['Vmat'], frameon = False, layer = 'scvi_normalized', use_raw = False)


#marker_genes = ['VAChT', 'VGlut', 'Gad1', 'Dop1R2', 'ey', 'Pka-C1', 'Vmat', 'Ilp3', 'Ilp2', 'CG10433']


#%% Map gene IDs to a cluster label

# # D11 wasn't found
# marker_genes = ['VAChT', 'VGlut', 'Gad1', 'Rx', 'Imp', 'Syp', 'pnt', 'svp', 'mamo', 'br','Eip93F', 'ct', 'dan','Optix', 'stg', 'bsh','chinmo']

# #marker_genes = ['Eaat1', 'Gad1', 'Eaat2', 'VGlut','VGAT','Gad1','VAChT','ChAT','Tbh','Tdc2','Tdc1','Trh','SerT', 'Vmat','Ddc','ple','DAT']
# # diff ways of visualizing expression of specific genes across clusters
# ax = sc.pl.dotplot(adata, marker_genes, groupby='leiden',swap_axes = True)
# ax = sc.pl.stacked_violin(adata, marker_genes, groupby='leiden', rotation=90,swap_axes = True)


#%% Exploratory cluster plots
# layers = counts or scvi_normalized

sc.pl.umap(adata, color = ['Gad1'], frameon = False, layer = 'scvi_normalized')#,vmax = ['p99','p99','p99','p99'])
gene_ids, clusters, obs, average_obs, obs_bool, fraction_obs = af.findAveGeneExpPerCluster(adata)
#%% Exploratory ave gene expression plots
gene_ids, clusters, obs, average_obs, obs_bool, fraction_obs = af.findAveGeneExpPerCluster(adata)

# Set thresholds for markers -- initial focus on GABA

# assign celltype to cells based on gene exp relative to thresholds
# split dataset by celltype
# recluster 
# Park paper labelled the EB & FB in the reclustered GABAergic cells --> recreate

results_table = af.makeResultsTable(adata)
geneScore_cluster = af.geneScoresPerCluster(results_table, gene_ids, clusters)

genesOfInterest = ['pnt','Imp', 'Rx']
af.plotGeneScoresAcrossClusters(genesOfInterest,geneScore_cluster)
#%%

gene = 'Gad1'

aboveThres = geneScore_cluster.loc['Gad1'][geneScore_cluster.loc['Gad1'] > 0].index.tolist()

adata_gaba = adata[adata.obs['leiden'].isin(aboveThres)]

adata_gaba.X = csr_matrix(adata_gaba.X)
num_genes = 'all'
adata_gaba, model_gaba, cluster_dic_gaba, markers_scvi_gaba = ppf.clusterData(adata_gaba, num_genes)

#%%
adata_gaba.write_h5ad('allconds_allgenes_normBeforeCocat_gaba.h5ad')
model_gaba.save('model_allconds_allgenes_normBeforeCocat_gaba.model')

#%% Make & visualize clusters - GABA
resolution = 2
adata_gaba, cluster_dic = ppf.defineClusters(adata_gaba, resolution)

#%% Identify & annotate subclusters

# Park paper lists these markers for GABAergic neurons (find citations):
    # Large-field EB Ring neurons: 'cv-c', 'Dh31', 'Octbeta2R', '5-HT7'
    # Small-field EB Ring neurons: 'cv-c', 'Dh31', 'Octbeta2R', NOT 5-HT7
    # Medial FB: 'cv-c', 'Dh31', 'sNPF', 'Octbeta2R'
    # Ventral & dorsal FB: 'cv-c', 'Dh31', 'sNPF', NOT Octbeta2R
    
gene_ids, clusters, obs, average_obs, obs_bool, fraction_obs = af.findAveGeneExpPerCluster(adata_gaba)  
results_table = af.makeResultsTable(adata_gaba)
geneScore_cluster = af.geneScoresPerCluster(results_table, gene_ids, clusters)

# look at single gene
# plt.scatter(x=clusters, y = geneScore_cluster.loc['cv-c'])

genesOfInterest = ['cv-c', 'Dh31', 'Octbeta2R', '5-HT7','chinmo','Rx']
af.plotGeneScoresAcrossClusters(genesOfInterest,geneScore_cluster)
sc.pl.umap(adata_gaba, color = genesOfInterest, frameon = False, layer = 'scvi_normalized',vmax = ['p99','p99','p99','p99','p99','p99'])

# thresholds = [16,0,-5,15,30]

# lfEBR = (set(geneScore_cluster.loc['cv-c'][geneScore_cluster.loc['cv-c'] > 2].index.tolist()) 
#         & set(geneScore_cluster.loc['Dh31'][geneScore_cluster.loc['Dh31'] > 1].index.tolist()) 
#         & set(geneScore_cluster.loc['Octbeta2R'][geneScore_cluster.loc['Octbeta2R'] > -5].index.tolist())
#         & set(geneScore_cluster.loc['5-HT7'][geneScore_cluster.loc['5-HT7'] > 8].index.tolist()))
# smEBR = (set(geneScore_cluster.loc['cv-c'][geneScore_cluster.loc['cv-c'] > 2].index.tolist()) 
#         & set(geneScore_cluster.loc['Dh31'][geneScore_cluster.loc['Dh31'] > 1].index.tolist()) 
#         & set(geneScore_cluster.loc['Octbeta2R'][geneScore_cluster.loc['Octbeta2R'] > -5].index.tolist())
#         & set(geneScore_cluster.loc['5-HT7'][geneScore_cluster.loc['5-HT7'] < 8].index.tolist()))
# mFB = 
# vdRB = 


def testFunction():
    a = 1
    b = 2
    return a, b
    
    
    





