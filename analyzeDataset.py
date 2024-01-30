# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 12:02:14 2024

@author: ewest

Run this file to perform basic analyses of a preprocessed and clustered 10xgenomics scRNA seq dataset looking at 
different average marker expression across clusters and performs subclustering of the GABAergic expressing clusters.
Copy & edit to look at different markers. 

clustering performed by scVI model: https://www.nature.com/articles/s41592-018-0229-2 
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

import analysisFunctions as af 
import clusterFunctions as cf
import scanpy as sc 
from scipy.sparse import csr_matrix
    
sc.set_figure_params(dpi=100, dpi_save=100)

#%% Load in preprocessed dataset

# change path if needed
adata = sc.read_h5ad(os.path.join(saveDir, 'allconds_allgenes_postCluster_normBeforeCocat.h5ad'))

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

# plot with cluster number showing
resolution = 1
adata, cluster_dic = cf.defineClusters(adata, resolution)

sc.pl.umap(adata, color = ['leiden'], frameon = False, legend_loc = 'on data', legend_fontsize = 'x-small')

# edit code below to look @ markers of interest
# Looks to have well clustered glia & IPCs
sc.pl.umap(adata, color = ['CG10433','Ilp2','Ilp3'], frameon = False, layer = 'scvi_normalized')
# chol, glut, gaba cells are in 3 diff major clusters but broadly scattered within those
sc.pl.umap(adata, color = ['VAChT','VGlut','Gad1'], frameon = False, layer = 'scvi_normalized')
# kenyon looks to largely be in 3 nearby clusters
sc.pl.umap(adata, color = ['ey','Dop1R2', 'Pka-C1'], frameon = False, layer = 'scvi_normalized')
# monoaminergic neurons also fairly well isolated
sc.pl.umap(adata, color = ['Vmat'], frameon = False, layer = 'scvi_normalized', use_raw = False,vmax = ['p99'])

# uncomment & edit to rename identified clusters
#af.renameClusters(adata, [36,28,47,2,18,25,29], ['s glia','c glia1','c glia2','sp glia1','sp glia2','astrocytes1','astrocytes2'], 1)

#%% Find average marker expression in each cluster 

gene_ids, clusters, obs, average_obs, obs_bool, fraction_obs = af.findAveGeneExpPerCluster(adata)
results_table = af.makeResultsTable(adata)
geneScore_cluster = af.geneScoresPerCluster(results_table, gene_ids, clusters)

#%% Exploratory cluster marker expression plots -- GABA example

# Set thresholds for markers

# assign celltype to cells based on gene exp relative to thresholds
# split dataset by celltype
# recluster 

genesOfInterest = ['Gad1','VAChT', 'Rx']
af.plotGeneScoresAcrossClusters(genesOfInterest,geneScore_cluster)

gene = 'Gad1'

# set threshold based on previous plot
aboveThres = geneScore_cluster.loc['Gad1'][geneScore_cluster.loc['Gad1'] > 0].index.tolist()
adata_gaba = adata[adata.obs['leiden'].isin(aboveThres)]

adata_gaba.X = csr_matrix(adata_gaba.X)
# set number of genes (ranked by variability) to keep & define clusters by
num_genes = 'all'
adata_gaba, model_gaba, cluster_dic_gaba, markers_scvi_gaba = cf.clusterData(adata_gaba, num_genes)

# Make & visualize clusters - GABA
resolution = 2
adata_gaba, cluster_dic = cf.defineClusters(adata_gaba, resolution)

#%% Save new dataset & model
adata_gaba.write_h5ad(os.path.join(saveDir, 'allconds_allgenes_normBeforeCocat_gaba.h5ad'))
model_gaba.save(os.path.join(saveDir, 'model_allconds_allgenes_normBeforeCocat_gaba.model'))

