#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 09:16:56 2024

@author: elenawesteinde
"""
def gatherFiles(rootDir):
    #Import modules 
    from os import listdir, system
    from os.path import isfile, join
    
    # Unzip files if necessary
    onlyfiles = [f for f in listdir(rootDir) if isfile(join(rootDir,f)) and ".gz" in f]
    
    # if zipped files exist unzip them 
    if onlyfiles:
        for f in onlyfiles:
            fullpath = join(rootDir,f)
            system("gunzip %s" %(fullpath,))
    
    # Find data files & cocatenate paths into list
    allData = []
    cond = []        
    for f in listdir(rootDir):
        if isfile(join(rootDir,f)) and (".h5" in f or ".tsv" in f):
            allData.append(join(rootDir,f))
        elif not isfile(join(rootDir,f)) and ".DS_Store" not in f:
            for file in listdir(join(rootDir,f)):
                if isfile(join(rootDir,f,file)) and (".h5" in file or ".tsv" in file):
                    allData.append(join(rootDir,f,file))
                    cond.append(f)
            
            
            
    return allData, cond

def doubletRemoval(rawdata,file):
    import scanpy as sc 
    import scvi
    adata = rawdata
    # filter out genes that are found in < 1 cells
    sc.pp.filter_genes(adata, min_cells = 1)
    # keep the top 2000 most variable genes
    sc.pp.highly_variable_genes(adata, n_top_genes = 2000, subset = True, flavor = 'seurat_v3')
    # setup model to detect doublets
    scvi.model.SCVI.setup_anndata(adata)
    vae = scvi.model.SCVI(adata)
    # train
    vae.train()
    solo = scvi.external.SOLO.from_scvi_model(vae)
    solo.train()
    # predict whether cells are likely doublets or not
    df = solo.predict()
    df['prediction'] = solo.predict(soft = False)
    df.index = df.index.map(lambda x: x[:-2])
    # find diff b/w doublet vs singlet likelihood scores
    df['dif'] = df.doublet - df.singlet
    # choose threshold at which to count things as doublets
    doublets = df[(df.prediction == 'doublet') & (df.dif > 1)]
    # load in raw file again
    adata = sc.read_10x_h5(file)
    # remove cells that are likely doublets
    adata.obs['doublet'] = adata.obs.index.isin(doublets.index)
    adata = adata[~adata.obs.doublet]
    
    return adata

def cleanData(adata):
    # Import modules
    # import numpy as np
    import scanpy as sc 

    # Doublet prediction -- look into if appropriate in a dataset where certain neuron pops consist of 2 cells
    # only keep genes that are found in at least 1 cell
    # not sure where the genes come from that aren't in any cell?
    sc.pp.filter_genes(adata,min_cells=1) 
    # calculate & plot percent of mt, ribo RNA & proteins
    adata.var['mt'] = adata.var.index.str.startswith('mt:')
    adata.var['rRNA'] = adata.var.index.str.contains('SrRNA|-rRNA')
    adata.var['rProt'] = adata.var.index.str.contains('RpL|RpS')

    # calculate quality control metrics, 
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','rRNA', 'rProt'], percent_top = None, log1p=False, inplace = True)
    # updates adata.var to contain gene counts & cell dropout metrics (how many cells was the gene present in)
    # adata.obs have various stats related to the mitochondrial & ribosomal genes

    # plot qc metrics
    #sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_rRNA', 'pct_counts_rProt'], jitter = 0.2, multi_panel = True)
    # use qc metrics to remove outliers
    # people set mitochondrial filter anywhere from 5-20%, high percentages indicate dying cells
    # can also set upper bound for num genes per cell, choosing the 99th percentile
    # upper_lim = np.quantile(adata.obs.n_genes_by_counts.values,.99)
    # adata = adata[adata.obs.n_genes_by_counts < upper_lim]


    # Thirst github assumes cells with > 4500 genes are likely doublets while cells with < 300 genes carry insufficient info. Removes those outside that range
    # max allowed UMIs is 25000 --> also assumes above this is likely doublets
    # removes cells with at least one of each metric:
    # < 300 features (genes)
    # > 4500 features
    # > 25 000 UMIs
    # > 15% mt-RNA
    # > 10% rRNA
    # > 15% ribosomal proteins 

    adata = adata[adata.obs.pct_counts_mt < 20]
    adata = adata[adata.obs.pct_counts_rRNA < 10]
    adata = adata[adata.obs.pct_counts_rProt < 15]
    adata = adata[(adata.obs.n_genes_by_counts > 300) | (adata.obs.n_genes_by_counts < 4500)]
    adata = adata[adata.obs.total_counts < 25000]

    #sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_rRNA', 'pct_counts_rProt'], jitter = 0.2, multi_panel = True)

    return adata

def preprocess(file):
    import scanpy as sc 
    import platform
    # load file & clean data
    rawdata = sc.read_10x_h5(file)
    #adata.obs: cell barcodes, no additional obs yet
    #adata.var: shows genes, no additional info yet
    # adata.X: numpy array of counts   
    doubletRemovedData = doubletRemoval(rawdata, file)
    adata = cleanData(doubletRemovedData) 
    if platform.system() == 'Windows':
        f = file.split('\\')[6]
    else:
        f = file.split('/')[8]
    adata.obs['Sample'] = f.split('_')[0]
    # save raw data prior to normalization 
    adata.layers['counts'] = adata.X.copy()
    # Normalize data
    sc.pp.normalize_total(adata,target_sum=1e4) 
    sc.pp.log1p(adata)
    # freeze data here before doing further modifications, keep 'raw' version
    adata.raw = adata
    # scores expression level of roX1 (X chromosome associated gene) in each cell to regress against later
    sc.tl.score_genes(adata,['lncRNA:roX1'], score_name = 'sex',use_raw = False)
    
    return adata 

def defineClusters(adata, res):
    import scanpy as sc
    # Make & visualize clusters

    sc.tl.leiden(adata, resolution=res) 
    # plot clusters, should see decent overlap between the diff samples, indicates
    # proper integration of samples --> no major batch effect remaining
    sc.pl.umap(adata, color = ['leiden', 'Sample'], frameon = False)

    cluster_nums = list(range(len(adata.obs['leiden'].unique())))
    cluster_names = [str(x) for x in cluster_nums]
    cluster_dic = dict(zip(cluster_nums, cluster_names))
    adata.uns['cluster_dic'] = cluster_dic
      
    # plot with cluster number showing
    sc.pl.umap(adata, color = ['leiden'], frameon = False, legend_loc = 'on data', legend_fontsize = 'x-small')
    
    return adata, cluster_dic

def identifyClusterMarkers(adata, model, find_markers):
    import scanpy as sc
    # get array of gene markers in clusters
    #markers = sc.get.rank_genes_groups_df(adata, None)
    # filter out markers with pvalues above 0.05 & logfold change of at least 0.5
    # logfoldchange represents degree to which genes are up or down regulated 
    #markers = markers[(markers.pvals_adj < 0.05) & (markers.logfoldchanges > 0.5)]
    # use model
    # find markers for clusters
    # updates a new layer
    sc.tl.rank_genes_groups(adata,'leiden',method='wilcoxon')
    # plots gene exp in clusters
    sc.pl.rank_genes_groups(adata,n_genes=20,sharey=False)
    
    if find_markers:        
        markers = sc.get.rank_genes_groups_df(adata, None)
        markers = markers[(markers.pvals_adj < 0.05) & (markers.logfoldchanges > 0.5)]
        
        markers_scvi = model.differential_expression(groupby = 'leiden')
        markers_scvi = markers_scvi[(markers_scvi['is_de_fdr_0.05']) & (markers_scvi.lfc_mean > 0.5)]
        
        marker_dict = {'markers': markers, 'markers_scvi': markers_scvi}
    else:
        marker_dict = {}
    
    return adata, marker_dict
    

def clusterData(adata, num_genes):
    # Setup & train model 
    # for scvi want number of cells to be at least half the number of genes you have, otherwise might need a diff model
    # or only x number most variable genes    
    import scanpy as sc
    import scvi
    if num_genes == 'all':
        num_genes = adata.shape[1]
    elif num_genes > adata.shape[1]:
        num_genes = adata.shape[1]
    # filter out genes that occur less than 3 times, causes error if too many counts are near 0
    sc.pp.filter_genes(adata, min_counts = 8)
    sc.pp.highly_variable_genes(adata, n_top_genes = num_genes, subset= True, layer = 'counts', flavor = 'seurat_v3', 
                                batch_key = 'Sample') # no batch key if one sample                           
    sc.pl.highly_variable_genes(adata)

    # then setup model
    scvi.model.SCVI.setup_anndata(adata, layer = 'counts',
                                  categorical_covariate_keys = ['Sample'],
                                  continuous_covariate_keys = ['pct_counts_mt', 'total_counts', 'pct_counts_rRNA', 'pct_counts_rProt'])
    model = scvi.model.SCVI(adata)
    model.train()

    # save some model outputs
    adata.obsm['X_scVI'] = model.get_latent_representation()
    adata.layers['scvi_normalized'] = model.get_normalized_expression(library_size = 1e4)

    sc.pp.neighbors(adata, use_rep = 'X_scVI') # use latent representation from model to calculate neighbours
    sc.tl.umap(adata)
    
    adata, cluster_dic = defineClusters(adata, 1) 
    adata, markers_dict = identifyClusterMarkers(adata, model, 0)
    
    
    
    return adata, model, cluster_dic, markers_dict

