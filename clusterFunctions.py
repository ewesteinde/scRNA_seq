# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 13:23:37 2024

@author: ewest

Functions to perform clustering of a 10xgenomics scRNA seq dataset

"""

def defineClusters(adata, res):
    import scanpy as sc
    # Make & visualize clusters

    sc.tl.leiden(adata, resolution=res) 
    # plot clusters, should see decent overlap between the diff samples, indicates
    # proper integration of samples --> no major batch effect remaining
    sc.pl.umap(adata, color = ['leiden', 'Sample'], frameon = False)

    cluster_nums = list(range(len(adata.obs['leiden'].unique())))
    cluster_names = [str(x) for x in cluster_nums]
    cluster_dic = dict(zip(cluster_names, cluster_names))
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
    
    sc.pl.umap(adata, color = ['leiden', 'Sample'], frameon = False)
        
    return adata, model, cluster_dic, markers_dict
