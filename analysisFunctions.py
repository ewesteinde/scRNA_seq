#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 16:56:15 2024

@author: elenawesteinde
"""

def findAveGeneExpPerCluster(adata):
    # Import modules
    import pandas as pd
    
    gene_ids = adata.raw.var.index.values
    clusters = adata.obs['leiden'].cat.categories
    obs = adata.raw[:,gene_ids].X.toarray()
    obs = pd.DataFrame(obs,columns = gene_ids,index=adata.obs['leiden'])
    average_obs = obs.groupby(level=0).mean()
    obs_bool = obs.astype(bool) # any exp vs none
    fraction_obs = obs_bool.groupby(level=0).sum()/obs_bool.groupby(level=0).count()
    average_obs.T.to_csv("average.csv")
    fraction_obs.T.to_csv("fraction.csv")
    
    return gene_ids, clusters, obs, average_obs, obs_bool, fraction_obs

def makeResultsTable(adata):
    import pandas as pd
    results = adata.uns['rank_genes_groups']
    groups = results['names'].dtype.names
    results_table = pd.DataFrame({group + '_' + key[:1]: results[key][group]
        for group in groups for key in ['names', 'scores']})
    
    return results_table


def geneScoresPerCluster(results_table, gene_ids, clusters):
    import numpy as np
    import pandas as pd
    genes_ordered = pd.DataFrame()
    scores_ordered = pd.DataFrame()
    for col in results_table.columns:
        group = col[0:-2]
        if 'n' in col:       
            genes_ordered[group] = results_table[col]
        if 's' in col:
            scores_ordered[group] = results_table[col]
            
    geneScore_cluster = pd.DataFrame(index = np.sort(gene_ids), columns = clusters)
    
    idx = list(range(len(gene_ids)))
    for col in genes_ordered:
        cluster_genes = np.array(genes_ordered[col])
        sort_index = np.argsort(cluster_genes)
        geneScore_cluster[col].iloc[idx] = scores_ordered[col][sort_index]
        
    return geneScore_cluster

def renameClusters(adata, toRename, newNames, toPlot):
    import scanpy as sc
    numClusters = len(adata.obs['leiden'].unique())
    if max(toRename) <= numClusters:
        cluster_dic = adata.uns['cluster_dic']
        cluster_dic.update({key:newNames[idx] for idx, key in enumerate(toRename)})        
        new_cluster_names = list(cluster_dic.values())
        adata.rename_categories('leiden', new_cluster_names)
        adata.uns['cluster_dic'] = cluster_dic
        sc.tl.rank_genes_groups(adata,'leiden',method='wilcoxon')
    else:
        print('Tried to rename a nonexistant cluster')
        
    if toPlot:
        import scanpy as sc
        sc.pl.umap(adata, color = ['leiden'], frameon = False, legend_loc = 'on data', legend_fontsize = 'x-small')
        
def plotGeneScoresAcrossClusters(genesOfInterest,geneScore_cluster):
    import math
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    numYticks = 10
    numSubplots = len(genesOfInterest)
    numCol = 4
    if numCol > numSubplots:
        numCol = numSubplots
    numRow = math.ceil(numSubplots/numCol)
    fig,ax = plt.subplots(numRow,numCol, figsize = (10,5))
    plt.tight_layout()
    fig.set_figwidth(10)
    
    
    count = 0   
    for row in range(numRow):
        for col in range(numCol):
            if count < numSubplots:
                gene = genesOfInterest[count]
                geneScores = geneScore_cluster.loc[gene]
                clusterNames = geneScore_cluster.columns
                orderedIdx = geneScores.argsort()
                orderedScores = geneScores[orderedIdx]
                orderedClusters = clusterNames[orderedIdx]
                if numRow == 1:
                    ax[col].scatter(x = orderedClusters, y = orderedScores, s = 3)
                    ax[col].set_title(gene)
                    ax[col].set_xticks([])               
                    yticks = mpl.ticker.MaxNLocator(numYticks)
                    ax[col].yaxis.set_major_locator(yticks)
                else:
                    ax[row,col].scatter(x = orderedClusters, y = orderedScores, s = 3)
                    ax[row,col].set_title(gene)
                    ax[row,col].set_xticks([])
                    yticks = mpl.ticker.MaxNLocator(numYticks)
                    ax[row,col].yaxis.set_major_locator(yticks)
            else:
                fig.delaxes(ax[row,col])                
            count += 1
    


    