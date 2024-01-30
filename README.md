# scRNA_seq
## Code used to process and analyze scRNA_seq data

Code to analyse *Drosophila melanogaster* single cell RNA sequencing datasets to identify specific cell populations from their expression profiles.  

In its current form this code generally follows along the early analytical steps described in the Park., *et al* paper https://pubmed.ncbi.nlm.nih.gov/35963239/ and their associated github https://github.com/sims-lab/FlyThirst/tree/main, with various modifications. It is a work in progress and being continuously updated. 

Three major steps are described:
1. preprocessing
2. Clustering
3. Marker expression analysis and population identification (in development)

### Preprocessing (doublet Removal & cleaning)

Mistakes can be made when isolating cells for single cell RNA sequencing resulting in 'doublets' where 2 (or more) cells are accidently counted as one. Identifying doublets is difficult but various algorithms exist to estimate which cells may be doublets and one can also simply set a threshold of # genes per cell or # Unique Molecular Identifiers (UMIs) where any cell above the threshold is assumed to be a doublet. These cells are then removed from the datasets prior to further processing. 

High mitochondrial RNA or ribosomal RNA counts are signs of an unhealthy cell and should also generally be removed from the dataset, The exact threshold is set depending on the dataset and experimental goals. 

Doublet prediction and clustering performed by scVI model: https://www.nature.com/articles/s41592-018-0229-2 

##### Basic QC metric pre cleaning:
![pre cleaning](https://github.com/ewesteinde/scRNA_seq/blob/main/exampleFigures/preCleaning.png)
##### Basic QC metric post cleaning:
![post cleaning](https://github.com/ewesteinde/scRNA_seq/blob/main/exampleFigures/postCleaning.png)

### Clustering

following these cleaning steps, each cell is normalized by total counts over all genes, such that every cell as the same total gene count after normalization and individual datasets are concatinated. Now clustering can occur to group cells more similar to eachother together. 
![cluster example]([Cluster pic](https://github.com/ewesteinde/scRNA_seq/blob/main/exampleFigures/allData_clusters.png))

Following concatinated (or alternatively, merging which only keeps genes found across every dataset) and clustering one should check for high overlap across the datasets (samples). If each dataset forms highly seperated clusters this indicates significant batch effects are affecting the data. 
![Sample merge example](https://github.com/ewesteinde/scRNA_seq/blob/main/exampleFigures/allData_samples.png)

clustering performed by scVI model: https://www.nature.com/articles/s41592-018-0229-2 

### Population Identification 

We can then check for expression of markers of interest (and that they look to match expression patterns seen in the original Park *et al* analysis)

##### cholinergic: VAChT, glutamatergic: VGlut, gabaergic: Gad1
![VAChT, VGlut, Gad1](https://github.com/ewesteinde/scRNA_seq/blob/main/exampleFigures/ACh_Glut_GABA.png)
##### kenyon: ey, Dop1R2, Pka-C1
![ey, Dop1R2, Pka-C1](https://github.com/ewesteinde/scRNA_seq/blob/main/exampleFigures/KenyonCells.png)
##### glia: CG10433,  IPCs: Ilp2, Ilp3
![CG10433, Ilp2, Ilp3](https://github.com/ewesteinde/scRNA_seq/blob/main/exampleFigures/glia_IPCs.png)

This pipeline is currently being compared, validated against, and expanded to other available datasets and analyses. 
