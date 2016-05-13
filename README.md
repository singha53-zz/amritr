for personal use only

# Table of Contents
1. [Batch correction](#batch-correction)
2. [Classification algorithms](#classification-algorithms)
3. [Network analyses](#network-analyses)
4. [Deconvolution Methods](#deconvolution-methods)
5. [Exploratory Data Analysis](#exploratory-data-analysis)
6. [Generalized Linear models](#generalized-linear-models)
7. [Time series analyses](#time-series-analyses)
8. [Linear algebra concepts](#linear-algebra-concepts)


## Batch correction
### ComBat
    Johnson WE, Li C, Rabinovic A (2007) Adjusting batch effects in microarray expression data using empirical Bayes methods. Biostatistics, 8 (1), 118-127

What is Combat?
Combat is an empiral bayes (EB) method used to adjust for batch effects (especially when the dataset is small)
batch effect: "systematic 'batch effects' or non-biological differences that make samples in different batches not directly comparable"
e.g. time, reagents used, chips --> e.g. in meta-analyses where samples from different studies are pooled

![](/Users/asingh/Documents/languages/R/amritr/figs/batchAdjustment.tiff)

  * Method 1: dChip software standardizes genes in each batch separately
  * Method 2: ComBat apply an EB method (useful for small sample sizes)
  * Method 3: using SVD to adjust for batch requires large small sizes

  * Combat uses the systematic biases common across all genes to make adjustments. This method assumes that the batch phenomena affects genes in similar ways (e.g. increased expression and variability)
    * estimate the mean (location) and variance (scale) of the batch effect by pooling information across genes in each batch and shrinking these batch estimates towards the overal mean of the batch effect estimates
    * the resulting estimates of the batch effect (mean and variance) are used to adjust the data for batch effects
    
**ComBat procedure**
Note: data is normalized and absent and noisy features have been filtered out

The resulting dataset contains m batches with $n_i$ samples within batch *i* for *i* = 1, ..., m for gene *g* = 1, ..., G. Assume the model:
  
  $Y_i_j_g$ = $/alpha$_g

```r





library(sva)



```


### Surrogate Variable Analysis (sva)
    Leek JT and Storey JD. (2007) Capturing heterogeneity in gene expression studies by ‘Surrogate Variable Analysis’. PLoS Genetics, 3: e161.

### Frozen Surrogate Variable Analysis (fsva)
    Parker HS, Bravo HC, Leek JT (2013) Removing batch effects for prediction problems with frozen surrogate variable analysis arXiv:1301.3947

## Classification algorithms
### Enet
## meta-analyses
### Random Forest
### sPLS-DA
### sPLS
#### Correlation circle
### OPLS-DA
### Deep Learning (H2o)

## Network analyses
### DINGO
### PANDA
### SNF

## Deconvolution Methods
### csSAM
### cellCODE
### CIBERSORT

## Exploratory Data Analysis
### Principal Component Analysis (PCA)
### Non-negative Matrix Factorization (NMF)
### Multi-dimensional scaling
### Clustering Algorithms
#### Hierarchical clustering
#### k-means

## Generalized Linear models
### linear regression
### robust linear regression
### weighted linear regression

## Time series analyses
## Mixed-effects models
## Linear mixed-effects model splines (lmms)

## Linear algebra concepts
### Singular Value Decomposition (SVD)
### covariance matrix
### robust covariable matrix
