# How to use the functions in this package
Amrit Singh  
December 3, 2016  



## biomarker pipeline (compare Basal vs. Her2)


```r
library(amritr); data("pathways")
library(mixOmics); data("breast.TCGA")
```

```
## Loading required package: MASS
```

```
## Loading required package: lattice
```

```
## Loading required package: ggplot2
```

```
## Warning: package 'ggplot2' was built under R version 3.3.2
```

```
## 
## Loaded mixOmics 6.1.1
## 
## Visit http://www.mixOmics.org for more details about our methods.
## Any bug reports or comments? Notify us at mixomics at math.univ-toulouse.fr or https://bitbucket.org/klecao/package-mixomics/issues
## 
## Thank you for using mixOmics!
```

```r
Y = breast.TCGA$data.train$subtype
X = breast.TCGA$data.train$mrna[which(Y != "LumA")]
Y = droplevels(Y[which(Y != "LumA")])

## run biomarker pipeline for a binary response
#biomarkerPipeline = function(X = X, Y = Y, topranked = 50, validation = "Mfold", M = 5, iter = 1, threads = 1, progressBar = TRUE, pathways = pathways)
```
