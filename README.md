# MicNet
```
library(mvtnorm)
library(Matrix)
library(matrixcalc)
library(corpcor)
library(cvTools)
library(huge)
library(edgeR)
library(rstan)
library(sfsmisc)
library(mclust)
library(coda)
library(cvTools)
library(cluster)
library(factoextra)
library(igraph)
library(ROCR)

library("devtools")
install_github("sahatava/MicNet")
library("MicNet")
```
# 1 - generate the dataset
## usage

```
out_generate = generate(K , N , d , graph_type)
```
## arguments
```
K-------> number of components
N-------> number of samples in each component
d-------> number of taxa for the synthetic sample-taxa count matrix
graph_type------>  "band", "cluster" or "scale_free"
```
## value
```
out_generate$M -----------------> synthetic sample-taxa count matrix(rows represent samples and columns represent taxa)
out_generate$real_precision[[i]]-----> i_th real precision matrix
```


# 2 - applying the MixGGM
## usage

```
out_MixGGM = MixGGM(M , K , penalty , init , rep )
```
## arguments
```
M--------> name of the csv file which includes the sample-taxa matrix
K-----------> number of components 
penalty------> method to select tuning parameter: "no_sparsity" , "CV" (cross validation) , "StARS" , "fixed" , "iterative" 
init---------> initialization can be "Kmeans" or "Random"
rep----------> number of repeats when initialization is "Random"
out---------> the type of output for each componenet: precision matrix, partial correlation, lambda, mixing coefficiant, clustering membership 
threshold---> a value between 0 and 1. This threshold is used to generate adjacency matrix from partial correlation
```
## value
```
out_MixGGM$precision[[i]]-----> i_th precision matrix(dXd)
out_MixGGM$partial[[i]]-------> i_th partial correlation matrix(dXd)
out_MixGGM$lambda[[i]]--------> i_th lambda matrix(NXd)
out_MixGGM$pi-----------------> mixing coefficiant
out_MixGGM$cluster------------> clustering membership
```


# 2 - applying the MixMCMC
## usage

```
out_MixMCMC = MixMPLN(M , K , penalty , init , rep )
```
## arguments
```
M--------> name of the csv file which includes the sample-taxa matrix
K-----------> number of components 
penalty------> method to select tuning parameter: "no_sparsity" , "CV" (cross validation) , "StARS" , "fixed" , "iterative" 
init---------> initialization can be "Kmeans" or "Random"
rep----------> number of repeats when initialization is "Random"
out---------> the type of output for each componenet: precision matrix, partial correlation, lambda, mixing coefficiant, clustering membership 
threshold---> a value between 0 and 1. This threshold is used to generate adjacency matrix from partial correlation

```
## value
```
out_MixMCMC$precision[[i]]-----> i_th precision matrix(dXd)
out_MixMCMC$partial[[i]]-------> i_th partial correlation matrix(dXd)
out_MixMCMC$lambda[[i]]--------> i_th lambda matrix(NXd)
out_MixMCMC$pi-----------------> mixing coefficiant
out_MixMCMC$cluster------------> clustering membership
```

# 3 - applying the MixMPLN
## usage

```
out_MixMPLN = MixMPLN(M , K , penalty , init , rep )
```
## arguments
```
M--------> name of the csv file which includes the sample-taxa matrix
K-----------> number of components 
penalty------> method to select tuning parameter: "no_sparsity" , "CV" (cross validation) , "StARS" , "fixed" , "iterative" 
init---------> initialization can be "Kmeans" or "Random"
rep----------> number of repeats when initialization is "Random"
out---------> the type of output for each componenet: precision matrix, partial correlation, lambda, mixing coefficiant, clustering membership 
threshold---> a value between 0 and 1. This threshold is used to generate adjacency matrix from partial correlation
```
## value
```
out_MixMPLN$precision[[i]]-----> i_th precision matrix(dXd)
out_MixMPLN$partial[[i]]-------> i_th partial correlation matrix(dXd)
out_MixMPLN$lambda[[i]]--------> i_th lambda matrix(NXd)
out_MixMPLN$pi-----------------> mixing coefficiant
out_MixMPLN$cluster------------> clustering membership```
