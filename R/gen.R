
library(edgeR)
library(mvtnorm)
library(Matrix)
library(matrixcalc)
library(corpcor)
library(glasso)
library(cvTools)
library(huge)
library(edgeR)
 
####################################################################
#' A function to  
#'
#' @param nrow 
#' @param sparse
#' @param min 
#' @param max
#' @return M
#' @export
generate_random_pd_sym_matrix<-function(nrow,sparse,min,max){
  if(sparse<0 || sparse>1 || nrow<=0){
    stop ("Error: Sparse needs to be between 0 and 1, and nrow needs to be positive integer")
  }
  M=matrix(0,nrow=nrow,ncol=nrow)
  x=runif(nrow*(nrow-1)/2, min=min,max=max)
  y=sample.int(nrow*(nrow-1)/2, size=sparse*nrow*(nrow-1)/2)
  x[y]=0
  k=1
  for(i in 1:(nrow-1)){
    for(j in (i+1):nrow){
      if(abs(x[k])>1e-2){
        M[i,j]=x[k]
      }else{
        M[i,j]=0
      }
      k=k+1
    }
  }
  M=M+t(M)
  M=M+diag(runif(nrow,min=0.5,max=1.5),nrow=nrow)
  cond=rank.condition(M,tol=1e-5)$condition
  while(cond>30 || !is.positive.definite(M,tol=1e-5)){
    M=M+diag(runif(nrow,min=0.05,max=0.5),nrow=nrow)
    cond=rank.condition(M,tol=1e-5)$condition
  }
  return(M)
}
 
#########################################################
##########################################################
#' A function to  
#' 
#' @param  X_samp
#' @return list(TMM=X_TMM)
#' @export
normalized_data <- function(X_samp){
  N = nrow(X_samp)
  t <-as.vector(calcNormFactors(as.matrix(t(X_samp)), method = "TMM"))
  X_TMM = X_samp
  for(i in 1:N){
    X_TMM[i,] = X_samp[i,]/t[i]
  }
 
  return(list(TMM=X_TMM))
}
##########################################################
#' A function to    
#' @param y 
#' @param MinValue
#' @param MaxValue 
#' @return x
#' @export
sampling <-  function(y, MinValue , MaxValue){
  N = nrow(y)
  d = ncol(y)
  x=matrix(NA,ncol=d,nrow=N) 
  for (i in 1:N){
    x[i,] = rmultinom(1, size = runif(1,MinValue , MaxValue), prob = y[i,])
  }
	 
  return(x)
}
##########################################################


#' A function to generate the synthetic data
#'
#' @param K number of componenet in the synthetic data
#' @param N number of samples in the synthetic sample-taxa count matrix
#' @param d number of taxa in the synthetic sample-taxa count matrix
#' @param sp sparsity level. value between 0 and 1. 
#' @param type "orig" original count data, "samp" sampled data , "TMM" after TMM normalization
#' @return
#' @export
generate <- function(K , N , d , graph_type){
     
real_precision = list()


library(devtools)
library(SpiecEasi)
data(amgut1.filt)
depths <- rowSums(amgut1.filt)
amgut1.filt.n  <- t(apply(amgut1.filt, 1, norm_to_total))
amgut1.filt.cs <- round(amgut1.filt.n * min(depths))

d <- ncol(amgut1.filt.cs)
n <- nrow(amgut1.filt.cs)
e <- d

 

set.seed(10010)
if( graph_type== 'band') { graph <- make_graph('band', d, e) }
if( graph_type== 'cluster') { graph <- make_graph('cluster', d, e) }
if( graph_type== 'scale_free') { graph <- make_graph('scale_free', d, e) }
 

graph_new = graph
#Prec = generate_random_pd_sym_matrix(nrow=d,sparse=sparse,min=-1,max=1)
Prec  <- graph2prec(graph)
Cor   <- cov2cor(prec2cov(Prec))
inv1 = Prec
X <- synth_comm_from_counts(amgut1.filt.cs, mar=2, distr='zinegbin', Sigma=Cor, n=10*d)
X1 <- sampling( X , 5000 , 10000)
 

 
new = sample(c(1:d), d , replace = FALSE, prob = NULL)
for(i in 1:d){for(j in 1:d){ graph_new[i,j] = graph[new[i],new[j]]}}
#print(myfr(graph , graph_new))
#Prec = generate_random_pd_sym_matrix(nrow=d,sparse=sparse,min=-1,max=1)
Prec  <- graph2prec(graph_new)
Cor   <- cov2cor(prec2cov(Prec))
inv2 = Prec
X <- synth_comm_from_counts(amgut1.filt.cs, mar=2, distr='zinegbin', Sigma=Cor, n=10*d)
X2 <- sampling( X , 5000 , 10000)
 

 
new = sample(c(1:d), d , replace = FALSE, prob = NULL)
for(i in 1:d){for(j in 1:d){ graph_new[i,j] = graph[new[i],new[j]]}}
#print(myfr(graph , graph_new))
#Prec = generate_random_pd_sym_matrix(nrow=d,sparse=sparse,min=-1,max=1)
Prec  <- graph2prec(graph_new)
Cor   <- cov2cor(prec2cov(Prec))
inv3 = Prec
X <- synth_comm_from_counts(amgut1.filt.cs, mar=2, distr='zinegbin', Sigma=Cor, n=10*d)
X3 <- sampling( X , 5000 , 10000)
 

if(K==1){X = X1[1:N,]
         real_precision[[1]] = inv1}
if(K==2){X = rbind(X1[1:N,] , X2[1:N,])
         real_precision[[1]] = inv1
         real_precision[[2]] = inv2}
if(K==3){X = rbind(X1[1:N,] , X2[1:N,] , X3[1:N,])
         real_precision[[1]] = inv1
         real_precision[[2]] = inv2
         real_precision[[3]] = inv3}	 
 
	 
 
  return( list("M"=X, "real_precision"=real_precision))
}
	
	 
     
         
	
 	 

 

 
 
 
 
 

