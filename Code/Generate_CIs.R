##########################################################
# Script to generate CIs around a subset of 100 rhats
##########################################################

rm(list = ls())
set.seed(1) # for reproducibility
library(ggplot2)
library(dplyr)
library(Rcpp)
library(doParallel)
library(doRNG)
source("./simulate_data.R")
sourceCpp("./hmmloglikelihood.cpp")
registerDoParallel(cores = detectCores()-2)
epsilon <- 0.001 # Fix epsilon throughout
CIs <- T
nboot <- 500 
num_CI <- 100 

## Mechanism to generate Ys given fs, distances, k, r, epsilon
simulate_Ys_hmm <- function(frequencies, distances, k, r, epsilon){
  Ys <- simulate_data(frequencies, distances, k = k, r = r, epsilon, rho = 7.4 * 10^(-7))
  return(Ys)
}

## Mechanism to compute MLE given fs, distances, Ys, epsilon
compute_rhat_hmm <- function(frequencies, distances, Ys, epsilon){
  ndata <- nrow(frequencies)
  ll <- function(k, r) loglikelihood_cpp(k, r, Ys, frequencies, distances, epsilon, rho = 7.4 * 10^(-7))
  optimization <- optim(par = c(50, 0.5), fn = function(x) - ll(x[1], x[2]))
  rhat <- optimization$par
  return(rhat)
}

load("../RData/hmmInput_freqs.RData")
dataset_names <- names(hmmInput_freqs)

### Construct bootstrap confidence intervals for 100 select comparisons per data set
for (idataset in 1:length(dataset_names)){
  
  load(sprintf('../RData/%s_mles.RData', dataset_names[idataset]))
  print(dataset_names[idataset])
  ncomp = nrow(mle_df)
  
  # Sample prop to r so as to sample few range
  inds = sample(1:ncomp, size = min(num_CI, ncomp), prob = mle_df$rhat_hmm)     
  dataset_ <- hmmInput_freqs[[dataset_names[idataset]]]
  mle_df_CI <- mle_df

  # sort by chromosome and position
  dataset_ <-  dataset_ %>% arrange(chrom, pos) 
  frequencies <- cbind(1- dataset_$fs, dataset_$fs)
  dataset_$dt <- c(diff( dataset_$pos), Inf)
  pos_change_chrom <- 1 + which(diff( dataset_$chrom) != 0) # find places where chromosome changes
  dataset_$dt[pos_change_chrom-1] <- Inf
  
  for(ind in inds){
    rhats_hmm_boot = as.numeric(foreach(iboot = 1:nboot, .combine = c) %dorng% {
      Ys_boot <- simulate_Ys_hmm(frequencies, distances = dataset_$dt, k = mle_df$khat_hmm[ind], r = mle_df$rhat_hmm[ind], epsilon)
      ll_hmm <- function(k, r) loglikelihood_cpp(k, r, Ys_boot, frequencies, dataset_$dt, epsilon, rho = 7.4 * 10^(-7))
      optimization <- optim(par = c(50, 0.5), fn = function(x) - ll_hmm(x[1], x[2]))
      optimization$par[2]
    })
    
    CI <- as.numeric(quantile(rhats_hmm_boot, probs = c(0.025, 0.975)))
    if(any(is.na(CI))){stop('Bootstrap CI is NA')}
    mle_df_CI[ind, 'CI_low'] = CI[1]
    mle_df_CI[ind, 'CI_high'] = CI[2]
  }
  
  mle_df_CI = mle_df_CI[inds,]
  print(mle_df_CI$rhat_hmm > mle_df_CI$CI_low & mle_df_CI$rhat_hmm < mle_df_CI$CI_high)
  save(mle_df_CI, file = paste0("../RData/", dataset_names[idataset], "_mles_CI.RData"))
}


# # Sample (0,1) range uniformly at random (results in many small)
# inds = sapply(runif(num_CI), function(y, x = mle_df$rhat_hmm){
#   inds_close = which(abs(x-y)==min(abs(x-y)))
#   sample(inds_close, 1) # If more than 1, sample uniformly at random
# })







