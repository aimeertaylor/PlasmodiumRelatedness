#######################################################################
# Script to simulate data for the plot that shows the sensitivity of IBS
# Need to return hs so that they can be added to the sim plot
#######################################################################
rm(list = ls())
library(ggplot2)
library(dplyr)
library(Rcpp)
library(doParallel)
library(doRNG)
source("./simulate_data.R")
sourceCpp("./hmmloglikelihood.cpp")
registerDoParallel(cores = detectCores()-2)
set.seed(1) # for reproducibility
MinSampleSize = 100 # Remove Gambia (too few markers)
load("../RData/hmmInput_freqs.RData") # Load all available data 
load("../RData/nsnps_nsamps.RData")
sites = names(which(sapply(hmmInput_freqs, ncol)>MinSampleSize))
efxd = 0.001 # Fix epsilon throughout
rfxd = 0.5 # Data generating relatedness
kfxd = 8 # The number of generations
nsim = 1000 # Number of sample pairs to simulate (5000 produces smooth but need to run over night)
kinit = kfxd
rinit = rfxd

## Mechanism to compute MLE given fs, distances, Ys, epsilon
compute_rhat_hmm <- function(frequencies, distances, Ys, epsilon){
  ndata <- nrow(frequencies)
  ll <- function(k, r) loglikelihood_cpp(k, r, Ys, frequencies, distances, epsilon, rho = 7.4 * 10^(-7))
  optimization <- optim(par = c(kinit, rinit), fn = function(x) - ll(x[1], x[2]))
  rhat <- optimization$par
  return(rhat)
}

## Mechanism to compute MLE given fs, Ys, epsilon
compute_rhat_iid <- function(frequencies, Ys, epsilon){
  ndata <- nrow(frequencies)
  distances <- rep(Inf, ndata)
  ll <- function(r) loglikelihood_cpp(1, r, Ys, frequencies, distances, epsilon, rho = 7.4 * 10^(-7))
  optimization <- optimize(f = function(x) - ll(x), interval = c(0, 1))
  rhat <- optimization$minimum
  return(rhat)
}

## Mechanism to generate Ys given fs, distances, k, r, epsilon
simulate_Ys_hmm <- function(frequencies, distances, k, r, epsilon){
  Ys <- simulate_data(frequencies, distances, k = k, r = r, epsilon, rho = 7.4 * 10^(-7))
  return(Ys)
}


# Fixed SNPs and positions 
num_SNPs = sapply(sites, function(x) nrow(hmmInput_freqs[[x]]))
mfxd = min(num_SNPs) # Lowest number of SNPs of sites 
site <- names(which.min(sapply(sites, function(x) nrow(hmmInput_freqs[[x]]))))
data_set <- hmmInput_freqs[[which.min(num_SNPs)]][,1:3]
subdata <- data_set[sample(mfxd, min(nsnps_nsamps[1,sites]), replace = F), ] # down sample full data set 
subdata <- subdata %>% arrange(chrom, pos)  # sort by chromosome and position  
subdata$dt <- c(diff(subdata$pos), Inf)
pos_change_chrom <- 1 + which(diff(subdata$chrom) != 0) # find places where chromosome changes
subdata$dt[pos_change_chrom-1] <- Inf
dfxd = subdata$dt

## Wrapper function to generate simulate data and calculate quantities of interest
wrapper_fun = function(site, subsample, model){
  
  # Extract data and subsample
  data_set <- hmmInput_freqs[[site]][,1:3] # extract frequencies and positions only 
  
  # Either subsample or not
  if(subsample){
    subdata <- data_set[sample(mfxd, min(nsnps_nsamps[1,sites]), replace = F), ] # down sample full data set 
    distances <- dfxd
  } else {
    subdata <- data_set 
    subdata <- subdata %>% arrange(chrom, pos)  # sort by chromosome and position  
    subdata$dt <- c(diff(subdata$pos), Inf)
    pos_change_chrom <- 1 + which(diff(subdata$chrom) != 0) # find places where chromosome changes
    subdata$dt[pos_change_chrom-1] <- Inf
    distances = subdata$dt
  }
  
  # Extract frequencies and distances 
  frequencies <- cbind(1-subdata$fs, subdata$fs)
  
  # Simulate the data 
  relatedness_measures = t(sapply(1:nsim, function(x){
    Ys = simulate_Ys_hmm(frequencies, distances, k = kfxd, r = rfxd, epsilon = efxd)
    if(model == 'hmm'){
      rhat = compute_rhat_hmm(frequencies, distances, Ys, efxd)[2]
    } else {
      rhat = compute_rhat_iid(frequencies, Ys, efxd)
    }
    ibs = mean(Ys[,1] == Ys[,2], na.rm = TRUE)
    results = c(ibs = ibs, rhat = rhat)
    return(results)
  }))
  
  # Calculate h constant
  h_constant = mean(frequencies^2 + (1-frequencies)^2) 
  
  # Return all
  return(list(relatedness_measures = relatedness_measures, h_constant = h_constant))
}


X1 = lapply(sites, wrapper_fun, subsample = T, model = 'hmm') 
X3 = lapply(sites, wrapper_fun, subsample = F, model = 'hmm')
names(X1) = names(X3) = sites
Xs = list(X1, X3) 
names(Xs) = c('HMM_redSNPs','HMM_allSNPs')
save(Xs, file = '../RData/IBS_sensitivity_sim.RData')

