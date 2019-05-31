### This script times inference for the iid model and hmm
### using data sets stored in tables_many_repeats_m.RData
### Data are simulated under the hmm using fixed r and k 

## Set up
rm(list = ls())
set.seed(1)
library(Rcpp)
library(doParallel)
library(doRNG)
require(kableExtra) # to print as table
source("./simulate_data.R")
sourceCpp("./hmmloglikelihood.cpp")
registerDoParallel(cores = detectCores()-2)
epsilon <- 0.001 # Fix epsilon throughout
set.seed(1) # for reproducibility
nrepeats <- 500 # number of repeats 
kfixed <- 8
rfixed <- 0.5 
rinit = rfixed
kinit = kfixed
load('../RData/tables_many_repeats_m.RData')
dataset_names = names(tables_many_repeats_m_dataset)


## Mechanism to generate Ys given fs, distances, k, r, epsilon
simulate_Ys_hmm <- function(frequencies, distances, k, r, epsilon){
  Ys <- simulate_data(frequencies, distances, k = k, r = r, epsilon, rho = 7.4 * 10^(-7))
  return(Ys)
}

## Mechanism to compute MLE given fs, Ys, epsilon
compute_rhat_iid_optimize <- function(frequencies, Ys, epsilon){
  ndata <- nrow(frequencies)
  distances <- rep(Inf, ndata)
  ll <- function(r) loglikelihood_cpp(1, r, Ys, frequencies, distances, epsilon, rho = 7.4 * 10^(-7))
  ptm <- proc.time()
  optimization <- optimize(f = function(x) - ll(x), interval = c(0, 1))
  time = proc.time() - ptm
  names(time) = paste0('iid_optimize_',names(time))
  return(time = time[1:3]) # Return user, self, elapsed
}

## Mechanism to compute MLE given fs, Ys, epsilon
compute_rhat_iid_optim <- function(frequencies, Ys, epsilon){
  ndata <- nrow(frequencies)
  distances <- rep(Inf, ndata)
  ll <- function(r) loglikelihood_cpp(1, r, Ys, frequencies, distances, epsilon, rho = 7.4 * 10^(-7))
  ptm <- proc.time()
  optimization <- optim(par = rinit, fn = function(x) - ll(x))
  rhat <- optimization$par
  time = proc.time() - ptm
  names(time) = paste0('iid_optim_',names(time))
  names(optimization$counts) = paste0('iid_optim_',names(optimization$counts))
  return(c(time[1:3], # Return user, self, elapsed
           optimization$counts[1], # Return calls to fn (calls to gr are all NA)
           'convergence_code' = optimization$convergence)) # Return convergence code
}

## Mechanism to compute MLE given fs, distances, Ys, epsilon
compute_rhat_hmm <- function(frequencies, distances, Ys, epsilon){
  ndata <- nrow(frequencies)
  ll <- function(k, r) loglikelihood_cpp(k, r, Ys, frequencies, distances, epsilon, rho = 7.4 * 10^(-7))
  ptm <- proc.time()
  optimization <- optim(par = c(kinit, rinit), fn = function(x) - ll(x[1], x[2]))
  rhat <- optimization$par
  time = proc.time() - ptm
  names(time) = paste0('hmm_',names(time))
  names(optimization$counts) = paste0('hmm_',names(optimization$counts))
  return(c(time[1:3], # Return user, self, elapsed
           optimization$counts[1], # Return calls to fn (calls to gr are all NA)
           'convergence_code' = optimization$convergence)) # Return convergence code
}


#===================================================
# Generate average times and save
#===================================================
Comp_results = lapply(dataset_names, function(iname){
  dataset_ <- tables_many_repeats_m_dataset[[iname]]
  frequencies <- cbind(1-dataset_$fs, dataset_$fs)
  Times <- foreach (irepeat = 1:nrepeats, .combine = rbind) %dorng% {
    Ys <- simulate_Ys_hmm(frequencies, dataset_$dt, kfixed, rfixed, epsilon)
    c(compute_rhat_iid_optimize(frequencies, Ys, epsilon), 
      compute_rhat_iid_optim(frequencies, Ys, epsilon),
      compute_rhat_hmm(frequencies, dataset_$dt, Ys, epsilon))}
  
  inds = grepl('convergence_code', colnames(Times))
  con_code_iid = table(Times[,which(inds)[1]]) 
  con_code_hmm = table(Times[,which(inds)[2]])
  list(Times = colMeans(Times[,!inds]), con_code_iid = con_code_iid, con_code_hmm = con_code_hmm)
})
names(Comp_results) = dataset_names
save(Comp_results,file = '../RData/Comp_results.RData')

