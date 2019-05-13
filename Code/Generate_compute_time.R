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
nrepeats <- 1000 # number of repeats 
kfixed <- 12 
rfixed <- 0.5 
load('../RData/tables_many_repeats_m.RData')
dataset_names = names(tables_many_repeats_m_dataset)

## Mechanism to generate Ys given fs, distances, k, r, epsilon
simulate_Ys_hmm <- function(frequencies, distances, k, r, epsilon){
  Ys <- simulate_data(frequencies, distances, k = k, r = r, epsilon, rho = 7.4 * 10^(-7))
  return(Ys)
}

## Mechanism to compute MLE given fs, Ys, epsilon
compute_rhat_iid <- function(frequencies, Ys, epsilon){
  ndata <- nrow(frequencies)
  distances <- rep(Inf, ndata)
  ll <- function(r) loglikelihood_cpp(1, r, Ys, frequencies, distances, epsilon, rho = 7.4 * 10^(-7))
  ptm <- proc.time()
  optimization <- optimize(f = function(x) - ll(x), interval = c(0, 1))
  time = proc.time() - ptm
  names(time) = paste0('iid_',names(time))
  return(time = time[1:3]) # Return user, self, elapsed
}

## Mechanism to compute MLE given fs, distances, Ys, epsilon
compute_rhat_hmm <- function(frequencies, distances, Ys, epsilon){
  ndata <- nrow(frequencies)
  ll <- function(k, r) loglikelihood_cpp(k, r, Ys, frequencies, distances, epsilon, rho = 7.4 * 10^(-7))
  ptm <- proc.time()
  optimization <- optim(par = c(kfixed, 0.5), fn = function(x) - ll(x[1], x[2]))
  rhat <- optimization$par
  time = proc.time() - ptm
  names(time) = paste0('hmm_',names(time))
  names(optimization$counts) = paste0('hmm_',names(optimization$counts))
  return(c(time[1:3], # Return user, self, elapsed
           optimization$counts[1])) # Return calls to fn (calls to gr are all NA)
}


#===================================================
# Generate average times and save
#===================================================
Mean_time = t(sapply(dataset_names, function(iname){
  dataset_ <- tables_many_repeats_m_dataset[[iname]]
  frequencies <- cbind(1-dataset_$fs, dataset_$fs)
  Times <- foreach (irepeat = 1:nrepeats, .combine = rbind) %dorng% {
    Ys <- simulate_Ys_hmm(frequencies, dataset_$dt, kfixed, rfixed, epsilon)
    c(compute_rhat_iid(frequencies, Ys, epsilon),compute_rhat_hmm(frequencies, dataset_$dt, Ys, epsilon))
  }
  colMeans(Times)}))
save(Mean_time,file = '../RData/Mean_time.RData' )



#===================================================
# Post-process
#===================================================
f_strategy = "Proportional to MAF" # Choose a strategy
inds_i = grepl(f_strategy, rownames(Mean_time))

# Separate models
inds_iid =  grepl('iid', colnames(Mean_time))
inds_hmm =  grepl('hmm', colnames(Mean_time))

# Reformat rownames 
X_iid = Mean_time[inds_i,inds_iid]
rownames(X_iid) = gsub('m','',do.call(rbind, strsplit(rownames(X_iid), split = '_'))[,1])

# Reformat rownames 
X_hmm = Mean_time[inds_i,inds_hmm]
rownames(X_hmm) = gsub('m','',do.call(rbind, strsplit(rownames(X_hmm), split = '_'))[,1])

# Print latex code
kable(format(X_iid,digits = 4,drop0trailing = F), format = 'latex')
kable(format(X_hmm,digits = 4,drop0trailing = F), format = 'latex')



