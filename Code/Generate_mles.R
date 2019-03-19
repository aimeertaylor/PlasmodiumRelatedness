##########################################################
# Use k start value of 12 to avoid to match k = 0.5
##########################################################
rm(list = ls())
set.seed(1)
library(ggplot2)
library(dplyr)
library(Rcpp)
library(doParallel)
library(doRNG)
source("./simulate_data.R")
sourceCpp("./hmmloglikelihood.cpp")
registerDoParallel(cores = detectCores()-2)
epsilon <- 0.001 # Fix epsilon throughout
set.seed(1) # for reproducibility

## Mechanism to compute MLE given fs, Ys, epsilon
compute_rhat_iid <- function(frequencies, Ys, epsilon){
  ndata <- nrow(frequencies)
  distances <- rep(Inf, ndata)
  ll <- function(r) loglikelihood_cpp(1, r, Ys, frequencies, distances, epsilon, rho = 7.4 * 10^(-7))
  optimization <- optimize(f = function(x) - ll(x), interval = c(0, 1))
  rhat <- optimization$minimum
  return(rhat)
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
for (idataset in 1:length(dataset_names)){
  print(dataset_names[idataset])
  data_set <- hmmInput_freqs[[dataset_names[idataset]]]
  print(dim(data_set))
}

for (idataset in 1:length(dataset_names)){
  print(dataset_names[idataset])
  data_set <- hmmInput_freqs[[dataset_names[idataset]]]
  print(dim(data_set))
  if (dim(data_set)[2] < 100){
    next
  }
  individual_names <- names(data_set)[4:ncol(data_set)]
  nindividuals <- length(individual_names)
  name_combinations <- matrix(nrow = nindividuals*(nindividuals-1)/2, ncol = 2)
  count <- 0
  for (i in 1:(nindividuals-1)){
    for (j in (i+1):nindividuals){
      count <- count + 1
      name_combinations[count,1] <- individual_names[i]
      name_combinations[count,2] <- individual_names[j]
    }
  }
  
  mle_df <- foreach(icombination = 1:nrow(name_combinations), .combine = rbind) %dorng% {
    individual1 <- name_combinations[icombination,1]
    individual2 <- name_combinations[icombination,2]
    # let's focus on one pair of individuals
    i1 <- which(individual1 == names(data_set))
    i2 <- which(individual2 == names(data_set))
    subdata <- data_set[,c(1,2,3,i1,i2)]
    names(subdata) <- c("chrom", "pos", "fs", "Yi", "Yj")
    # sort by chromosome and position
    subdata <- subdata %>% arrange(chrom, pos) 
    frequencies <- cbind(1-subdata$fs, subdata$fs)
    subdata$dt <- c(diff(subdata$pos), Inf)
    pos_change_chrom <- 1 + which(diff(subdata$chrom) != 0) # find places where chromosome changes
    subdata$dt[pos_change_chrom-1] <- Inf
    krhat_hmm <- compute_rhat_hmm(frequencies, subdata$dt, cbind(subdata$Yi, subdata$Yj), epsilon)
    rhat_iid <- compute_rhat_iid(frequencies, cbind(subdata$Yi, subdata$Yj), epsilon)
    ll_iid <- function(r) loglikelihood_cpp(1, r, cbind(subdata$Yi, subdata$Yj), frequencies, rep(Inf, dim(frequencies)[1]), epsilon, rho = 7.4 * 10^(-7))
    ll_hmm <- function(k, r) loglikelihood_cpp(k, r, cbind(subdata$Yi, subdata$Yj), frequencies, subdata$dt, epsilon, rho = 7.4 * 10^(-7))
    # compute second order derivative of -log-likelihood at MLE, divided by large sample size
    FIM_iid <- numDeriv::hessian(function(x) -ll_iid(x), rhat_iid)/(dim(frequencies)[1])
    FIM_hmm <- numDeriv::hessian(function(x) -ll_hmm(x[1], x[2]), krhat_hmm)/(dim(frequencies)[1])
    var_hmm <- try(solve(FIM_hmm))
    if (inherits(var_hmm, "try-error")){
      var_hmm <- matrix(NA, 2, 2)
    }
    var_iid <- 1/(FIM_iid[1,1])
    var_k_hmm <- var_hmm[1,1]
    var_r_hmm <- var_hmm[2,2]
    data.frame(individual1 = individual1, individual2 = individual2, 
               rhat_iid = rhat_iid, khat_hmm = krhat_hmm[1], rhat_hmm = krhat_hmm[2],
               var_iid = var_iid, var_k_hmm = var_k_hmm, var_r_hmm = var_r_hmm)
  }
  save(mle_df, file = paste0("./RData/", dataset_names[idataset], "_mles.RData"))
}

