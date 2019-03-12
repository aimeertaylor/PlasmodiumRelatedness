### This script creates tables using both the independence model and the HMM 
### (including genotyping errors and allowing multiallelic loci)
### for various sample size, various rs, various fs, various Ks

## Set up
rm(list = ls()) # for reproducibility
set.seed(1)
library(dplyr)
library(Rcpp)
library(doParallel)
library(doRNG)
library(MCMCpack) # For dirichlet
source("./simulate_data.R")
sourceCpp("./hmmloglikelihood.cpp")
registerDoParallel(cores = detectCores()-2)
epsilon <- 0.001 # Fix epsilon throughout
kfixed <- 12
rfixed <- 0.5
mfixed <- 24
rvaried <- seq(from = 0.01, to = 0.99, length.out = 5)
kvaried <- c(1,10,100)
Kvaried <- seq(2,10,2)
alphas = c(100, 1) # High versus low entropy
nrepeats <- 500 # number of times data are simulated for a given set of parameters 
nboot <- 500 # number of bootstrap repeats

## Mechanism to generate Ys given fs, r, epsilon
simulate_Ys_iid <- function(frequencies, r, epsilon){
  Ys <- simulate_data(frequencies, rep(Inf, nrow(frequencies)), 1, r, epsilon, rho = 7.4 * 10^(-7))
  return(Ys)
}

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
  optimization <- optimize(f = function(x) - ll(x), interval = c(0, 1))
  rhat <- optimization$minimum
  return(rhat)
}

## Mechanism to compute MLE given fs, distances, Ys, epsilon
compute_rhat_hmm <- function(frequencies, distances, Ys, epsilon){
  ndata <- nrow(frequencies)
  ll <- function(k, r) loglikelihood_cpp(k, r, Ys, frequencies, distances, epsilon, rho = 7.4 * 10^(-7))
  optimization <- optim(par = c(kfixed, 0.5), fn = function(x) - ll(x[1], x[2]))
  rhat <- optimization$par
  return(rhat)
}

# Mechanism to generate fs given K_t and alpha
simulate_fs = function(K, alpha, samplesize){
  fs = rdirichlet(n = samplesize, rep(alpha, K))
  fs <- t(apply(fs, 1, function(v) (v + 1e-8)/sum(v + 1e-8)))
  return(fs)
}

# Generate distances on time only since m is fixed
load("../RData/hmmInput_freqs.RData")
data_ <- hmmInput_freqs$TM_WGS[,c(1,2,3)] # Using Thai frequencies and positions as a template
indices <- sample(x = 1:nrow(data_), size = mfixed, replace = FALSE)
dataset_ <- data_[sort(indices),]
dataset_$dt <- c(diff(dataset_$pos), Inf)
pos_change_chrom <- 1 + which(diff(dataset_$chrom) != 0) # find places where chromosome changes
dataset_$dt[pos_change_chrom-1] <- Inf
Distances <- dataset_$dt
names(Distances) <- mfixed

#----------------------------------------------------------------
# 1) Different Kt and rs given fixed m = 24 and k = 12
#----------------------------------------------------------------
# Save stores 
tables_many_repeats_Kr_coverage_iid <- array(dim = c(length(Kvaried), length(rvaried), length(alphas)),
                                             dimnames = list(Kvaried, rvaried, alphas))
tables_many_repeats_Kr_coverage_hmm <- tables_many_repeats_Kr_coverage_iid
tables_many_repeats_Kr_length_iid <- array(dim = c(length(Kvaried), length(rvaried), length(alphas)),
                                           dimnames = list(Kvaried, rvaried, alphas))
tables_many_repeats_Kr_length_hmm <- tables_many_repeats_Kr_length_iid
tables_many_repeats_Kr_dataset <- list('dts' = Distances)

system.time(
  for(ilayer in 1:length(alphas)){
    for (irow in 1:length(Kvaried)){
      frequencies <- simulate_fs(Kvaried[irow], alphas[ilayer], mfixed)
      freq_name <- sprintf('K%s_alpha:%s', Kvaried[irow], alphas[[ilayer]])
      tables_many_repeats_Kr_dataset[[freq_name]] <- frequencies
      for (icolumn in 1:length(rvaried)){
        rhats_ <- foreach (irepeat = 1:nrepeats, .combine = rbind) %dorng% {
          Ys <- simulate_Ys_hmm(frequencies, dataset_$dt, kfixed, rvaried[icolumn], epsilon)
          ndata <- nrow(frequencies)
          distances <- rep(Inf, ndata)
          ll_iid <- function(r) loglikelihood_cpp(1, r, Ys, frequencies, distances, epsilon, rho = 7.4 * 10^(-7))
          optimization <- optimize(f = function(x) - ll_iid(x), interval = c(0, 1))
          rhat_iid <- optimization$minimum
          #
          rhats_iid_boot <- rep(0, nboot)
          for (iboot in 1:nboot){
            Ys_boot <- simulate_Ys_iid(frequencies, rhat_iid, epsilon)
            ll_iid <- function(r) loglikelihood_cpp(1, r, Ys_boot, frequencies, distances, epsilon, rho = 7.4 * 10^(-7))
            optimization <- optimize(f = function(x) - ll_iid(x), interval = c(0, 1))
            rhats_iid_boot[iboot] <- optimization$minimum
          }
          ci_iid <- as.numeric(quantile(rhats_iid_boot, probs = c(0.025, 0.975)))
          
          ll_hmm <- function(k, r) loglikelihood_cpp(k, r, Ys, frequencies, dataset_$dt, epsilon, rho = 7.4 * 10^(-7))
          optimization <- optim(par = c(kfixed, 0.5), fn = function(x) - ll_hmm(x[1], x[2]))
          rhat_hmm <- optimization$par
          rhats_hmm_boot <- rep(0, nboot)
          for (iboot in 1:nboot){
            Ys_boot <- simulate_Ys_hmm(frequencies, dataset_$dt, rhat_hmm[1], rhat_hmm[2], epsilon)
            ll_hmm <- function(k, r) loglikelihood_cpp(k, r, Ys_boot, frequencies, dataset_$dt, epsilon, rho = 7.4 * 10^(-7))
            optimization <- optim(par = c(kfixed, 0.5), fn = function(x) - ll_hmm(x[1], x[2]))
            rhats_hmm_boot[iboot] <- optimization$par[2]
          }
          ci_hmm <- as.numeric(quantile(rhats_hmm_boot, probs = c(0.025, 0.975)))
          data.frame(len_iid = ci_iid[2] - ci_iid[1], cover_iid = (ci_iid[1] < rvaried[icolumn]) & (ci_iid[2] > rvaried[icolumn]),
                     len_hmm = ci_hmm[2] - ci_hmm[1], cover_hmm = (ci_hmm[1] < rvaried[icolumn]) & (ci_hmm[2] > rvaried[icolumn]))
        }
        tables_many_repeats_Kr_length_iid[irow, icolumn, ilayer] <- mean(rhats_$len_iid, na.rm = TRUE)
        tables_many_repeats_Kr_length_hmm[irow, icolumn, ilayer] <- mean(rhats_$len_hmm, na.rm = TRUE)
        tables_many_repeats_Kr_coverage_iid[irow, icolumn, ilayer] <- mean(rhats_$cover_iid, na.rm = TRUE)
        tables_many_repeats_Kr_coverage_hmm[irow, icolumn, ilayer] <- mean(rhats_$cover_hmm, na.rm = TRUE)
      }
    }
  }
)
save(tables_many_repeats_Kr_length_iid,
     tables_many_repeats_Kr_length_hmm,
     tables_many_repeats_Kr_coverage_iid,
     tables_many_repeats_Kr_coverage_hmm,
     tables_many_repeats_Kr_dataset,
     file = '../RData/tables_many_repeats_Kr_lengthcoverage_bootstrap.RData' )

tables_many_repeats_Kr_length_iid
tables_many_repeats_Kr_length_hmm
tables_many_repeats_Kr_coverage_iid
tables_many_repeats_Kr_coverage_hmm

#----------------------------------------------------------------
# 2) Different sample sizes and ks given fixed K = 2 and r = 0.5
#----------------------------------------------------------------
# Save stores 
tables_many_repeats_Kk_coverage_iid <- array(dim = c(length(Kvaried), length(kvaried), length(alphas)),
                                              dimnames = list(Kvaried, kvaried, names(alphas)))
tables_many_repeats_Kk_coverage_hmm <- tables_many_repeats_Kk_coverage_iid
tables_many_repeats_Kk_length_iid <- array(dim = c(length(Kvaried), length(kvaried), length(alphas)),
                                            dimnames = list(Kvaried, kvaried, names(alphas)))
tables_many_repeats_Kk_length_hmm <- tables_many_repeats_Kk_length_iid
tables_many_repeats_Kk_dataset <- list('dts' = Distances)

system.time(
  for(ilayer in 1:length(alphas)){
    for (irow in 1:length(Kvaried)){
      # draw frequencies 'once'
      # draw frequencies 'once'
      frequencies <- simulate_fs(Kvaried[irow], alphas[ilayer], mfixed)
      freq_name <- sprintf('K%s_alpha:%s', Kvaried[irow], alphas[[ilayer]])
      tables_many_repeats_Kk_dataset[[freq_name]] <- frequencies
      for (icolumn in 1:length(kvaried)){
        rhats_ <- foreach (irepeat = 1:nrepeats, .combine = rbind) %dorng% {
          Ys <- simulate_Ys_hmm(frequencies, dataset_$dt, kvaried[icolumn], rfixed, epsilon)
          ndata <- nrow(frequencies)
          distances <- rep(Inf, ndata)
          ll_iid <- function(r) loglikelihood_cpp(1, r, Ys, frequencies, distances, epsilon, rho = 7.4 * 10^(-7))
          optimization <- optimize(f = function(x) - ll_iid(x), interval = c(0, 1))
          rhat_iid <- optimization$minimum
          #
          rhats_iid_boot <- rep(0, nboot)
          for (iboot in 1:nboot){
            Ys_boot <- simulate_Ys_iid(frequencies, rhat_iid, epsilon)
            ll_iid <- function(r) loglikelihood_cpp(1, r, Ys_boot, frequencies, distances, epsilon, rho = 7.4 * 10^(-7))
            optimization <- optimize(f = function(x) - ll_iid(x), interval = c(0, 1))
            rhats_iid_boot[iboot] <- optimization$minimum
          }
          ci_iid <- as.numeric(quantile(rhats_iid_boot, probs = c(0.025, 0.975)))
          
          ll_hmm <- function(k, r) loglikelihood_cpp(k, r, Ys, frequencies, dataset_$dt, epsilon, rho = 7.4 * 10^(-7))
          optimization <- optim(par = c(kfixed, 0.5), fn = function(x) - ll_hmm(x[1], x[2]))
          rhat_hmm <- optimization$par
          rhats_hmm_boot <- rep(0, nboot)
          for (iboot in 1:nboot){
            Ys_boot <- simulate_Ys_hmm(frequencies, dataset_$dt, rhat_hmm[1], rhat_hmm[2], epsilon)
            ll_hmm <- function(k, r) loglikelihood_cpp(k, r, Ys_boot, frequencies, dataset_$dt, epsilon, rho = 7.4 * 10^(-7))
            optimization <- optim(par = c(kfixed, 0.5), fn = function(x) - ll_hmm(x[1], x[2]))
            rhats_hmm_boot[iboot] <- optimization$par[2]
          }
          ci_hmm <- as.numeric(quantile(rhats_hmm_boot, probs = c(0.025, 0.975)))
          data.frame(len_iid = ci_iid[2] - ci_iid[1], cover_iid = (ci_iid[1] < rfixed) & (ci_iid[2] > rfixed),
                     len_hmm = ci_hmm[2] - ci_hmm[1], cover_hmm = (ci_hmm[1] < rfixed) & (ci_hmm[2] > rfixed))
        }
        tables_many_repeats_Kk_length_iid[irow, icolumn, ilayer] <- mean(rhats_$len_iid, na.rm = TRUE)
        tables_many_repeats_Kk_length_hmm[irow, icolumn, ilayer] <- mean(rhats_$len_hmm, na.rm = TRUE)
        tables_many_repeats_Kk_coverage_iid[irow, icolumn, ilayer] <- mean(rhats_$cover_iid, na.rm = TRUE)
        tables_many_repeats_Kk_coverage_hmm[irow, icolumn, ilayer] <- mean(rhats_$cover_hmm, na.rm = TRUE)
      }
    }
  }
)
save(tables_many_repeats_Kk_length_iid,
     tables_many_repeats_Kk_length_hmm,
     tables_many_repeats_Kk_coverage_iid,
     tables_many_repeats_Kk_coverage_hmm,
     tables_many_repeats_Kk_dataset,
     file = '../RData/tables_many_repeats_Kk_lengthcoverage_bootstrap.RData' )

tables_many_repeats_Kk_length_iid
tables_many_repeats_Kk_length_hmm
tables_many_repeats_Kk_coverage_iid
tables_many_repeats_Kk_coverage_hmm