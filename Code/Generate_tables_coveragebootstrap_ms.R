### This script creates tables using both the independence model and the HMM 
### (including genotyping errors and allowing multiallelic loci)
### for various sample size, various rs, various fs, various Ks


## Set up
rm(list = ls())
set.seed(1) # for reproducibility
library(dplyr)
library(Rcpp)
library(doParallel)
library(doRNG)
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
samplesizes <- c(96*c(0.25, 1:5)) 
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


#============================================================
# Different sample sizes fixed K = 2 (Thai), specifically: 
# 1) Different sample sizes and rs given fixed K = 2 and k = 12
# 2) Different sample sizes and ks given fixed K = 2 and r = 0.5
#============================================================
# Mechanism to generate f and K: Option 1 
load("../RData/hmmInput_freqs.RData")
data_ <- hmmInput_freqs$TM_WGS[,c(1,2,3)] # Using Thai frequencies and positions as a template
head(data_)
nSNPs <- nrow(data_)
all_frequencies <- cbind(1-data_$fs, data_$fs)
data_$dt <- c(diff(data_$pos), Inf)
pos_change_chrom <- 1 + which(diff(data_$chrom) != 0) # find places where chromosome changes
data_$dt[pos_change_chrom-1] <- Inf
maf <- pmin(data_$fs, 1 - data_$fs)

# Sample allele frequencies with probability "probs" and positions uniformally at random
simulate_ <- function(samplesize, probs, distances = NULL){
  # 1) Sample data set
  indices <- sample(x = 1:nrow(data_), size = samplesize, replace = FALSE)
  dataset_ <- data_[sort(indices),]
  # 2) Format distances
  if(!is.null(distances)){  
    if(length(distances) != samplesize){stop('Distances have the wrong length')}
    dataset_$dt <- distances
  } else { # If no distances then return distances based on data set
    dataset_$dt <- c(diff(dataset_$pos), Inf)
    pos_change_chrom <- 1 + which(diff(dataset_$chrom) != 0) # find places where chromosome changes
    dataset_$dt[pos_change_chrom-1] <- Inf
  }
  # 3) Sample frequencies: no need to sort 
  dataset_$fs <- data_$fs[sample(x = 1:nrow(data_), size = samplesize, prob = probs, replace = FALSE)]
  return(dataset_)
}

# Generate distances s.t. they are the same over different fs_strategy
Distances <- lapply(samplesizes, function(x){simulate_(x, probs = rep(1,nSNPs), distances = NULL)$dt})
names(Distances) <- samplesizes

fs_strategy <- list("Proportional to MAF" = maf)

#----------------------------------------------------------------
# 1) Different sample sizes and rs given fixed K = 2 and k = 12
#----------------------------------------------------------------
tables_many_repeats_mr_coverage_iid <- array(dim = c(length(samplesizes), length(rvaried), length(fs_strategy)),dimnames = list(samplesizes, rvaried, names(fs_strategy)))
tables_many_repeats_mr_coverage_hmm <- tables_many_repeats_mr_coverage_iid
tables_many_repeats_mr_length_iid <- array(dim = c(length(samplesizes), length(rvaried), length(fs_strategy)),dimnames = list(samplesizes, rvaried, names(fs_strategy)))
tables_many_repeats_mr_length_hmm <- tables_many_repeats_mr_length_iid
tables_many_repeats_mr_dataset <- list()


system.time(
  for(ilayer in 1:length(fs_strategy)){
    for (irow in 1:length(samplesizes)){
      # draw frequencies / positions 'once'
      dataset_ <- simulate_(samplesizes[irow],
                            probs = fs_strategy[[ilayer]],
                            distances = Distances[[irow]])
      frequencies <- cbind(1-dataset_$fs, dataset_$fs)
      dataset_name <- sprintf('m%s_fs:%s', samplesizes[irow], names(fs_strategy)[[ilayer]])
      tables_many_repeats_mr_dataset[[dataset_name]] <- dataset_
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
        tables_many_repeats_mr_length_iid[irow, icolumn, ilayer] <- mean(rhats_$len_iid, na.rm = TRUE)
        tables_many_repeats_mr_length_hmm[irow, icolumn, ilayer] <- mean(rhats_$len_hmm, na.rm = TRUE)
        tables_many_repeats_mr_coverage_iid[irow, icolumn, ilayer] <- mean(rhats_$cover_iid, na.rm = TRUE)
        tables_many_repeats_mr_coverage_hmm[irow, icolumn, ilayer] <- mean(rhats_$cover_hmm, na.rm = TRUE)
      }
    }
  }
)

save(tables_many_repeats_mr_length_iid,
     tables_many_repeats_mr_length_hmm,
     tables_many_repeats_mr_coverage_iid,
     tables_many_repeats_mr_coverage_hmm,
     tables_many_repeats_mr_dataset,
     file = '../RData/tables_many_repeats_mr_lengthcoverage_bootstrap.RData' )

tables_many_repeats_mr_length_iid
tables_many_repeats_mr_length_hmm
tables_many_repeats_mr_coverage_iid
tables_many_repeats_mr_coverage_hmm



#----------------------------------------------------------------
# 2) Different sample sizes and ks given fixed K = 2 and r = 0.5
#----------------------------------------------------------------
tables_many_repeats_mk_coverage_iid <- array(dim = c(length(samplesizes), length(kvaried), length(fs_strategy)),dimnames = list(samplesizes, kvaried, names(fs_strategy)))
tables_many_repeats_mk_coverage_hmm <- tables_many_repeats_mk_coverage_iid
tables_many_repeats_mk_length_iid <- array(dim = c(length(samplesizes), length(kvaried), length(fs_strategy)),dimnames = list(samplesizes, kvaried, names(fs_strategy)))
tables_many_repeats_mk_length_hmm <- tables_many_repeats_mk_length_iid
tables_many_repeats_mk_dataset <- list()


system.time(
  for(ilayer in 1:length(fs_strategy)){
    for (irow in 1:length(samplesizes)){
      # draw frequencies / positions 'once'
      dataset_ <- simulate_(samplesizes[irow],
                            probs = fs_strategy[[ilayer]],
                            distances = Distances[[irow]])
      frequencies <- cbind(1-dataset_$fs, dataset_$fs)
      dataset_name <- sprintf('m%s_fs:%s', samplesizes[irow], names(fs_strategy)[[ilayer]])
      tables_many_repeats_mk_dataset[[dataset_name]] <- dataset_
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
        tables_many_repeats_mk_length_iid[irow, icolumn, ilayer] <- mean(rhats_$len_iid, na.rm = TRUE)
        tables_many_repeats_mk_length_hmm[irow, icolumn, ilayer] <- mean(rhats_$len_hmm, na.rm = TRUE)
        tables_many_repeats_mk_coverage_iid[irow, icolumn, ilayer] <- mean(rhats_$cover_iid, na.rm = TRUE)
        tables_many_repeats_mk_coverage_hmm[irow, icolumn, ilayer] <- mean(rhats_$cover_hmm, na.rm = TRUE)
      }
    }
  }
)

save(tables_many_repeats_mk_length_iid,
     tables_many_repeats_mk_length_hmm,
     tables_many_repeats_mk_coverage_iid,
     tables_many_repeats_mk_coverage_hmm,
     tables_many_repeats_mk_dataset,
     file = '../RData/tables_many_repeats_mk_lengthcoverage_bootstrap.RData')

tables_many_repeats_mk_length_iid
tables_many_repeats_mk_length_hmm
tables_many_repeats_mk_coverage_iid
tables_many_repeats_mk_coverage_hmm


