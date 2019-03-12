### This script creates tables using both the i.i.d model and the HMM 
### (including genotyping errors and allowing multiallelic loci)
### 
### The tables show the average mean square errors (RMSEs)
### for various sample size, various rs, various fs, various Ks
###
### Terminology: f refers to a frequency of an allele at each position in the genome, 
### where there are K possible alleles at a given position 
### r refers to relatedness


## Set up
rm(list = ls())
set.seed(1)
library(ggplot2)
library(dplyr)
library(Rcpp)
library(doParallel)
library(doRNG)
library(MCMCpack) # For dirichlet
source("./simulate_data.R")
sourceCpp("./hmmloglikelihood.cpp")
registerDoParallel(cores = detectCores()-2)
epsilon <- 0.001 # Fix epsilon throughout
rgrid <- seq(from = 0.01, to = 0.99, length.out = 5)
samplesizes = c(96*c(0.25, 1:5)) # add maximum sample size when working with real data
set.seed(1) # for reproducibility
nrepeats <- 1000 # number of repeats 
kfixed <- 12 
rfixed <- 0.5 
littleks <- c(1,5,10,50,1e3)

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


## Mechanism to generate f and K
## Option 1: get f w/wo positions for K = 2 by sampling from a real data set with probability proportional to probs (MAF vs uniformaly at random)
## Option 2: simulate by fixing K (alt: draw from distribution on integers), then drawing fs from a Dirichlet distribution (high vs and low entropy)
##           get positions uniformly at random from real data (alt: draw from positions of core reference genome)
## We use Option 1 for exploring CIs for different sample sizes and rs given fixed K = 2
## We use Option 2 for exploring CIs for different Ks and rs given fixed m = 96 

#====================================================
# Different sample sizes and rs given fixed K = 2
#====================================================
# Mechanism to generate f and K: Option 1 
load("../RData/hmmInput_freqs.RData")
data_ <- hmmInput_freqs$TM_WGS[,c(1,2,3)] # Using Thai frequencies and positions as a template
# data_ <- hmmInput_freqs$Colombia[,c(1,2,3)]
head(data_)
K <- 2
nSNPs <- nrow(data_)
all_frequencies <- cbind(1-data_$fs, data_$fs)
any(all_frequencies == 0 | all_frequencies == 1) # Check all polymorphic
data_$dt <- c(diff(data_$pos), Inf)
pos_change_chrom <- 1 + which(diff(data_$chrom) != 0) # find places where chromosome changes
data_$dt[pos_change_chrom-1] <- Inf
maf <- pmin(data_$fs, 1 - data_$fs)
hist(maf, col = 'gray', freq = F)
curve(dbeta(x, 1.2, 20), from = 0, to = 0.5, add = TRUE, col = 'blue') 
legend('top', fill = c('gray', 'blue'), bty = 'n',legend = c('Thai MAF','Beta(1.2, 20)'))


# Sample allele frequencies with probability "probs" and positions uniformally at random
simulate_ <- function(samplesize, probs, distances = NULL){
  # 1) Sample data set
  indices <- sample(x = 1:nrow(data_), size = samplesize, replace = FALSE)
  dataset_ <- data_[sort(indices),]
  # 2) Format distances
  if(!is.null(distances)){ 
    if(length(distances) != samplesize){stop('Distances have the wrong length')}
    dataset_$dt <- distances
  } else {
    dataset_$dt <- c(diff(dataset_$pos), Inf)
    pos_change_chrom <- 1 + which(diff(dataset_$chrom) != 0) # find places where chromosome changes
    dataset_$dt[pos_change_chrom-1] <- Inf
  }
  # 3) Sample frequencies: no need to sort 
  dataset_$fs <- data_$fs[sample(x = 1:nrow(data_), size = samplesize, prob = probs, replace = FALSE)]
  return(dataset_)
}


## First approach: repeat many MLE calculations and use empirical variance of obtained MLEs
## Specify array with row per samplesize, column per r, layer per fs_strategy 
## Result should be indifferent to fs_strategy at sample size that uses all available data
## 
fs_strategy <- list("Proportional to MAF" = maf, 
                    "Uniformly at random" = rep(1,nSNPs))
tables_many_repeats_m_iid <- array(dim = c(length(samplesizes), length(rgrid), length(fs_strategy)),
                                   dimnames = list(samplesizes, rgrid, names(fs_strategy)))
tables_many_repeats_m_hmm <- tables_many_repeats_m_iid
tables_many_repeats_m_dataset <- list()
tables_many_repeats_m_thetas <- list()

# Generate distances s.t. they are the same over different fs_strategy
Distances <- lapply(samplesizes, function(x){simulate_(x, probs = rep(1,nSNPs), distances = NULL)$dt})
names(Distances) <- samplesizes

system.time(
  for(ilayer in 1:length(fs_strategy)){
    for (irow in 1:length(samplesizes)){
      # draw frequencies / positions 'once'
      dataset_ <- simulate_(samplesizes[irow],
                            probs = fs_strategy[[ilayer]],
                            distances = Distances[[irow]])
      frequencies <- cbind(1-dataset_$fs, dataset_$fs)
      dataset_name <- sprintf('m%s_fs:%s', samplesizes[irow], names(fs_strategy)[[ilayer]])
      tables_many_repeats_m_dataset[[dataset_name]] <- dataset_
      thetastore <- list()
      for (icolumn in 1:length(rgrid)){
        rhats_ <- foreach (irepeat = 1:nrepeats, .combine = rbind) %dorng% {
          Ys <- simulate_Ys_hmm(frequencies, dataset_$dt, kfixed, rgrid[icolumn], epsilon)
          c(compute_rhat_iid(frequencies, Ys, epsilon), compute_rhat_hmm(frequencies, dataset_$dt, Ys, epsilon))
        }
        tables_many_repeats_m_iid[irow, icolumn, ilayer] <- sqrt(mean((rhats_[,1] - rgrid[icolumn])^2))
        tables_many_repeats_m_hmm[irow, icolumn, ilayer] <- sqrt(mean((rhats_[,3] - rgrid[icolumn])^2))
        thetastore[[as.character(rgrid[icolumn])]] <- matrix(rhats_, ncol = 3)
      }
      tables_many_repeats_m_thetas[[dataset_name]] <- thetastore
    }
  }
)
save(tables_many_repeats_m_iid,
     tables_many_repeats_m_hmm,
     tables_many_repeats_m_dataset,
     tables_many_repeats_m_thetas,
     file = '../RData/tables_many_repeats_m.RData' )
load('../RData/tables_many_repeats_m.RData')
# show some of the results
tables_many_repeats_m_iid[,,1]
tables_many_repeats_m_hmm[,,1]



# ## Second approach: compute asymptotic variance of MLE using lots of data 
# ## Not JUSTIFIED
# ## and use asymptotic approximation to get variance for different finite sample sizes
# ## not justified when r is at the boundary (at 0 or 1)
# 
# large_sample_size <- 1e3
# table_asymptotics_m_iid <- array(dim = c(length(samplesizes), length(rgrid)), dimnames = list(samplesizes, rgrid))
# table_asymptotics_m_hmm <- table_asymptotics_m_iid
# for (icolumn in 1:length(rgrid)){
#   asymptvariances <- foreach(irep = 1:nrepeats, .combine = rbind) %dorng% {
#     dataset_ <- simulate_(large_sample_size, probs = fs_strategy[[1]]) # Resample distances each time
#     frequencies <- cbind(1-dataset_$fs, dataset_$fs)
#     Ys_iid <- simulate_Ys_iid(frequencies, rgrid[icolumn], epsilon)
#     Ys_hmm <- simulate_Ys_hmm(frequencies, dataset_$dt, kfixed, rgrid[icolumn], epsilon)
#     # define log likelihood
#     ll_iid <- function(r) loglikelihood_cpp(1, r, Ys_iid, frequencies, rep(Inf, large_sample_size), epsilon, rho = 7.4 * 10^(-7))
#     mle_iid <- compute_rhat_iid(frequencies, Ys_iid, epsilon)
#     # # Check cannot do because of kfixed
#     ll_hmm <- function(k, r) loglikelihood_cpp(k, r, Ys_hmm, frequencies, dataset_$dt, epsilon, rho = 7.4 * 10^(-7))
#     mle_hmm <- compute_rhat_hmm(frequencies, dataset_$dt, Ys_hmm, epsilon)
#     # compute second order derivative of -log-likelihood at MLE, divided by large sample size
#     FIM_iid <- numDeriv::hessian(function(x) -ll_iid(x), mle_iid)/large_sample_size
#     FIM_hmm <- numDeriv::hessian(function(x) -ll_hmm(x[1], x[2]), mle_hmm)/large_sample_size
#     var_hmm <- try(solve(FIM_hmm)) # Solve is matrix inversion
#     if (inherits(var_hmm, "try-error")){
#       var_hmm <- matrix(NA, 2, 2)
#     }
#     c(solve(FIM_iid)[1,1], var_hmm[2,2])
#   }
#   ## Comment: we can get NA when trying to invert the FIM
#   ## especially when the MLE is at the boundary of the parameter space, i.e rhat = 0 or = 1.
#   meanasymptvariances <- colMeans(asymptvariances, na.rm = TRUE)
#   for (irow in 1:length(samplesizes)){
#     sdrhats_iid <- sqrt(meanasymptvariances[1]/samplesizes[irow])
#     sdrhats_hmm <- sqrt(meanasymptvariances[2]/samplesizes[irow])
#     table_asymptotics_m_iid[irow, icolumn] <- 1.96 * sdrhats_iid * 2 # Final 2 ensure full CI width
#     table_asymptotics_m_hmm[irow, icolumn] <- 1.96 * sdrhats_hmm * 2 # Final 2 ensure full CI width
#   }
# }
# save(table_asymptotics_m_iid, table_asymptotics_m_hmm,
#      file = '~/Dropbox/IBD_IBS/RData/tables_asymptotics_m.RData' )
# load(file = '~/Dropbox/IBD_IBS/RData/tables_asymptotics_m.RData')
# table_asymptotics_m_iid
# tables_many_repeats_m_iid[,,1]

# # check visually for iid case
# par(mfrow = c(1,3))
# for (icol in 1:3){
#   plot(x = samplesizes, y = table_asymptotics_m_iid[,icol], type = "l", main = paste0("r = ", rgrid[icol]), ylim = c(0.1, 1))
#   lines(x = samplesizes, y = tables_many_repeats_m_iid[,icol,1], lty = 2)
# }

# check visually for hmm
# table_asymptotics_m_hmm
# tables_many_repeats_m_hmm[,,1]
# check visually for iid case
# par(mfrow = c(1,3))
# for (icol in 1:3){
#   plot(x = samplesizes, y = table_asymptotics_m_hmm[,icol], type = "l", main = paste0("r = ", rgrid[icol]), ylim = c(0.1, 2))
#   lines(x = samplesizes, y = tables_many_repeats_m_hmm[,icol,1], lty = 2)
# }


#====================================================
# Different (little) ks and ms given fixed r = 0.5 
# Use the same data set as in for 
#====================================================
fs_strategy <- list("Proportional to MAF" = maf)
tables_many_repeats_k_hmm <- array(dim = c(length(samplesizes), length(littleks), length(fs_strategy)),
                                   dimnames = list(samplesizes, littleks, names(fs_strategy)))
tables_many_repeats_k_iid <- array(dim = c(length(samplesizes), length(littleks), length(fs_strategy)),
                                   dimnames = list(samplesizes, littleks, names(fs_strategy)))

system.time(
  for(ilayer in 1:length(fs_strategy)){
    for (irow in 1:length(samplesizes)){
      # draw frequencies / positions 'once'
      dataset_name <- sprintf('m%s_fs:%s', samplesizes[irow], names(fs_strategy)[[ilayer]])
      dataset_ <- tables_many_repeats_m_dataset[[dataset_name]] 
      frequencies <- cbind(1-dataset_$fs, dataset_$fs)
      for (icolumn in 1:length(littleks)){
        # compute many MLEs under the hmm for these frequencies and positions
        rhats_ <- matrix(foreach (irepeat = 1:nrepeats, .combine = rbind) %dorng% {
          Ys <- simulate_Ys_hmm(frequencies, dataset_$dt, littleks[icolumn], rfixed, epsilon)
          c(compute_rhat_iid(frequencies, Ys, epsilon), compute_rhat_hmm(frequencies, dataset_$dt, Ys, epsilon))
        }, ncol = 3)
        tables_many_repeats_k_iid[irow, icolumn, ilayer] <- sqrt(mean((rhats_[,1] - rfixed)^2))
        tables_many_repeats_k_hmm[irow, icolumn, ilayer] <- sqrt(mean((rhats_[,3] - rfixed)^2))
      }
    }
  }
)
save(tables_many_repeats_k_iid, tables_many_repeats_k_hmm, file = '~/Dropbox/IBD_IBS/RData/tables_many_repeats_littlek.RData' )
load('../RData/tables_many_repeats_littlek.RData')
# show some of the results
cbind(tables_many_repeats_k_iid[,,1], tables_many_repeats_k_hmm[,,1])


#====================================================
# Different Ks and rs given fixed m = 96 
#====================================================
alphas = c(100, 1) # High versus low entropy
Ks = seq(2,20,2)
samplesize = 96 # Fix at 96 polyallelic markers

# Mechanism to generate f given K: Option 2
simulate_fs = function(K, alpha, samplesize){
  fs = rdirichlet(n = samplesize, rep(alpha, K))
  fs <- t(apply(fs, 1, function(v) (v + 1e-8)/sum(v + 1e-8)))
  return(fs)
}

tables_many_repeats_K_iid <- array(dim = c(length(Ks), length(rgrid), length(alphas)),
                                   dimnames = list(Ks, rgrid, alphas))
tables_many_repeats_K_hmm <- tables_many_repeats_K_iid
tables_many_repeats_K_dataset <- list()
dts <- Distances[['96']] # Fix the distances 
tables_many_repeats_K_dataset[['dts']] <- dts
system.time(
  for(ilayer in 1:length(alphas)){
    for (irow in 1:length(Ks)){
      # draw frequencies 'once'
      frequencies <- simulate_fs(Ks[irow], alphas[ilayer], samplesize)
      freq_name <- sprintf('K%s_alpha:%s', Ks[irow], alphas[[ilayer]])
      tables_many_repeats_K_dataset[[freq_name]] <- frequencies
      for (icolumn in 1:length(rgrid)){
        rhats_ <- matrix(foreach (irepeat = 1:nrepeats, .combine = rbind) %dorng% {
          Ys <- simulate_Ys_hmm(frequencies, dts, kfixed, rgrid[icolumn], epsilon)
          c(compute_rhat_iid(frequencies, Ys, epsilon), compute_rhat_hmm(frequencies, dts, Ys, epsilon))
        }, ncol = 3)
        tables_many_repeats_K_iid[irow, icolumn, ilayer] <- sqrt(mean((rhats_[,1] - rgrid[icolumn])^2))
        tables_many_repeats_K_hmm[irow, icolumn, ilayer] <- sqrt(mean((rhats_[,3] - rgrid[icolumn])^2))
      }
    }
  }
)
save(tables_many_repeats_K_dataset, tables_many_repeats_K_iid, tables_many_repeats_K_hmm,
     file = '../RData/tables_many_repeats_bigK.RData' )
load(file = '../RData/tables_many_repeats_bigK.RData')

tables_many_repeats_K_iid[,,2]
tables_many_repeats_K_hmm[,,2]


#====================================================
# Different ms and Ks given fixed r  
# addition of 6 and 12 resulted in errors. Not sure why
#====================================================
alphas = c(100, 1) # High versus low entropy
Ks = seq(2,20,4)

# Mechanism to generate f given K: Option 2
simulate_fs = function(K, alpha, samplesize){
  fs = rdirichlet(n = samplesize, rep(alpha, K))
  fs <- t(apply(fs, 1, function(v) (v + 1e-8)/sum(v + 1e-8)))
  return(fs)
}


tables_many_repeats_Km_iid <- array(dim = c(length(samplesizes), length(Ks), length(alphas)),
                                    dimnames = list(samplesizes, Ks, alphas))
tables_many_repeats_Km_hmm <- tables_many_repeats_Km_iid
tables_many_repeats_Km_dataset <- list()

system.time(
  for(samplesize in samplesizes){
    dts <- Distances[[as.character(samplesize)]] # Fix the distances 
    for(ilayer in 1:length(alphas)){
      for (icolumn in 1:length(Ks)){
        # draw frequencies 'once'
        frequencies <- simulate_fs(Ks[icolumn], alphas[ilayer], samplesize)
        freq_name <- sprintf('K%s_alpha:%s_%s', Ks[icolumn], alphas[[ilayer]], samplesize)
        tables_many_repeats_Km_dataset[[freq_name]] <- frequencies
        rhats_ <- matrix(foreach (irepeat = 1:nrepeats, .combine = rbind) %dorng% {
          Ys <- simulate_Ys_hmm(frequencies, dts, kfixed, rfixed, epsilon)
          c(compute_rhat_iid(frequencies, Ys, epsilon), compute_rhat_hmm(frequencies, dts, Ys, epsilon))
        }, ncol = 3)
        tables_many_repeats_Km_iid[as.character(samplesize), icolumn, ilayer] <- sqrt(mean((rhats_[,1] - rfixed)^2))
        tables_many_repeats_Km_hmm[as.character(samplesize), icolumn, ilayer] <- sqrt(mean((rhats_[,3] - rfixed)^2))
      }
    }
  }
)
save(tables_many_repeats_Km_dataset, tables_many_repeats_Km_iid, tables_many_repeats_Km_hmm,
     file = '../RData/tables_many_repeats_bigKm.RData' )
load(file = '../RData/tables_many_repeats_bigKm.RData')

tables_many_repeats_Km_iid[,,1]
tables_many_repeats_Km_hmm[,,1]








