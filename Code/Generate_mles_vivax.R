# This script is adapted for IBD/IBS from the original script saved in the vivax project  
# It currently includes the method for formating the data, generating mles and CIs as well
rm(list = ls())
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
nboot = 500
  
## Mechanism to compute MLE given fs, Ys, epsilon
compute_rhat_iid <- function(frequencies, Ys, epsilon){
  ndata <- nrow(frequencies)
  distances <- rep(Inf, ndata)
  ll <- function(r) loglikelihood_cpp(1, r, Ys, frequencies, distances, epsilon, rho = 7.4 * 10^(-7))
  optimization <- optimize(f = function(x) - ll(x), interval = c(0, 1))
  rhat <- optimization$minimum
  return(rhat)
}

## Mechanism to generate Ys given fs, r, epsilon
simulate_Ys_iid <- function(frequencies, r, epsilon){
  Ys <- simulate_data(frequencies, rep(Inf, nrow(frequencies)), 1, r, epsilon, rho = 7.4 * 10^(-7))
  return(Ys)
}

#=========================================================
# Load data and format
#=========================================================
load('../Data_for_relatedness.RData')

# Create a matrix of frequencies
markers <- names(Fs_Combined)
max_cardin <- max(sapply(Fs_Combined, length)) # To set dimension of frequency matrix
frequencies = array(0, dim = c(length(markers), max_cardin), 
                    dimnames = list(markers, 1:max_cardin))
for(x in markers){fs <- Fs_Combined[[x]]; frequencies[x, names(fs)] <- fs}
MS_freq = frequencies
rowSums(frequencies) # Check
dimnames(frequencies) = NULL

# Recode MS repeat lengths as ranks
MS_data = monoclonal_data[, markers]
MS_rank = monoclonal_data[, markers]
for(i in 1:ncol(MS_data)){
  x = MS_data[,i]
  cardin = length(unique(x))
  X = sort(unique(x))
  for(j in 1:cardin){
    x[x==X[j]] <- j-1
  }
  MS_rank[,i] <- x
}

# Save some processed data for IBS file 
rownames(MS_rank) = monoclonal_data$Episode_Identifier
save(MS_rank, MS_freq, file = "../RData/MS_vivax.RData") 

# Get all unique individuals
individual_names <- unique(monoclonal_data$ID) 
nindividuals <- length(individual_names)
print(nindividuals)
name_combinations <- matrix(nrow = nindividuals*(nindividuals-1)/2, ncol = 2)
count <- 0
for (i in 1:(nindividuals-1)){
  for (j in (i+1):nindividuals){
    count <- count + 1
    name_combinations[count,1] <- individual_names[i]
    name_combinations[count,2] <- individual_names[j]
  }
}

system.time( # 2.860 for 1000; 13 for 5000
  mle_df <- foreach(icombination = 1:nrow(name_combinations),  
                    .combine = rbind) %dorng% {
                      # Individual names
                      individual1 <- name_combinations[icombination,1] 
                      individual2 <- name_combinations[icombination,2]
                      # Search individual names in matrix
                      ind_indivd1 <- grepl(individual1, monoclonal_data$ID)
                      ind_indivd2 <- grepl(individual2, monoclonal_data$ID)
                      # Choose one of there episodes at random
                      epi1 <- sample(monoclonal_data$Episode_Identifier[ind_indivd1], 1)
                      epi2 <- sample(monoclonal_data$Episode_Identifier[ind_indivd2], 1)
                      # Create logical vector for the chosen episodes
                      ind_epi1 <- monoclonal_data$Episode_Identifier == epi1
                      ind_epi2 <- monoclonal_data$Episode_Identifier == epi2
                      
                      # Extract MS data 
                      Ys = cbind(as.numeric(MS_rank[ind_epi1, markers]), 
                                 as.numeric(MS_rank[ind_epi2, markers]))
                      
                      # Record number of makers
                      num_ms_typed = sum(rowSums(is.na(Ys)) == 0)
                      
                      # Estimate relatedness
                      rhat_iid <- compute_rhat_iid(frequencies, Ys, epsilon)
       
                      # Construct parametric boot strap
                      rhats_boot = rep(NA, nboot)
                      for(b in 1:nboot){
                        Ys_boot <- simulate_Ys_iid(frequencies, rhat_iid, epsilon)  
                        rhats_boot[b] <- compute_rhat_iid(frequencies, Ys_boot, epsilon)  
                      }
                      if(any(is.na(rhats_boot))){stop('Bootstraps NA')}
                      CI <- as.numeric(quantile(rhats_boot, probs = c(0.025, 0.975)))
                      
                      # Store results
                      data.frame(individual1 = individual1, 
                                 individual2 = individual2, 
                                 epi1 =  epi1, # Record episode selected 
                                 epi2 =  epi2, # Record episode selected 
                                 rhat_iid = rhat_iid, 
                                 rhat_iid_CI_low = CI[1],
                                 rhat_iid_CI_high = CI[2],
                                 num_ms_typed = num_ms_typed)
                    })

save(mle_df, file = "../RData/mle_df_vivax.RData")

load("../RData/mle_df_vivax.RData")
