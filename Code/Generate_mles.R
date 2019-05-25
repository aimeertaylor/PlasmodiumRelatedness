##########################################################
# Generate MLEs for all pairs in published processed data 
##########################################################

rm(list = ls())
set.seed(1) # for reproducibility
library(dplyr) # for arrange()
library(Rcpp)
library(doParallel)
library(doRNG)
source("./simulate_data.R")
sourceCpp("./hmmloglikelihood.cpp")
registerDoParallel(cores = detectCores()-2)
epsilon <- 0.001 # Fix epsilon throughout
kinit = 50 # To initiate numerical optimization
rinit = 0.5 # To initiate numerical optimization
load("../RData/hmmInput_freqs.RData") # Data 
dataset_names <- names(hmmInput_freqs) # Dataset names


## Mechanism to compute MLE given fs, Ys, epsilon
compute_rhat_iid <- function(frequencies, Ys, epsilon){
  ndata <- nrow(frequencies)
  distances <- rep(Inf, ndata)
  ll <- function(r) loglikelihood_cpp(1, r, Ys, frequencies, distances, epsilon, rho = 7.4 * 10^(-7))
  optimization <- optimize(f = function(x) - ll(x), interval = c(0, 1))
  rhat <- optimization$minimum
  return(rhat)
}

## Mechanism to compute MLE given fs, Ys, epsilon and distances
compute_rhat_hmm <- function(frequencies, distances, Ys, epsilon){
  ndata <- nrow(frequencies)
  ll <- function(k, r) loglikelihood_cpp(k, r, Ys, frequencies, distances, epsilon, rho = 7.4 * 10^(-7))
  optimization <- optim(par = c(kinit, rinit), fn = function(x) - ll(x[1], x[2]))
  rhat <- optimization$par
  return(rhat)
}

#=====================================================
# Extract pairwise comparisons and then calculate MLEs
#=====================================================
for (idataset in 1:length(dataset_names)){
  
  data_set <- hmmInput_freqs[[dataset_names[idataset]]]
  print(paste(dataset_names[idataset], dim(data_set))) 
  if (dim(data_set)[2] < 100){next} # Skip datasets with fewer than 100 samples
  
  # Sort data by chromosome and position (in cases otherwise unsorted)
  data_set <- data_set %>% arrange(chrom, pos) 
  data_set$dt <- c(diff(data_set$pos), Inf)
  pos_change_chrom <- 1 + which(diff(data_set$chrom) != 0) # find places where chromosome changes
  data_set$dt[pos_change_chrom-1] <- Inf
  
  # Extract distances and frequencies since same for all pairs within a dataset
  distances = data_set$dt  
  frequencies <- cbind(1-data_set$fs, data_set$fs)

  # Create list of all pairwise comparisons
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
  
  
  # For each pair in parallel
  mle_df <- foreach(icombination = 1:nrow(name_combinations), .combine = rbind) %dorng% {
    
    # Extract data for pair of individuals
    individual1 <- name_combinations[icombination,1]
    individual2 <- name_combinations[icombination,2]
    i1 <- which(individual1 == names(data_set))
    i2 <- which(individual2 == names(data_set))
    subdata <- data_set[,c(i1,i2)]
    names(subdata) <- c("Yi", "Yj")
    
    # Compute MLEs
    krhat_hmm <- compute_rhat_hmm(frequencies, distances, cbind(subdata$Yi, subdata$Yj), epsilon)
    rhat_iid <- compute_rhat_iid(frequencies, cbind(subdata$Yi, subdata$Yj), epsilon)
    
    # Store MLES
    data.frame(individual1 = individual1, individual2 = individual2, 
               rhat_iid = rhat_iid, khat_hmm = krhat_hmm[1], rhat_hmm = krhat_hmm[2])
  }
  save(mle_df, file = paste0("../RData/", dataset_names[idataset], "_mles.RData"))
}

