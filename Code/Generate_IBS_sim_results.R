rm(list = ls())
set.seed(1)
library(dplyr)
library(Rcpp)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
source("./simulate_data.R")
epsilon = 0 # For IBS simple model 
nsim = 1000
ms = c(24, 96, 192) 
rs = c(0, 0.5)
corrected = c(F, T)

# Function to calculate IBS
calculate_IBShat <- function(x,y){
  z <- mean(x == y, na.rm = TRUE)
  return(z)
}

## Mechanism to generate Ys given fs, r, epsilon
simulate_Ys_iid <- function(frequencies, r, epsilon){
  Ys <- simulate_data(frequencies, rep(Inf, nrow(frequencies)), 1, r, epsilon, rho = 7.4 * 10^(-7))
  return(Ys)
}

# Load frequencies and positions from the Thai WGS data set 
load("../RData/hmmInput_freqs.RData")
data_ <- hmmInput_freqs$TM_WGS[,c(1,2,3)] # Using Thai frequencies and positions as a template
nSNPs <- nrow(data_)
mafs <- apply(cbind(data_$fs,1-data_$fs), 1, min)
data_$dt <- c(diff(data_$pos), Inf)
pos_change_chrom <- 1 + which(diff(data_$chrom) != 0) # find places where chromosome changes
data_$dt[pos_change_chrom-1] <- Inf


## Sample frequencies and distances (don't need the distances)
tables_ibs_dataset <- lapply(ms, function(m){
  distances = sample(data_$dt, m)
  frequencies = sample(data_$fs, m, prob = mafs)
  data.frame(distances = distances, frequencies = frequencies)
})
names(tables_ibs_dataset) = ms

# Generate stores 
tables_ibs <- list()

# IBS_simulation: convergence to expectation
for(m in ms){
  
  X = tables_ibs_dataset[[as.character(m)]]
  frequencies = cbind(X$frequencies, 1-X$frequencies) 
  h_constant = mean(rowSums(frequencies^2)) # Expected heterozygosity for outbred diploid
  
  for(r in rs){
    ibs_ibsc = foreach(irepeat = 1:nsim, .combine = rbind) %dorng% {
      Ys <- simulate_Ys_iid(frequencies, r, epsilon)    
      ibs <- calculate_IBShat(Ys[,1],Ys[,2])
      ibsc <- max(0, (ibs - h_constant)/(1 - h_constant)) # truncate
      data.frame(ibs = ibs, ibsc = ibsc)
    }
    tables_ibs[[as.character(m)]][[as.character(r)]] = ibs_ibsc
  }}


save(tables_ibs, tables_ibs_dataset, 
     file = '../RData/tables_ibs_results.RData')



