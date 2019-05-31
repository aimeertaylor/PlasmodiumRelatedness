########################################################
# Generate IBS for falciparum data sets 
########################################################
rm(list = ls())
load('../RData/hmmInput_freqs.RData')
load('../RData/allele_frequencies.RData')
load('../RData/nsnps_nsamps.RData')
MinSampleSize = 60 # Threshold on number of samples may change
sites = names(which(nsnps_nsamps[2,] > MinSampleSize)) # Extract sites to plot

# h constant and effective cardinality 
hs = sapply(allele_frequencies, function(x){mean(x^2 + (1-x)^2, na.rm = T)})

# Function to calculate IBS
calculate_IBShat <- function(x,y){
  z <- mean(x == y, na.rm = TRUE)
  return(z)
}

# Function to calculate IBSc
calculate_IBSc <- function(ibs, h_constant){
  z = max(0, (ibs - h_constant)/(1 - h_constant))
  return(z)
}

for (site in sites){
  print(site)
  data_set <- hmmInput_freqs[[site]]
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
  
  ibs_s <- foreach(icombination = 1:nrow(name_combinations), .combine = rbind) %dorng% {
    individual1 <- name_combinations[icombination,1]
    individual2 <- name_combinations[icombination,2]
    # let's focus on one pair of individuals
    i1 <- which(individual1 == names(data_set))
    i2 <- which(individual2 == names(data_set))
    ibs <- calculate_IBShat(data_set[,i1], data_set[,i2])
    ibsc <- calculate_IBSc(ibs, h_constant = hs[site])
    data.frame(ibs = ibs, ibsc = ibsc)
  }
  save(ibs_s, file = sprintf('../RData/%s_ibs.RData', site))
}
    


