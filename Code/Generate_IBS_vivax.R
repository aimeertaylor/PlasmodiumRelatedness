# Generate IBS for vivas only 

rm(list = ls())
load("../RData/MS_vivax.RData")
h_m = mean(rowSums(MS_freq^2))
mKeff = mean(1/rowSums(MS_freq^2))

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

inds_no_missing = !apply(MS_rank, 1, function(x)any(is.na(x)))
data_set <- MS_rank[inds_no_missing, ]
individual_names <- as.character(rownames(data_set))
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
  i1 <- which(individual1 == individual_names)
  i2 <- which(individual2 == individual_names)
  ibs <- calculate_IBShat(data_set[i1,-1], data_set[i2,-1])
  ibsc <- calculate_IBSc(ibs, h_constant = h_m)
  data.frame(ibs = ibs, ibsc = ibsc)
}
save(ibs_s, file = '../vivax_ibs.RData')



