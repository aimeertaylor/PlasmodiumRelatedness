rm(list = ls())

path <- '~/Documents/BroadLaptop/'
Colombia_hmmInput <- read.delim(sprintf('%sColombianBarcode/TxtData/hmmInput.txt', path))
TM_Barcode_hmmInput <- read.delim(sprintf('%s/TM_border/QuantLocalPfConnIBD/TxtData/Barcode93.txt', path))
TM_WGS_hmmInput <- read.delim(sprintf('%s/TM_border/QuantLocalPfConnIBD/TxtData/WGS.txt', path))
WKenya_hmmInput <- read.delim('~/Dropbox/IBD_IBS/TxtData/hmmInput_WKenya.txt')
Gambia_hmmInput <- read.delim('~/Dropbox/IBD_IBS/TxtData/hmmInput_Gambia.txt')
Kilifi_hmmInput <- read.delim('~/Dropbox/IBD_IBS/TxtData/hmmInput_Kilifi.txt')

# Collate all hmmInput
hmmInput <- list(Colombia = Colombia_hmmInput,
                 TM_Barcode = TM_Barcode_hmmInput, 
                 TM_WGS = TM_WGS_hmmInput,
                 WKenya = WKenya_hmmInput, 
                 Gambia = Gambia_hmmInput, 
                 Kilifi = Kilifi_hmmInput)

# Calculate allele frequencies
allele_frequencies <- sapply(hmmInput, function(x){x[x == -1] <- NA;  apply(x[, -(1:2)], 1, mean, na.rm = TRUE)})
nsnps_nsamps <- sapply(hmmInput, function(x){dim(x[,-(1:2)])})

# Create a store of hmmInput and frequencies combined
hmmInput_freqs = sapply(hmmInput, function(x){
  x[x == -1] <- NA
  fs = apply(x[, -(1:2)], 1, mean, na.rm = TRUE)
  y = cbind(x[,1:2], fs, x[,-(1:2)])
})

#================================================================================================
# Save results
#================================================================================================
save(allele_frequencies, file = '../RData/allele_frequencies.RData')
save(nsnps_nsamps, file = '../RData/nsnps_nsamps.RData')
save(hmmInput_freqs, file = '../RData/hmmInput_freqs.RData')

hs = sapply(allele_frequencies, function(x){mean(x^2 + (1-x)^2, na.rm = T)})
keffs = sapply(allele_frequencies, function(x){mean(1/(x^2 + (1-x)^2))})
round(keffs, 2)
round(nsnps_nsamps)
round(hs,2)
