##############################################################################################
# Freely available data from The Gambia and Kenya
# Omedo I, Mogeni P, Bousema T, Rockett K, Amambua-Ngwa A, Oyier I, et al. 
# Micro-epidemiological structuring of Plasmodium falciparum parasite populations in regions with 
# varying transmission intensities in Africa. Wellcome Open Res. 2017;2.

# SNPs typed, all from same umbrella panel (see manuscript)
# Kilifi, 158 to 226 SNPs per sample
# The Gambia, 131 SNPs in 143 samples 
# Rachuonyo south 111 in 2744 samples
#
# Exclusion criteria for both samples and SNPs of Omedo et al.:
# Samples with >40% of SNP failed excluded,
# Among included samples, SNPs with >30% of samples failed excluded
# SNPs: "Alleles were called as successful if they were above an intensity 
# cut-off value ranging between 0.5 and 1.0"
# Multiallelic SNPs: "For those SNPs that were above the cut-off, 
# allelic intensity ratios ranging between 0 and 1 were used to 
# classify them as homozygous or heterozygous."
# 
# Mixed SNP approach of Omedo et al.: "Where mixed parasite populations were identified, 
# we took the majority SNP calls at each position to indicate the dominant genotype"
#
#---------------------------------------------------------------------------------------------
# Seems some samples have been tested twice at half the SNPs rather than at all SNPs
# I don't think they were all typed at the same SNPs. 
#---------------------------------------------------------------------------------------------
##############################################################################################
rm(list = ls())

#---------------------------------------------------------------------------------------------
# Dataset 1: Genotyping results for 276 single nucleotide polymorphisms (SNPs) typed 
# in 5199 Plasmodium falciparum samples collected from The Gambia and Kenya.
#---------------------------------------------------------------------------------------------
SNPDataRaw <- read.delim('~/Documents/BroadLaptop/Kenya/OriginalData/FreelyAvailableOnline/TheGambiaKenya/Dataset1_InformationongenesandSNPstypedinP.falciparumparasitepopulations.txt')
head(SNPDataRaw)
str(SNPDataRaw) # 810330 obs. of  17 variables

# 1 Check that all Xs are failures: Yes
all(SNPDataRaw$pass_fail[SNPDataRaw$result == 'X'] == 0) 
# 2 Check that all failures are Xs: Yes
all(SNPDataRaw$result[SNPDataRaw$pass_fail == 0] == 'X') 
# 3 Check for pass/fail for NAs: no NAs
unique(SNPDataRaw$pass_fail)
# 4 Check result for NAs: no NAs, just Xs
unique(SNPDataRaw$result)

Sample_names <- as.character(unique(SNPDataRaw$sample_id))
SNP_names <- as.character(unique(SNPDataRaw$assay_code)) 
nsamp <- length(Sample_names) # 5026 samples (5199 in paper)
nSNPs <- length(SNP_names) # 272 SNPs (276 SNPs in paper)

# Break TheGambiaKenya down into different sites
sites = c('Gambia', 'Kilifi', 'Rachuonyo')
samples_gambia <- as.character(unique(SNPDataRaw[SNPDataRaw$study_location == 'gambia', 'sample_id']))
samples_kilifi <- as.character(unique(SNPDataRaw[SNPDataRaw$study_location == 'kilifi', 'sample_id']))
samples_rachuonyo <- as.character(unique(SNPDataRaw[SNPDataRaw$study_location == 'rachuonyo', 'sample_id']))
sample_names <- list(Gambia = samples_gambia, Kilifi = samples_kilifi, Rachuonyo = samples_rachuonyo)

Samples_with_duplicate_aft_exclusion_criteria = list()

for(site in sites){
  
  site_samples <- sample_names[[site]] # unlist(sample_names) # For all samples 
  SNPData <- array(dim = c(length(site_samples),nSNPs), 
                   dimnames = list(paste('sid', site_samples, sep = ''), SNP_names))
  Samples_with_repeat_SNPs <- NULL
  pb <- txtProgressBar(min = 0, max = length(site_samples), style = 3)
  
  # Unpackage the data (don't overwrite failures in order to distinguish failures from not typed)
  for(Sample_name in site_samples){ 
    setTxtProgressBar(pb, which(site_samples == Sample_name))
    ind <- Sample_name == as.character(SNPDataRaw$sample_id)
    Subset <- SNPDataRaw[ind, ]
    Row <- paste('sid', Sample_name, sep = '')
    Cols <-  as.character(Subset$assay_code)
    X <- table(Subset$assay_code) # Report any sample with two or more results per SNP
    if(any(X > 1)){Samples_with_repeat_SNPs <- c(Samples_with_repeat_SNPs, Row)}
    SNPData[Row, Cols] <- as.character(Subset$result) # remove later, keep for now
  }
  dim(SNPData) 
  
  # Remove SNPs that are NA across all samples (i.e. untyped at site)
  ind_snps_completely_na <- colMeans(is.na(SNPData)) == 1
  SNPData <- SNPData[, !ind_snps_completely_na]
  dim(SNPData)
  
  # Replace failed by NA 
  SNPData[SNPData == 'X'] <- NA 
  
  # Remove sample if > 40% missing (following Omedo2017a, Omedo2017b)
  ind_sample_40pct_plus_missing = rowMeans(is.na(SNPData)) > 0.4
  SNPData <- SNPData[!ind_sample_40pct_plus_missing, ]
  dim(SNPData)
  
  # Remove SNPs where >30% remaining samples failed (following Omedo2017a)
  ind_snp_30pct_plus_missing <- colMeans(is.na(SNPData)) > 0.3
  SNPData <- SNPData[,!ind_snp_30pct_plus_missing]
  dim(SNPData)

  # Of the samples remaining samples, do any have duplicate SNP calls? If so, remove
  ind_duplicate = rownames(SNPData) %in% Samples_with_repeat_SNPs
  Samples_with_duplicate_aft_exclusion_criteria[[site]] = rownames(SNPData)[ind_duplicate]
  SNPData = SNPData[!ind_duplicate,]
  dim(SNPData)

  # Remove samples with > 5% het (matches Thai barcode multiclonal cutoff)
  Het_threshold = 0.05
  Het_calls <- t(apply(SNPData, 1, function(x){sapply(strsplit(x, split = ''), length)}))
  Het_calls[is.na(SNPData)] = NA # Make sure missing are recognised as such for het rate calc. 
  Het_rate <- rowMeans(Het_calls > 1, na.rm = T) # Number of het calls per sample
  hist(Het_rate); abline(v = mean(Het_rate), col = 'blue')
  SNPData = SNPData[Het_rate <= Het_threshold, ]
  max(rowSums(Het_calls[Het_rate <= Het_threshold,] > 1, na.rm = T)) # check max het call remaining
  dim(SNPData)
  
  # Replace remaining het calls by NA
  SNPData[Het_calls[Het_rate <= Het_threshold,] > 1] <- NA # replace het calls in the SNPData with NA
  unique(as.vector(SNPData)) # Check no mixed (either got rid of or replaced by missing)
  
  # Convert data into zeros and ones 
  SNPData01 <- apply(SNPData, 2, function(x){
    y <- unique(x)
    z <- y[!is.na(y)]
    if(length(z)){ # Added for checking allele frequency spectra
      for(i in 1:length(z)){
        x <- gsub(z[i], as.character(i-1), x) 
      } 
    }
    x <- as.numeric(x)})
  rownames(SNPData01) <- rownames(SNPData)
  dim(SNPData01)
  
  # #++++++++++++++
  # # Added to inspect spectrum
  # freq = colMeans(SNPData01, na.rm = T)
  # MAF = pmin(1-freq, freq)
  # hist(MAF, col = 'gray')
  # #++++++++++++++
  
  # Remove any SNPs that are monomorphic in the remaining data
  ind_mono = colMeans(SNPData01, na.rm = T) == 0 | colMeans(SNPData01, na.rm = T) == 1
  SNPData01 <- SNPData01[,!ind_mono]
  print(sprintf('%s has %s samples with %s SNPs', site, nrow(SNPData01), ncol(SNPData01)))
  
  # Unlike with Western Kenyan data, don't need SNPinfo since chr and pos in SNPDataRaw
  chrom <- SNPDataRaw$chr_valid
  pos <- SNPDataRaw$coord_valid
  names(chrom) <- as.character(SNPDataRaw$assay_code)
  names(pos) <- as.character(SNPDataRaw$assay_code)
  
  # HMM format 
  HMMInput <- data.frame(chrom = chrom[colnames(SNPData01)], 
                         pos = pos[colnames(SNPData01)], 
                         t(SNPData01))
  colnames(HMMInput) <- c("chrom", "pos", rownames(SNPData01))
  
  # Re-order the SNPs
  old_order <- colnames(SNPData01)
  new_order <- vector('list', length(sort(unique(chrom))))
  for(chr in sort(unique(HMMInput$chrom))){
    ind <- HMMInput$chrom == chr
    X <- sort(as.numeric(HMMInput$pos[ind]), ind = TRUE)
    new_order[[chr]] <- old_order[ind][X$ix]
  }
  HMMInput <- HMMInput[unlist(new_order),]
  
  # Save data
  write.table(HMMInput, file = sprintf('~/Dropbox/IBD_IBS/TxtData/hmmInput_%s.txt', site),
              quote = FALSE, row.names = FALSE, sep = '\t')
}

Samples_with_duplicate_aft_exclusion_criteria

     