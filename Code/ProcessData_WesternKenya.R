##############################################################################################
# Freely available data from western Kenya 
# Omedo,I. et al. (2017) Geographic-genetic analysis of Plasmodium falciparum parasite 
# populations from surveys of primary school children in Western Kenya. Wellcome Open Res.
# https://wellcomeopenresearch.org/articles/2-29/v2

# Text from article: 
# Our positive control criteria were to include samples where at least 60% of SNP typing was successful 
# and, among these, to include SNPs that were successfully typed in at least 60% of all samples. 
# The selection criterion for successful typing was based on individually defined SNP intensity 
# values (R) ranging from 0 to 1. SNPs with intensity values <0.1 were considered low quality 
# and were categorized as failed and excluded from further analyses. In addition, allelic 
# intensity ratios (θ) nearing 0 or 1 were used to classify SNP positions as homozygous, and 
#intensity ratios of intermediate values were used to classify SNPs as heterozygous, 
# representing mixed parasite populations in a single sample. Where mixed parasite populations were identified, 
# we took the dominant genotype forward for further analysis, as represented by the majority SNP calls. 
# Applying these inclusion criteria, we restricted our analyses to 83 SNPs and 1809 samples
##############################################################################################
rm(list = ls())

#---------------------------------------------------------------------------------------------
# Dataset 1: Genotyping results for 111 single nucleotide polymorphisms (SNPs) typed 
# in 2486 Plasmodium falciparum samples collected from primary school children during 
# a parasitological survey in western Kenya in 2009 and 2010.
#---------------------------------------------------------------------------------------------
SNPDataRaw <- read.csv('~/Documents/BroadLaptop/Kenya/OriginalData/FreelyAvailableOnline/SchoolChildrenKenya/Dataset1_Genotypingresultsfor111SNPsin2486samples.csv')
head(SNPDataRaw)
str(SNPDataRaw) # 275946 obs. of  11 variables

# 1 Check that all Xs are failures: Yes
all(SNPDataRaw$pass_fail[SNPDataRaw$result == 'X'] == 0) 
# 2 Check that all failures are Xs: Yes
all(SNPDataRaw$result[SNPDataRaw$pass_fail == 0] == 'X') 
# 3 Check for pass/fail for NAs: no NAs
unique(SNPDataRaw$pass_fail)
# 4 Check result for NAs: no NAs, just Xs
unique(SNPDataRaw$result)

Sample_names <- as.character(unique(SNPDataRaw$sample_id))
SNP_names <- as.character(unique(SNPDataRaw$assay_code)) # 111 SNPs
nsamp <- length(Sample_names) # 2486 as reported in the paper
nSNPs <- length(SNP_names) # 111 as reported in the paper

# Make a store of the data 
SNPData <- array(dim = c(nsamp,nSNPs), dimnames = list(paste('sid', Sample_names, sep = ''), SNP_names))
Samples_with_repeat_SNPs <- NULL # To check if any duplicate calls
pb <- txtProgressBar(min = 0, max = nsamp, style = 3)

# Unpackage the data 
for(Sample_name in Sample_names){
  setTxtProgressBar(pb, which(Sample_names == Sample_name))
  ind <- Sample_name == as.character(SNPDataRaw$sample_id)
  Subset <- SNPDataRaw[ind, ]
  Row <- paste('sid', Sample_name, sep = '')
  Cols <-  as.character(Subset$assay_code)
  X <- table(Subset$assay_code) # Report any sample with two or more results per SNP
  if(any(X > 1)){Samples_with_repeat_SNPs <- c(Samples_with_repeat_SNPs, Row)}
  SNPData[Row, Cols] <- as.character(Subset$result)
}
dim(SNPData)

# Check for samples with duplicate SNP calls: no samples
Samples_with_repeat_SNPs

# Check to see if any untyped SNPs: no
any(colMeans(is.na(SNPData)) > 0)

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

# Remove samples with > 5% het (matches Thai barcode multiclonal cutoff)
Het_calls <- t(apply(SNPData, 1, function(x){sapply(strsplit(x, split = ''), length)}))
Het_calls[is.na(SNPData)] = NA # Make sure missing are recognised as such for het rate calc. 
Het_rate <- rowMeans(Het_calls > 1, na.rm = T) # Number of het calls per sample
hist(Het_rate); abline(v = mean(Het_rate), col = 'blue')
SNPData = SNPData[Het_rate <= 0.05, ]
max(rowSums(Het_calls[Het_rate <= 0.05,] > 1, na.rm = T)) # check max het call remaining
dim(SNPData)

# Replace remaining het calls by NA
SNPData[Het_calls[Het_rate <= 0.05,] > 1] <- NA # replace het calls in the SNPData with NA
unique(as.vector(SNPData)) # Check no mixed (either got rid of or replaced by missing)

# Convert data into zeros and ones 
SNPData01 <- apply(SNPData, 2, function(x){
  y <- unique(x)
  z <- y[!is.na(y)]
  for(i in 1:length(z)){
    x <- gsub(z[i], as.character(i-1), x) 
  }
  x <- as.numeric(x)})
rownames(SNPData01) <- rownames(SNPData)
dim(SNPData01)

# Remove any SNPs that are monomorphic in the remaining data
ind_mono = colMeans(SNPData01, na.rm = T) == 0 | colMeans(SNPData01, na.rm = T) == 1
SNPData01 <- SNPData01[,!ind_mono]
print(sprintf('%s samples with %s SNPs', nrow(SNPData01), ncol(SNPData01)))

# Extract pos and chrom of 111 SNPs
# Load SNPInfo in order to extract pos and chrom
SNPInfo <- read.table('~/Documents/BroadLaptop/Kenya/OriginalData/FreelyAvailableOnline/SchoolChildrenKenya/Supplementary_Table_1-Sequenom_assay_design.txt', 
                      sep = '\t', header = TRUE)
str(SNPInfo) # 111 obs. of  18 variables 
any(!sum(SNP_names %in% as.character(SNPInfo$assay_code))) # Check SNPs match: yes!
chrom <- SNPInfo$chr_valid
pos <- SNPInfo$coord_valid
names(chrom) <- as.character(SNPInfo$assay_code)
names(pos) <- as.character(SNPInfo$assay_code)

# HMM format 
HMMInput <- data.frame(chrom = chrom[colnames(SNPData01)], 
                       pos = pos[colnames(SNPData01)], t(SNPData01))
colnames(HMMInput) <- c("chrom", "pos", rownames(SNPData))

# Re-order the SNPs
old_order <- colnames(SNPData01)
new_order <- vector('list', length(sort(unique(chrom))))
for(chr in sort(unique(HMMInput$chrom))){
  ind <- HMMInput$chrom == chr
  X <- sort(as.numeric(HMMInput$pos[ind]), ind = TRUE)
  new_order[[chr]] <- old_order[ind][X$ix]
}
HMMInput <- HMMInput[unlist(new_order),]

# Save Data
write.table(HMMInput, file = '~/Dropbox/IBD_IBS/TxtData/hmmInput_WKenya.txt',
            quote = FALSE, row.names = FALSE, sep = '\t')


