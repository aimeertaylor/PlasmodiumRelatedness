########################################################
# Plots of all Plasmodium data sets: 
# makers, their frequencies and the distances between
########################################################
rm(list = ls())
par(mfrow = c(3,2))
require(RColorBrewer)
require(dplyr)
cols = brewer.pal(n = 4, name = "Dark2")
PDF = T
min_sample_no = 0

load('../RData/hmmInput_freqs.RData')
load('../RData/allele_frequencies.RData')
load('../RData/nsnps_nsamps.RData') 

# To select sites with more than X samples
sites = names(which(nsnps_nsamps[2,] > min_sample_no)) # Select sites 

# Names 
Names = c(Colombia = 'Colombia', 
          TM_Barcode = 'Thai 93-SNP',
          TM_WGS = 'Thai WGS', 
          WKenya = 'Western Kenya', 
          Gambia = 'The Gambia', 
          Kilifi = 'Kilifi')

#=========================================
# Plot actual data - plot as png as large 
#=========================================

png(file = '../Plots/Data_plot.png', 
    width = 2100, height = 3000)
par(mfrow = c(3,2), family = 'serif', mar = c(7,5,5,2))
for(site in sites){
  SNPData = as.matrix(hmmInput_freqs[[site]])[,-(1:3)]
  num_SNPs = nrow(SNPData) # Over all samples
  min_SNPs = min(colSums(!is.na(SNPData)))
  avg_SNPs = mean(colSums(!is.na(SNPData)))
  order = names(sort(colMeans(is.na(SNPData))))
  image(SNPData[,order], ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', col = cols)
  title(ylab = sprintf('Sample ID (%s samples)', ncol(SNPData)),
        line = 1.5, 
        main = Names[site], cex.lab = 5, cex.main = 5)
  title(xlab = sprintf('SNP ID (SNPs analysed: total %s, min %s, mean %s)', 
                       num_SNPs, min_SNPs, round(avg_SNPs)), cex.lab = 5, line = 3)
}
dev.off()



if(PDF){pdf('../Plots/Plot_real_data.pdf', 
            width = 10, height = 7)}

#=========================================
# Plot frequencies 
#=========================================
par(mfcol = c(2,3), mar = c(5,5,4,2), family = 'serif')
for(site in sites){
  X = pmin(1-allele_frequencies[[site]], allele_frequencies[[site]])
  hist(X, col = 'gray', breaks = 30, ylab = '', freq = T, 
       xlab = '', 
       xlim = c(0,0.5), panel.first= grid(), 
       main = Names[site], 
       las = 1)
  title(ylab = 'Number of SNPs', line = 3)
  title(xlab = 'Minor allele frequency estimate', line = 2)
}

#=========================================
# Plot distances
#=========================================
par(mfcol = c(2,3), mar = c(5,5,4,2), family = 'serif')
for(site in sites){
  dataset_ = hmmInput_freqs[[site]]
  dataset_ <-  dataset_ %>% arrange(chrom, pos) 
  dataset_$dt <- c(diff( dataset_$pos), Inf)
  pos_change_chrom <- 1 + which(diff( dataset_$chrom) != 0) # find places where chromosome changes
  dataset_$dt[pos_change_chrom-1] <- Inf
  
  X = log10(dataset_$dt)
  hist(X, col = 'gray', breaks = 30, ylab = '', freq = F, 
       xlab = '', main = Names[site], panel.first= grid(), 
       las = 1)
  title(ylab = 'Density', line = 3)
  title(xlab = expression(log[10]~'distance'), line = 3)
}
if(PDF){dev.off()}