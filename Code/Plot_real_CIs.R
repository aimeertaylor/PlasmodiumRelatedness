##################################################
# Plotting rs with CIs 
##################################################

# Set up 
rm(list = ls())
set.seed(6)
load('../RData/nsnps_nsamps.RData') # For plots of real CIs
load('../RData/allele_frequencies.RData')
sites <- colnames(nsnps_nsamps)
PDF = T
Vivax = T

Names = c(Colombia = 'Colombia', 
          TM_Barcode = 'Thailand 93-SNP',
          TM_WGS = 'Thailand WGS', 
          WKenya = 'Western Kenya', 
          Gambia = 'The Gambia', 
          Kilifi = 'Kilifi',
          vivax = 'Thailand MS')

# Ordering data sets by info content proxy 
Keffs = sapply(allele_frequencies, function(x){
  Kteffs = 1/(x^2 + (1-x)^2)
  mean(Kteffs)})
mKeffs = sort(Keffs[sites]*nsnps_nsamps[1,sites])
sites = names(mKeffs) # Sites reordered according to mKeff

if(Vivax){
  load("../RData/MS_vivax.RData")
  Keff_vivax = mean(1/rowSums(MS_freq^2))
  mKeff_vivax = nrow(MS_freq) * Keff_vivax
}


# Plot non vivax function
plot_non_vivax = function(site){
  
  load(sprintf('../RData/%s_mles_CI.RData', site))
  Order = sort.int(mle_df_CI$rhat_hmm, index.return = T)$ix
  num_plot = length(Order)
  
  plot(NULL, ylim = c(0,1), xlim = c(1,num_plot),
       bty = 'n', yaxt = 'n',
       panel.first = grid(), main = Names[site],
       ylab = expression('Relatedness'~italic(r)), 
       xlab = expression('Index as ranked by'~italic(r)))
  axis(side = 2, las = 2)
  title(main = bquote(italic(m)[max]%*%bar(italic(K))~"'"[italic(m)[max]]==.(round(mKeffs[site]))), line = 0.2)
  
  # Plot polygones
  polygon(c(1:num_plot, rev(1:num_plot)),
          y = c(mle_df_CI$CI_low[Order],rev(mle_df_CI$CI_high[Order])),
          border = NA, col = adjustcolor('gray',alpha.f = 0.5))
  
  # # Plot segments
  # segments(x0 = 1:num_plot, x1 = 1:num_plot,  
  #          y0 = mle_df_CI$CI_low[Order], y1 = mle_df_CI$CI_high[Order], 
  #          col = 'darkgray')
  
  # Plot MLEs
  points(mle_df_CI$rhat_hmm[Order], pch = 20, cex = 0.2)
}




#====================================================
# Plotting rs
#====================================================
if(PDF){pdf(file = '../Plots/Plot_real_CIs.pdf', height = 8, width = 9)}
par(mfrow = c(3,3), mar = c(5,5,4,1), family = 'serif')

# Plot rhats with CIs
if(!Vivax){
  for(site in sites){plot_non_vivax(site)}
} else {
  for(site in names(which(mKeffs < mKeff_vivax))){plot_non_vivax(site)}
  
  #-----------------------
  # Add vivax
  #-----------------------
  load('../RData/mle_df_vivax.RData') # to change 
  num_plot = 100
  
  # Set 100 with probability r
  inds = sample(1:nrow(mle_df), num_plot, prob = mle_df$rhat_iid)
  mle_df = mle_df[inds,] # Filter
  Order = sort.int(mle_df$rhat_iid, index.return = T)$ix
  
  plot(NULL, ylim = c(0,1), xlim = c(1,num_plot),
       bty = 'n', yaxt = 'n',
       panel.first = grid(), main = 'Thailand MS', 
       ylab = expression('Relatedness'~italic(r)), 
       xlab = expression('Index as ranked by'~italic(r)))
  axis(side = 2, las = 2)
  title(main = bquote(italic(m)[max]%*%bar(italic(K))~"'"[italic(m)[max]]==.(round(mKeff_vivax))), line = 0.3)

  # Plot polygones
  polygon(c(1:num_plot, rev(1:num_plot)),
          y = c(mle_df$rhat_iid_CI_low[Order],rev(mle_df$rhat_iid_CI_high[Order])),
          border = NA, col = adjustcolor('gray',alpha.f = 0.5))
  
  # # Plot segments
  # segments(x0 = 1:num_plot, x1 = 1:num_plot,  
  #          y0 = mle_df$rhat_iid_CI_low[Order], y1 = mle_df$rhat_iid_CI_high[Order], 
  #          col = 'darkgray')
  points(mle_df$rhat_iid[Order], pch = 20, cex = 0.2) # Plot MLEs
  
  for(site in names(which(mKeffs > mKeff_vivax))){plot_non_vivax(site)}
}

if(PDF){dev.off()}







