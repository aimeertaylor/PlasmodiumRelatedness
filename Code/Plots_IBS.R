# Plots containing IBS (includes sensitivity simulation)

rm(list = ls())
library(RColorBrewer)
library(sm)
load('../RData/allele_frequencies.RData')
load('../RData/nsnps_nsamps.RData')
load('../RData/vivax_ibs.RData') # Add ibs for vivas
load("../RData/MS_vivax.RData") # For MS_freq
source('./vioplot2.R') # For modified violin plots
MinSampleSize = 60 # Threshold on number of samples may change
PDF = T
Vivax = T # Set to true to include vivax

# Sites 
Names = c(Colombia = 'Colombia', 
          TM_Barcode = 'Thailand 93-SNP',
          TM_WGS = 'Thailand WGS', 
          WKenya = 'Western Kenya', 
          Gambia = 'The Gambia', 
          Kilifi = 'Kilifi',
          vivax = 'Thailand MS')

# Extract sites to plot
if(!Vivax){
  sites = names(which(nsnps_nsamps[2,] > MinSampleSize))  
} else {
  sites = c(names(which(nsnps_nsamps[2,] > MinSampleSize)), 'vivax')
}

# For plotting
cols = brewer.pal(n = 4, name = "Dark2")
names(cols) = c('rhat_hmm', 'IBS', 'rhat_iid', 'IBSc')
cols_sites = brewer.pal(n = length(sites), name = "Dark2")
names(cols_sites) = sites
subplot_rows = 2
subplot_cols = length(sites) / subplot_rows

# h constant and effective cardinality 
hs = sapply(allele_frequencies, function(x){mean(x^2 + (1-x)^2, na.rm = T)})
if(Vivax){hs = c(hs, 'vivax' = mean(rowSums(MS_freq^2)))} 

if(PDF){pdf('../Plots/Plots_IBS.pdf', height = 8, width = 7)}

#================================================
# Histograms of IBSc real for the supplementary
#================================================
par(family = 'serif', mfrow = c(3,3), 
    mar = c(4,4,2,2), pty = 's')
for(site in sites){
  load(sprintf('../RData/%s_ibs.RData', site))
  hist(ibs_s$ibsc,breaks = 30, plot = TRUE, las = 2, xlim = c(0,1),
       col = cols["IBSc"], cex.axis = 0.75, las = 1, freq = F, 
       main = Names[site], 
       xlab = '', panel.first = grid(), 
       ylab = 'Density')
  title(line = 3, xlab = expression(widehat(IBS)[m]))
  box()
}

#================================================
# Violin / histogram plot of IBS and rhat real
#================================================
VIO = T
cols = brewer.pal(n = 4, name = "Dark2")
Order = sort.int(hs[sites], index.return = T)
par(mfrow = c(2,1), pty = 'm', mar = c(4,7.5,1,1), family = 'serif')

plot(NULL, xlim = c(0,1), ylim = c(1, length(sites)+1), panel.first = grid(), 
     ylab = '', xlab = '', yaxt = 'n')
title(xlab = expression(widehat(IBS)[italic(m)]), line = 1.2)
axis(side = 2, at = seq(1.5,length(sites)+1,1), Names[names(Order$x)], 
     las = 2, tick = F,line = 0)
count = 1


# Histogram instead of violin plot

for(site in names(Order$x)){
  
  load(sprintf('../RData/%s_ibs.RData', site))
  z <- ibs_s$ibs
  if(VIO){
    vioplot2(z, horizontal = T, h = 0.06,
             at = count, side = 'right', add = T, col = cols[2],
             own_line = Order$x[site])
    points(x = mean(z), y = count-0.15, pch = 17, cex = 0.75)
  } else {
    x <- hist(z, plot = F,breaks = 30)
    x$density = (x$density - min(x$density))/((max(x$density) - min(x$density))*1.8) + count
    ind = x$density != count
    eps = mean(diff(x$mids))/2.5
    rect(xleft = (x$mids-eps)[ind], 
         xright = (x$mids+eps)[ind], 
         ybottom = rep(count, sum(ind)), 
         ytop = x$density[ind], 
         col = cols[2], 
         border = cols[2])
    points(x = Order$x[site], y = rep(count-0.15,1), pch = c(17), cex = 0.75)
  }
  count = count + 1
}

plot(NULL, xlim = c(0,1), ylim = c(1, length(sites)+1), panel.first = grid(), 
     ylab = '', xlab = '', yaxt = 'n')
title(xlab = expression(widehat(italic(r))[italic(m)]), line = 1.2)
axis(side = 2, at = seq(1.5,length(sites)+1,1), Names[names(Order$x)], 
     las = 2, tick = F,line = 0)
count = 1

for(site in names(Order$x)){
  
  if(site == 'vivax'){load(sprintf('../RData/mle_df_vivax.RData', site))
  } else {
    load(sprintf('../RData/%s_mles.RData', site))
  }   
  
  if(site == 'vivax'){
    z <- mle_df$rhat_iid
  } else {
    z <- mle_df$rhat_hmm
  }
  
  # z = z[z< 0.5] Explore outliers
  
  if(VIO){
    vioplot2(x = z, horizontal = T, h = 0.06, # Surpress line
             at = count, side = 'right', add = T, col = cols[1], 
             colline = NA) 
  } else {
    # Histogram instead of violin plot
    x <- hist(z, plot = F,breaks = 30)
    x$density = (x$density - min(x$density))/((max(x$density) - min(x$density))*1.8) + count
    ind = x$density != count
    eps = mean(diff(x$mids))/2.5
    rect(xleft = (x$mids-eps)[ind], 
         xright = (x$mids+eps)[ind], 
         ybottom = rep(count, sum(ind)), 
         ytop = x$density[ind], 
         col = cols[1], 
         border = cols[1])
  }
  points(x = mean(z), y = count-0.15, pch = 17, cex = 0.75)
  count = count + 1
}


#=========================================
# Generate plot of IBS inferior to IBD
#=========================================
load('../RData/IBS_sensitivity_sim.RData')
Order = sort.int(sapply(Xs[[1]], function(x){x$h_constant}), index.return = T)
lower = min(sapply(Xs[[1]], function(x){min(x$relatedness_measures[,'ibs'])}))
rfxd = 0.5 # Data generating relatedness

# Extract h constants
h_constants_tab = sapply(Xs, function(x){sapply(x, function(y){round(y$h_constant, 2)})})

for(i in 1:length(Xs)){
  
  X = Xs[[i]]
  Order = sort.int(sapply(X, function(x){x$h_constant}), index.return = T)
  lower = min(sapply(X, function(x){min(x$relatedness_measures[,'ibs'])}))
  par(mfrow = c(2,1), pty = 'm', mar = c(4,7.5,1,1), family = 'serif')
  lower = min(sapply(X, function(x){min(x$relatedness_measures[,'ibs'])}))
  sites = names(Order$x)
  
  
  plot(NULL, xlim = c(0,1), ylim = c(0.5, length(sites)+1), panel.first = grid(), 
       ylab = '', xlab = '', yaxt = 'n')
  title(xlab = expression(widehat(IBS)[italic(m)]), line = 1.2)
  axis(side = 2, at = seq(1.5,length(sites)+1,1), Names[names(Order$x)], 
       las = 2, tick = F,line = 0)
  count = 1
  
  for(site in sites){
    
    z <- X[[site]]$relatedness_measures[,'ibs']
    
    if(VIO){
      vioplot2(z, horizontal = T, h = 0.06,
               at = count, side = 'right', add = T, col = cols[2],
               own_line = Order$x[site] + (1-Order$x[site])*rfxd)
      points(x = mean(z), y = count-0.15, pch = 17, cex = 0.75)
    } else {
      x <- hist(z, plot = F,breaks = 30)
      x$density = (x$density - min(x$density))/((max(x$density) - min(x$density))*1.8) + count
      ind = x$density != count
      eps = mean(diff(x$mids))/2.5
      rect(xleft = (x$mids-eps)[ind], 
           xright = (x$mids+eps)[ind], 
           ybottom = rep(count, sum(ind)), 
           ytop = x$density[ind], 
           col = cols[2], 
           border = cols[2])
      points(x =  Order$x[site] + (1-Order$x[site])*rfxd, y = count-0.1, pch = 17, cex = 0.75)
    }
    count = count + 1
  }
  
  
  plot(NULL, xlim = c(0,1), ylim = c(0.5, length(sites)+1), panel.first = grid(), 
       ylab = '', xlab = '', yaxt = 'n')
  title(xlab = expression(widehat(italic(r))[italic(m)]), line = 1.2)
  axis(side = 2, at = seq(1.5,length(sites)+1,1), Names[names(Order$x)], 
       las = 2, tick = F,line = 0)
  count = 1
  
  for(site in sites){
    
    z = X[[site]]$relatedness_measures[,'rhat']
    
    if(VIO){
      vioplot2(z, horizontal = T, h = 0.06, 
               at = count, side = 'right', add = T, col = cols[1],
               own_line = NULL, colMed = cols[1], drawRect = F)
    } else {
      # Histogram instead of violin plot
      x <- hist(z, plot = F,breaks = 30)
      x$density = (x$density - min(x$density))/((max(x$density) - min(x$density))*1.8) + count
      ind = x$density != count
      eps = mean(diff(x$mids))/2.5
      rect(xleft = (x$mids-eps)[ind], 
           xright = (x$mids+eps)[ind], 
           ybottom = rep(count, sum(ind)), 
           ytop = x$density[ind], 
           col = cols[1], 
           border = cols[1])
    }
    points(x = mean(z), y = count-0.15, pch = 17, cex = 0.75)
    count = count + 1
  }
}

if(PDF){dev.off()}






