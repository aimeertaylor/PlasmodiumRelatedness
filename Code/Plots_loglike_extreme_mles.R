#=====================================
# Understanding extreme values of r and k 
# and sensitivity to start values
#------------------------------------
# Summary: large k estimates due nearly flat likelihoods
# Baysian approach with informative priors would regularise
# Missing data may bias relatedness estimates towards one
#
# When r is close to either 0 or 1
# maximum of log-likelihood of r is indifferent to k
# and the log-likelihood of k is practically flat
# so we simply return the start value 
# (which was 50 for these mles)
# Except examples from "Rachuonyo"
# log-likelihood of r is not flat, and
# so not sensitive to starting values
#=====================================
rm(list = ls())
library(Rcpp)
sourceCpp("./hmmloglikelihood.cpp")
PLOT = T
if(PLOT){pdf("../Plots/Plots_loglike_extreme_mles.pdf")}

## Mechanism to compute MLE given fs, distances, Ys, epsilon (added str_vals)
str_vals <- c(1,0.5) # Starting values used in Generate_mles.R (we should change to 12 when re-run for bootstrap)
compute_rhat_hmm <- function(frequencies, distances, Ys, epsilon, str_vals){
  ndata <- nrow(frequencies)
  ll <- function(k, r) loglikelihood_cpp(k, r, Ys, frequencies, distances, epsilon, rho = 7.4 * 10^(-7))
  optimization <- optim(par = str_vals, fn = function(x) - ll(x[1], x[2]))
  rhat <- optimization$par
  return(rhat)
}

## Function for plotting log likelihood
ll_hmm <- function(r,k) loglikelihood_cpp(k, r, Ys, frequencies, subdata$dt, epsilon = 0.001, rho = 7.4 * 10^(-7))

## Data for plugging into likelihood functions
load('../RData/hmmInput_freqs.RData')

# Site summaries for selecting sites
load('../RData/nsnps_nsamps.RData') 
sites = names(which(nsnps_nsamps[2,] > 100)) # Select sites 

ks = c(1,10,100,1000)
rs = seq(0,0.999,length.out = 4)
ks_many = seq(1,1000,1)
rs_many <- seq(0,0.999,0.001) # X axis 
cols_ks = rainbow(length(ks), end = 0.8);names(cols_ks) = ks
cols_rs = rainbow(length(ks), end = 0.8);names(cols_rs) = rs


#=====================================================
# Extreme values of r 
#=====================================================
Input_sensitivity_rxtreme = list()

for(site in 'Colombia'){
  
  # Get data for plugging into the liklihood function 
  subdata <- hmmInput_freqs[[site]]
  
  # Get mles for selecting extremes
  load(sprintf('../RData/%s_mles.RData', site))
  frequencies <- cbind(1-subdata$fs, subdata$fs)
  subdata$dt <- c(diff(subdata$pos), Inf)
  pos_change_chrom <- 1 + which(diff(subdata$chrom) != 0) # find places where chromosome changes
  subdata$dt[pos_change_chrom-1] <- Inf
  
  # Obtain indices of extremes 
  inds = c('min r' = which(mle_df$rhat_hmm == min(mle_df$rhat_hmm))[1], 
           'max r' = which(mle_df$rhat_hmm == max(mle_df$rhat_hmm))[1], 
           'r approx 0.5' = which(mle_df$rhat_hmm > 0.48 & mle_df$rhat_hmm < 0.52)[1])
  any(mle_df$rhat_hmm == 0.5) # No evidence exact start being returned. 
  
  # Create stores to evaluate sensitivity of estimates to start values 
  k_input_sensitivity = array(dim = c(length(ks),length(inds)), 
                              dimnames = list(ks, names(inds)))
  r_input_sensitivity = array(dim = c(length(rs),length(inds)), 
                              dimnames = list(rs, names(inds)))
  
  par(mfrow = c(3,2), pty = 'm', family = 'serif', mar = c(3,6,1,1), oma = c(1,2,2,1))
  
  for(i in names(inds)){ 
    
    # Get data and original mle values
    
    ind = inds[i] # Get the index of a data comparison with extreme values 
    Ys <- as.matrix(subdata[,c(as.character(mle_df$individual1)[ind], as.character(mle_df$individual2)[ind])])
    mle_dfs <- as.matrix(mle_df[ind, c('rhat_hmm', 'khat_hmm')]) 
    m_ij = sum(rowSums(is.na(Ys)) == 0) # Number of non-missing SNP comparisons
    print(sprintf('%s %s %s', site, i, m_ij))
    
    #----------------------------------------------
    # Fix k and plot log likelihood of r
    #----------------------------------------------
    plot(y = sapply(rs_many, ll_hmm, k = mle_dfs[,'khat_hmm']), x = rs_many, type = 'l', xaxt = 'n', panel.first = grid(), 
         las = 2, ylab = '', xlab = '')
    axis(side = 1, at = c(0,1,0.5), tick = F, line = -1)
    abline(v = mle_dfs[,'rhat_hmm'], col = 'black', lty = 5)
    title(xlab = expression(italic(r)), line = 1)
    title(ylab = expression('log likelihood of'~italic(r)), line = 3.5)
    for(k in ks){ # Add likelihhoods for different k 
      lines(y = sapply(rs_many, ll_hmm, k), x = rs_many, col = cols_ks[as.character(k)])
      krhat_hmm_str <- compute_rhat_hmm(frequencies, subdata$dt, Ys, epsilon = 0.001, str_vals = c(k,0.5)) 
      k_input_sensitivity[as.character(k), i] <- krhat_hmm_str[1]
    }
    lines(y = sapply(rs_many, ll_hmm, k = mle_dfs[,'khat_hmm']), x = rs_many, col = 'black', lty = 5) # Re-add real
    legend(col = c(cols_ks, "#000000FF"), 'bottomright', inset = 0.07, bg = 'white',
           legend = format(c(names(cols_ks), round(mle_dfs[,'khat_hmm'])), digits = 0, drop0trailing = F), 
           cex = 0.75, pch = 20, title = expression(italic(k)))
    
    #----------------------------------------------
    # Fix r and plot log likelihood of k
    #----------------------------------------------
    plot(NULL, xaxt = 'n', panel.first = grid(), 
         las = 2, ylab = '', xlab = '',
         xlim = range(ks_many), 
         ylim = range(sapply(rs_many, ll_hmm, k = mle_dfs[,'khat_hmm'])))
    axis(side = 1, at = ks, tick = F, line = -1)
    abline(v = mle_dfs[,'khat_hmm'], col = 'black', lty = 5)
    title(xlab = expression(italic(k)), line = 1)
    title(ylab = expression('log likelihood of'~italic(k)), line = 3.5)
    for(r in rs){
      lines(y = sapply(ks_many, ll_hmm, r = r), x = ks_many, col = cols_rs[as.character(r)])
      krhat_hmm_str <- compute_rhat_hmm(frequencies, subdata$dt, Ys, epsilon = 0.001, str_vals = c(1,r))
      r_input_sensitivity[as.character(r), i] <- krhat_hmm_str[2]
    }
    lines(y = sapply(ks_many, ll_hmm, r = mle_dfs[,'rhat_hmm']), x = ks_many, col = 'black', lty = 5) # Re-add real
    legend(col = c(cols_rs, "#000000FF"), 'bottomright', inset = 0.07, bg = 'white',
           legend = format(c(as.numeric(names(cols_rs)), round(mle_dfs[,'rhat_hmm'],2)), digits = 2, drop0trailing = F), 
           cex = 0.75, pch = 20, title = expression(italic(r)))
  }
  
  # Add some text to plots
  mtext(text = site, outer = T, side = 3, line = 0)
  mtext(text = c(expression(italic(widehat(r))[italic(m)]%~~%0.5),
                 expression('maximum'~italic(widehat(r))[italic(m)]),
                 expression('minimum'~italic(widehat(r))[italic(m)])), outer = T, side = 2, line = 0, at = c(0.15,0.5,0.85))
  
  # Collate results in terms of sensitivity to starting values
  Input_sensitivity_rxtreme[[site]] = list(k_input_sensitivity = k_input_sensitivity, 
                                           r_input_sensitivity = r_input_sensitivity)
}


#=====================================================
# Extreme values of k 
#=====================================================
Input_sensitivity_kxtreme = list()

for(site in 'Colombia'){
  
  # Get data for plugging into the liklihood function 
  subdata <- hmmInput_freqs[[site]]
  
  # Get mles for selecting extremes
  load(sprintf('../RData/%s_mles.RData', site))
  frequencies <- cbind(1-subdata$fs, subdata$fs)
  subdata$dt <- c(diff(subdata$pos), Inf)
  pos_change_chrom <- 1 + which(diff(subdata$chrom) != 0) # find places where chromosome changes
  subdata$dt[pos_change_chrom-1] <- Inf
  
  # Obtain indices of extremes 
  inds = c('max k' = which(mle_df$khat_hmm == max(mle_df$khat_hmm))[1], 
           'min k' = which(mle_df$khat_hmm == min(mle_df$khat_hmm))[1])
  
  ks[which.max(ks)] = max(mle_df$khat_hmm)
  ks_many[which.max(ks_many)] = max(mle_df$khat_hmm)
  
  # Create stores to evaluate sensitivity of estimates to start values 
  k_input_sensitivity = array(dim = c(length(ks),length(inds)), 
                              dimnames = list(ks, names(inds)))
  r_input_sensitivity = array(dim = c(length(rs),length(inds)), 
                              dimnames = list(rs, names(inds)))
  
  par(mfrow = c(2,2), pty = 'm', family = 'serif', mar = c(3,6,1,1), oma = c(1,1,2,1))
  
  for(i in names(inds)){ 
    
    # Get data and original mle values
    ind = inds[i] # Get the index of a data comparison with extreme values 
    Ys <- as.matrix(subdata[,c(as.character(mle_df$individual1)[ind], as.character(mle_df$individual2)[ind])])
    mle_dfs <- as.matrix(mle_df[ind, c('rhat_hmm', 'khat_hmm')]) 
    
    #----------------------------------------------
    # Fix k and plot log likelihood of r
    #----------------------------------------------
    plot(y = sapply(rs_many, ll_hmm, k = mle_dfs[,'khat_hmm']), x = rs_many, type = 'l', xaxt = 'n', panel.first = grid(), 
         las = 2, ylab = '', xlab = '')
    axis(side = 1, at = c(0,1,0.5), tick = F, line = -1)
    abline(v = mle_dfs[,'rhat_hmm'], col = 'black', lty = 5)
    title(xlab = expression(italic(r)), line = 1)
    title(ylab = expression('log likelihood of'~italic(r)), line = 3.5)
    for(k in ks){ # Add likelihhoods for different k 
      lines(y = sapply(rs_many, ll_hmm, k), x = rs_many, col = cols_ks[as.character(k)])
      krhat_hmm_str <- compute_rhat_hmm(frequencies, subdata$dt, Ys, epsilon = 0.001, str_vals = c(k,0.5)) 
      k_input_sensitivity[as.character(k), i] <- krhat_hmm_str[1]
    }
    lines(y = sapply(rs_many, ll_hmm, k = mle_dfs[,'khat_hmm']), x = rs_many, col = 'black', lty = 5) # Re-add real
    legend(col = c(cols_ks, "#000000FF"), 'bottomright', inset = 0.07,
           legend = format(c(names(cols_ks), round(mle_dfs[,'khat_hmm'])), digits = 0, drop0trailing = F), 
           cex = 0.75, pch = 20, title = expression(italic(k)))
    
    #----------------------------------------------
    # Fix r and plot log likelihood of k
    #----------------------------------------------
    plot(NULL, xaxt = 'n', panel.first = grid(), 
         las = 2, ylab = '', xlab = '',
         xlim = range(ks_many), 
         ylim = range(sapply(rs_many, ll_hmm, k = mle_dfs[,'khat_hmm'])))
    axis(side = 1, at = ks, tick = F, line = -1)
    abline(v = mle_dfs[,'khat_hmm'], col = 'black', lty = 5)
    title(xlab = expression(italic(k)), line = 1)
    title(ylab = expression('log likelihood of'~italic(k)), line = 3.5)
    for(r in rs){
      lines(y = sapply(ks_many, ll_hmm, r = r), x = ks_many, col = cols_rs[as.character(r)])
      krhat_hmm_str <- compute_rhat_hmm(frequencies, subdata$dt, Ys, epsilon = 0.001, str_vals = c(1,r))
      r_input_sensitivity[as.character(r), i] <- krhat_hmm_str[2]
    }
    lines(y = sapply(ks_many, ll_hmm, r = mle_dfs[,'rhat_hmm']), x = ks_many, col = 'black', lty = 5) # Re-add real
    legend(col = c(cols_rs, "#000000FF"), 'bottomright', inset = 0.07, 
           legend = format(c(as.numeric(names(cols_rs)), round(mle_dfs[,'rhat_hmm'],2)), digits = 2, drop0trailing = F), 
           cex = 0.75, pch = 20, title = expression(italic(r)))
  }
  
  # Add some text to plots
  mtext(text = site, outer = T, side = 3, line = 0)
  mtext(text = c(expression('minimum'~italic(k)), 
                 expression('maximum'~italic(k))), outer = T, side = 2, line = 0, at = c(0.25,0.75))
  
  # Collate results in terms of sensitivity to starting values
  Input_sensitivity_kxtreme[[site]] = list(k_input_sensitivity = k_input_sensitivity, 
                                           r_input_sensitivity = r_input_sensitivity)
}
if(PLOT){dev.off()}


#========================================================
# Summarise sensitivity to starting values
#--------------------------------------------------------
# For all sites besides "Rachuonyo", 
# r insensitive to start values
# For all sites, k sensitive to start values
#========================================================
sd_cutoff_k <- 10
sd_cutoff_r <- 0.1 

# For min r, r approx 0.5, max r 
# Sensitivity of k to starting values 
cbind(t(sapply(Input_sensitivity_rxtreme, function(X){
  apply(X[['k_input_sensitivity']],2,sd) > sd_cutoff_k  
})),t(sapply(Input_sensitivity_kxtreme, function(X){
  apply(X[['k_input_sensitivity']],2,sd) > sd_cutoff_r  
})))

# Sensitivity of r to starting values 
cbind(t(sapply(Input_sensitivity_rxtreme, function(X){
  apply(X[['r_input_sensitivity']],2,sd) > sd_cutoff_r  
})), t(sapply(Input_sensitivity_kxtreme, function(X){
  apply(X[['r_input_sensitivity']],2,sd) > sd_cutoff_r  
})))




#========================================================
# Exploring Thai WGS comparison whose r approx. 0.5 but 
# k > 600 (outlier identified in Plot_real_mles.R), 
# compared with a comparison with similar r but typical k. 
# These plots do not feature in the main manuscript
#--------------------------------------------------------
# In conclusion: larger than average k is compatible with 
# fragmented IBS pattern of the outlier. However its 
# excessively large absolute size appears to be due to 
# numerical optimisation of a very slowly increasing
# likelihood. As such, a MAP estimate given informative 
# prior would likely regularise such large k.  
#========================================================
if(PLOT){png('../Plots/TM_WGS_r0.5_outlier.png', 
             width = 2400, height = 3000, res = 300)} # plot as png as otherwise large
library(RColorBrewer) # For colours
library(plotrix) # For plotting data radially
cols = brewer.pal(3, 'Dark2')
subdata <- hmmInput_freqs[['TM_WGS']] # Get data 
frequencies <- cbind(1-subdata$fs, subdata$fs)
subdata$dt <- c(diff(subdata$pos), Inf)
pos_change_chrom <- 1 + which(diff(subdata$chrom) != 0) # find places where chromosome changes
subdata$dt[pos_change_chrom-1] <- Inf

load('../RData/TM_WGS_mles.RData') # get mles 
inds_r = mle_df$rhat_hmm > 0.4 & mle_df$rhat_hmm < 0.6 # inds based on r estimates
#plot(y = mle_df$khat_hmm[inds_r], x = mle_df$rhat_hmm[inds_r]) # plot to identify anomalie
ind_rk = c(which(inds_r & mle_df$khat_hmm < median(mle_df$khat_hmm[inds_r]))[1], 
           which(inds_r & mle_df$khat_hmm > 600)[1])# inds to pick out anomalie

#--------------------------------------------
# Plot log-likelihoods
#--------------------------------------------
par(mfcol = c(3,2), pty = 'm', family = 'serif', mar = c(3,6,1,1), oma = c(1,3,2,1))
for(ind in ind_rk){
  
  # Get data and mles for specific example
  Ys <- as.matrix(subdata[,c(as.character(mle_df$individual1)[ind], 
                             as.character(mle_df$individual2)[ind])])
  mle_dfs <- as.matrix(mle_df[ind, c('rhat_hmm', 'khat_hmm')]) 
  
  #----------------------------------------------
  # First plot the data directly
  #----------------------------------------------
  IBS = 1*(Ys[,1] == Ys[,2]) # Code as identical by state or not
  cexs = c('0' = 0.2, '1' = 0.05) # To see 0s better
  offsets = c('0' = 10, '1' = -10) # To see 0s better
  plt.lns <- seq(1, length(IBS), 1)
  angles <- seq(0, 5*360, length=length(IBS))%%360
  polar.plot(plt.lns, polar.pos=angles+offsets[as.character(IBS)], 
             labels ='', show.grid.labels = FALSE, show.grid = F, line.col = NA, 
             rp.type = "symbols", point.symbols=20, point.col = cols[IBS+1], cex = cexs[as.character(IBS)])
  legend('topright', legend = c(sprintf('0 (%s pct)', round(mean(IBS == 0, na.rm = T),2)),
                                sprintf('1 (%s pct)', round(mean(IBS == 1, na.rm = T),2)),
                                sprintf('NA (%s pct)', round(mean(is.na(IBS)),2))), 
         col = c(cols[1:2],NA), pch = 20, bty = 'n', pt.cex = c(0.2*10,0.1*10,NA), 
         cex = 0.75, title = expression('IBS'[t]))
  
  #----------------------------------------------
  # Fix k and plot log likelihood of r
  #----------------------------------------------
  plot(y = sapply(rs_many, ll_hmm, k = mle_dfs[,'khat_hmm']), x = rs_many, type = 'l', xaxt = 'n', panel.first = grid(), 
       las = 2, ylab = '', xlab = '')
  axis(side = 1, at = c(0,1,0.5), tick = F, line = -1)
  abline(v = mle_dfs[,'rhat_hmm'], col = 'black', lty = 5)
  title(xlab = expression(italic(r)), line = 1)
  title(ylab = expression('log likelihood of'~italic(r)), line = 3.5)
  lines(y = sapply(rs_many, ll_hmm, k = mle_dfs[,'khat_hmm']), x = rs_many, col = 'black', lty = 5) # Re-add real
  legend(col = "#000000FF", 'bottomright', inset = 0.07,
         legend = format(mle_dfs[,'khat_hmm'], digits = 0, drop0trailing = F), 
         cex = 0.75, pch = 20, title = expression(italic(k)))
  
  #----------------------------------------------
  # Fix r and plot log likelihood of k
  #----------------------------------------------
  plot(NULL, xaxt = 'n', panel.first = grid(), 
       las = 2, ylab = '', xlab = '',
       xlim = range(ks_many), 
       ylim = range(sapply(rs_many, ll_hmm, k = mle_dfs[,'khat_hmm'])))
  axis(side = 1, at = ks, tick = F, line = -1)
  abline(v = mle_dfs[,'khat_hmm'], col = 'black', lty = 5)
  title(xlab = expression(italic(k)), line = 1)
  title(ylab = expression('log likelihood of'~italic(k)), line = 3.5)
  lines(y = sapply(ks_many, ll_hmm, r = mle_dfs[,'rhat_hmm']), x = ks_many, col = 'black', lty = 5) # Re-add real
  legend(col = "#000000FF", 'bottomright', inset = 0.07, 
         legend = format(mle_dfs[,'rhat_hmm'], digits = 2, drop0trailing = F), 
         cex = 0.75, pch = 20, title = expression(italic(r)))
}
mtext(text = c('Representative comparison','Outlier comparison'),outer = T, side = 3, line = 0, at = c(0.25, 0.75))
mtext(text = c('log liklihood plots', 'Data plot'), outer = T, side = 2, line = 1, at = c(0.40, 0.80))
if(PLOT){dev.off()}


#=========================================================
# Are k and r estimates correlated with NA values? 
#---------------------------------------------------------
# No striking patterns in k to worry about besides Rachuonyo 
# Relatedness estimates biased towards one when data missing  
#=========================================================
par(mfrow = c(3,2))
for(site in sites){
  
  # Get data for plugging into the liklihood function
  load(sprintf('../RData/%s_mles.RData', site))
  subdata <- hmmInput_freqs[[site]]
  inds = sample(size = 1000, x = nrow(mle_df)) # sample a collection of values 
  
  # Calculate average missingness
  m_succ = sapply(inds, function(ind){
      Ys <- as.matrix(subdata[,c(as.character(mle_df$individual1)[ind], 
                                 as.character(mle_df$individual2)[ind])])
      nrow(Ys) - sum(rowSums(is.na(Ys)) > 1)})
  
  # Plot missingness agains extreme k values
  plot(y = mle_df$khat_hmm[inds], x = m_succ, main = site, xlim = c(1,max(m_succ)), 
       pch = 20, col = adjustcolor('black',alpha.f = 0.5), 
       ylab = expression(hat(italic(k))), xlab = 'Number of loci with data')
  
  # Plot missingness agains extreme r values
  plot(y = mle_df$rhat_hmm[inds], x = m_succ, main = site, xlim = c(1,max(m_succ)), 
       pch = 20, col = adjustcolor('black',alpha.f = 0.5), 
       ylab = expression(hat(italic(r))), xlab = 'Number of loci with data')
}

#=========================================================
# What happens when no to increasing data etc 
# Depends on fs and whether IBSt = 1 or 1 for t = 1:m 
# If IBSt = 1 and m_ij is small -> bias towards r = 1
# Instability in k estimates is not a problem of small m_ij
#=========================================================
set.seed(1)
subdata <- hmmInput_freqs[['TM_WGS']] # Get data 
frequencies <- cbind(1-subdata$fs, subdata$fs)
subdata$dt <- c(diff(subdata$pos), Inf)
pos_change_chrom <- 1 + which(diff(subdata$chrom) != 0) # find places where chromosome changes
subdata$dt[pos_change_chrom-1] <- Inf

m <- 96
Ys <- cbind(rep(NA,m),rep(NA,m)) # No data
for(i in 1:m){
  Ys[i,] <- sample(c(0,1), 2, replace = T, prob = frequencies[i,]) # r should converge to 0
  compute_rhat_hmm(frequencies, subdata$dt, Ys, epsilon = 0.001, str_vals = c(50,0.5))
}













