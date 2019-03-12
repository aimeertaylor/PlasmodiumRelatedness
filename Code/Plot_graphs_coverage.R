#====================================================
# Plot MSE for different m and r
# Modifications: colour by width of by 
#====================================================
rm(list = ls())
dev.off()
library(RColorBrewer) # For color mapper functions
par_df <- par(no.readonly = T)
PDF = T
if(PDF){pdf(file = '../Plots/Plot_graphs_coverage.pdf', 
            height = 8, width = 7)}
cov_cutoff = 0.95

# Load data 
par(mfrow = c(2,2), family = 'serif')
load(file = '../RData/tables_many_repeats_Kr_lengthcoverage_bootstrap.RData')
load(file = '../RData/tables_many_repeats_Kk_lengthcoverage_bootstrap.RData')
load(file = '../RData/tables_many_repeats_mr_lengthcoverage_bootstrap.RData')
load(file = '../RData/tables_many_repeats_mk_lengthcoverage_bootstrap.RData')

#################################################
# Coverage plots
#################################################
# Function to plot graphs
Plot_coverage = function(XLAB, LEG_title, HMM = F, INDEP = F){
  
  # Extract info
  X0 = Tables_cov$hmm; X1 = Tables_cov$iid
  
  # Get variables
  Vs = as.numeric(rownames(X0)) # Xaxis vaiable
  vs = as.numeric(colnames(X1)) # Yaxis variable
  cols_vs = rainbow(length(vs), end = 0.8); names(cols_vs) = vs # Variable

  # Empty plot
  plot(NULL, xlim = range(Vs), ylim = c(0.4,1), panel.first = grid(), ylog = T, 
       xaxt = 'n', bty = 'n', xlab = XLAB, ylab = 'Coverage')
  legend('bottomleft', lwd = 2, col = cols_vs, 
         legend = format(vs, digits = 2, drop0trailing = F), 
         bty = 'n', title = LEG_title, cex = 1, inset = 0.05)
  axis(side = 1, Vs)
  abline(h = cov_cutoff, lty = 5)
  if(INDEP){Title <- "Independence model"}; if(HMM){Title <- "HMM"}
  title(main = Title)
  
  # Plot coverage
  for(j in 1:ncol(X0)){
    if(INDEP){lines(y = X1[,j], x = Vs, col = cols_vs[j], ylog = T, pch = 20, type = 'b')}
    if(HMM){lines(y = X0[,j], x = Vs, col = cols_vs[j], ylog = T, pch = 20, type = 'b')}
  }
}


#================================================
# Big K: affect of alpha when variable r
# Little difference so just use alpha 1 in ms
#================================================
Tables_len = list(hmm = tables_many_repeats_Kr_length_hmm[,,1],  
                  iid = tables_many_repeats_Kr_length_iid[,,1])
Tables_cov = list(hmm = tables_many_repeats_Kr_coverage_hmm[,,1], 
                  iid = tables_many_repeats_Kr_coverage_iid[,,1])
Plot_coverage(XLAB = expression(italic(K)), LEG_title = expression(italic(r)), HMM = T)
Plot_coverage(XLAB = expression(italic(K)), LEG_title = expression(italic(r)), INDEP = T)

Tables_len = list(hmm = tables_many_repeats_Kr_length_hmm[,,2],  
                  iid = tables_many_repeats_Kr_length_iid[,,2])
Tables_cov = list(hmm = tables_many_repeats_Kr_coverage_hmm[,,2], 
                  iid = tables_many_repeats_Kr_coverage_iid[,,2])
Plot_coverage(XLAB = expression(italic(K)), LEG_title = expression(italic(r)), HMM = T)
Plot_coverage(XLAB = expression(italic(K)), LEG_title = expression(italic(r)), INDEP = T)

#================================================
# Big K: affect of alpha when variable k
# Little difference so just use alpha 1 in ms
#================================================
Tables_len = list(hmm = tables_many_repeats_Kk_length_hmm[,,1],  
                  iid = tables_many_repeats_Kk_length_iid[,,1])
Tables_cov = list(hmm = tables_many_repeats_Kk_coverage_hmm[,,1], 
                  iid = tables_many_repeats_Kk_coverage_iid[,,1])
Plot_coverage(XLAB = expression(italic(K)), LEG_title = expression(italic(k)), HMM = T)
Plot_coverage(XLAB = expression(italic(K)), LEG_title = expression(italic(k)), INDEP = T)

Tables_len = list(hmm = tables_many_repeats_Kk_length_hmm[,,2],  
                  iid = tables_many_repeats_Kk_length_iid[,,2])
Tables_cov = list(hmm = tables_many_repeats_Kk_coverage_hmm[,,2], 
                  iid = tables_many_repeats_Kk_coverage_iid[,,2])
Plot_coverage(XLAB = expression(italic(K)), LEG_title = expression(italic(k)), HMM = T)
Plot_coverage(XLAB = expression(italic(K)), LEG_title = expression(italic(k)), INDEP = T)


#================================================
# Big K: summary plots
#================================================
Tables_len = list(hmm = tables_many_repeats_Kr_length_hmm[,,1],  
                  iid = tables_many_repeats_Kr_length_iid[,,1])
Tables_cov = list(hmm = tables_many_repeats_Kr_coverage_hmm[,,1], 
                  iid = tables_many_repeats_Kr_coverage_iid[,,1])
Plot_coverage(XLAB = expression(italic(K)), LEG_title = expression(italic(r)), HMM = T)
Plot_coverage(XLAB = expression(italic(K)), LEG_title = expression(italic(r)), INDEP = T)

Tables_len = list(hmm = tables_many_repeats_Kk_length_hmm[,,1],  
                  iid = tables_many_repeats_Kk_length_iid[,,1])
Tables_cov = list(hmm = tables_many_repeats_Kk_coverage_hmm[,,1], 
                  iid = tables_many_repeats_Kk_coverage_iid[,,1])
Plot_coverage(XLAB = expression(italic(K)), LEG_title = expression(italic(k)), HMM = T)
Plot_coverage(XLAB = expression(italic(K)), LEG_title = expression(italic(k)), INDEP = T)


#================================================
# m: summary plots
#================================================
Tables_len = list(hmm = tables_many_repeats_mr_length_hmm[,,'Proportional to MAF'],  
                  iid = tables_many_repeats_mr_length_iid[,,'Proportional to MAF'])
Tables_cov = list(hmm = tables_many_repeats_mr_coverage_hmm[,,'Proportional to MAF'], 
                  iid = tables_many_repeats_mr_coverage_iid[,,'Proportional to MAF'])
Plot_coverage(XLAB = expression('Number of markers,'~italic(m)), 
              LEG_title = expression('Relatedness,'~italic(r)), HMM = T)
Plot_coverage(XLAB = expression('Number of markers,'~italic(m)), 
              LEG_title = expression('Relatedness,'~italic(r)), INDEP = T)

Tables_len = list(hmm = tables_many_repeats_mk_length_hmm[,,'Proportional to MAF'],  
                  iid = tables_many_repeats_mk_length_iid[,,'Proportional to MAF'])
Tables_cov = list(hmm = tables_many_repeats_mk_coverage_hmm[,,'Proportional to MAF'], 
                  iid = tables_many_repeats_mk_coverage_iid[,,'Proportional to MAF'])
Plot_coverage(XLAB = expression('Number of markers,'~italic(m)), 
              LEG_title = expression('Switch rate parameter,'~italic(k)), HMM = T)
Plot_coverage(XLAB = expression('Number of markers,'~italic(m)), 
              LEG_title = expression('Switch rate parameter,'~italic(k)), INDEP = T)



#################################################
# CI length plots
#################################################
# Function to plot graphs
Plot_CI_length = function(XLAB, LEG_title, HMM = F, INDEP = F){
  
  # Get variables
  Vs = as.numeric(rownames(Tables_len$hmm)) # Xaxis vaiable
  vs = as.numeric(colnames(Tables_len$hmm)) # Yaxis variable
  cols_vs = rainbow(length(vs), end = 0.8); names(cols_vs) = vs # Variable
  
  # Extract info
  X0 = Tables_len$hmm
  Z0 = Tables_cov$hmm
  X1 = Tables_len$iid
  Z1 = Tables_cov$iid
  
  # Empty plot
  plot(NULL, xlim = range(Vs), ylim = c(0,1), panel.first = grid(), ylog = T, 
       xaxt = 'n', bty = 'n', xlab = XLAB, ylab = 'Confidence interval width')
  legend('topright', lwd = 2, col = cols_vs, legend = vs, bty = 'n', title = LEG_title, cex = 0.5)
  axis(side = 1, Vs)
  if(INDEP){Title <- "Independence model"}; if(HMM){Title <- "HMM"}
  title(main = Title)
  
  # Plot CI lengths
  for(j in 1:ncol(X0)){
    if(INDEP){lines(y = X1[,j], x = Vs, col = cols_vs[j], ylog = T, pch = 20, type = 'b')}
    if(HMM){lines(y = X0[,j], x = Vs, col = cols_vs[j], ylog = T, pch = 20, type = 'b')}
  }
  
  # Add stars if coverage cutoff met
  if(INDEP){
    for(i in 1:nrow(X1)){
      inds = Z1[i,] > cov_cutoff
      points(y = max(X1[i,])+seq(0.12,0.2,length.out = length(vs))[inds],
             x = rep(Vs[i], sum(inds)), col = cols_vs[inds], pch = '*')
    }
  }
  
  # Add stars if coverage cutoff met
  if(HMM){
    for(i in 1:nrow(X0)){
      inds = Z0[i,] > cov_cutoff
      points(y = max(X0[i,])+seq(0.12,0.2,length.out = length(vs))[inds],
             x = rep(Vs[i], sum(inds)), col = cols_vs[inds], pch = '*')
    }
  }
  
}

#================================================
# Big K: affect of alpha when variable r
# Little difference so just use alpha 1 in ms
#================================================
Tables_len = list(hmm = tables_many_repeats_Kr_length_hmm[,,1],  
                  iid = tables_many_repeats_Kr_length_iid[,,1])
Tables_cov = list(hmm = tables_many_repeats_Kr_coverage_hmm[,,1], 
                  iid = tables_many_repeats_Kr_coverage_iid[,,1])
Plot_CI_length(XLAB = expression(italic(K)), LEG_title = expression(italic(r)), HMM = T)
Plot_CI_length(XLAB = expression(italic(K)), LEG_title = expression(italic(r)), INDEP = T)

Tables_len = list(hmm = tables_many_repeats_Kr_length_hmm[,,2],  
                  iid = tables_many_repeats_Kr_length_iid[,,2])
Tables_cov = list(hmm = tables_many_repeats_Kr_coverage_hmm[,,2], 
                  iid = tables_many_repeats_Kr_coverage_iid[,,2])
Plot_CI_length(XLAB = expression(italic(K)), LEG_title = expression(italic(k)), HMM = T)
Plot_CI_length(XLAB = expression(italic(K)), LEG_title = expression(italic(k)), INDEP = T)

#================================================
# Big K: affect of alpha when variable k
# Little difference so just use alpha 1 in ms
#================================================
Tables_len = list(hmm = tables_many_repeats_Kk_length_hmm[,,1],  
                  iid = tables_many_repeats_Kk_length_iid[,,1])
Tables_cov = list(hmm = tables_many_repeats_Kk_coverage_hmm[,,1], 
                  iid = tables_many_repeats_Kk_coverage_iid[,,1])
Plot_CI_length(XLAB = expression(italic(K)), LEG_title = expression(italic(k)), HMM = T)
Plot_CI_length(XLAB = expression(italic(K)), LEG_title = expression(italic(k)), INDEP = T)

Tables_len = list(hmm = tables_many_repeats_Kk_length_hmm[,,2],  
                  iid = tables_many_repeats_Kk_length_iid[,,2])
Tables_cov = list(hmm = tables_many_repeats_Kk_coverage_hmm[,,2], 
                  iid = tables_many_repeats_Kk_coverage_iid[,,2])
Plot_CI_length(XLAB = expression(italic(K)), LEG_title = expression(italic(k)), HMM = T)
Plot_CI_length(XLAB = expression(italic(K)), LEG_title = expression(italic(k)), INDEP = T)


#================================================
# Big K: summary plots
#================================================
Tables_len = list(hmm = tables_many_repeats_Kr_length_hmm[,,1],  
                  iid = tables_many_repeats_Kr_length_iid[,,1])
Tables_cov = list(hmm = tables_many_repeats_Kr_coverage_hmm[,,1], 
                  iid = tables_many_repeats_Kr_coverage_iid[,,1])
Plot_CI_length(XLAB = expression(italic(K)), LEG_title = expression(italic(r)), HMM = T)
Plot_CI_length(XLAB = expression(italic(K)), LEG_title = expression(italic(r)), INDEP = T)

Tables_len = list(hmm = tables_many_repeats_Kk_length_hmm[,,1],  
                  iid = tables_many_repeats_Kk_length_iid[,,1])
Tables_cov = list(hmm = tables_many_repeats_Kk_coverage_hmm[,,1], 
                  iid = tables_many_repeats_Kk_coverage_iid[,,1])
Plot_CI_length(XLAB = expression(italic(K)), LEG_title = expression(italic(k)), HMM = T)
Plot_CI_length(XLAB = expression(italic(K)), LEG_title = expression(italic(k)), INDEP = T)


#================================================
# m: summary plots
#================================================
Tables_len = list(hmm = tables_many_repeats_mr_length_hmm[,,'Proportional to MAF'],  
                  iid = tables_many_repeats_mr_length_iid[,,'Proportional to MAF'])
Tables_cov = list(hmm = tables_many_repeats_mr_coverage_hmm[,,'Proportional to MAF'], 
                  iid = tables_many_repeats_mr_coverage_iid[,,'Proportional to MAF'])
Plot_CI_length(XLAB = expression(italic(m)), LEG_title = expression(italic(r)), HMM = T)
Plot_CI_length(XLAB = expression(italic(m)), LEG_title = expression(italic(r)), INDEP = T)

Tables_len = list(hmm = tables_many_repeats_mk_length_hmm[,,'Proportional to MAF'],  
                  iid = tables_many_repeats_mk_length_iid[,,'Proportional to MAF'])
Tables_cov = list(hmm = tables_many_repeats_mk_coverage_hmm[,,'Proportional to MAF'], 
                  iid = tables_many_repeats_mk_coverage_iid[,,'Proportional to MAF'])
Plot_CI_length(XLAB = expression(italic(m)), LEG_title = expression(italic(k)), HMM = T)
Plot_CI_length(XLAB = expression(italic(m)), LEG_title = expression(italic(k)), INDEP = T)


if(PDF){dev.off()}
