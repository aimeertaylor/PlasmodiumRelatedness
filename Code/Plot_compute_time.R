### This script produces tables and plots of compute time
## Set up
rm(list = ls())
PDF = T
load('../RData/Comp_results.RData')
Mean_results = t(sapply(Comp_results, function(x) x[['Times']]))
rownames(Mean_results) = names(Comp_results)
f_strategy = "Proportional to MAF" # Choose a strategy
inds_i = grepl(f_strategy, rownames(Mean_results))

# Separate models
inds_iid =  grepl('iid', colnames(Mean_results))
inds_hmm =  grepl('hmm', colnames(Mean_results))

# Reformat rownames 
X_iid = Mean_results[inds_i,inds_iid]
rownames(X_iid) = gsub('m','',do.call(rbind, strsplit(rownames(X_iid), split = '_'))[,1])

# Reformat rownames 
X_hmm = Mean_results[inds_i,inds_hmm]
rownames(X_hmm) = gsub('m','',do.call(rbind, strsplit(rownames(X_hmm), split = '_'))[,1])

# Print latex code
kable(format(X_iid,digits = 4,drop0trailing = F), format = 'latex')
kable(format(X_hmm,digits = 4,drop0trailing = F), format = 'latex')

# Generate a plot
X = as.numeric(rownames(X_iid))

# Plot independence for self
plot(NULL, xlim = range(X), ylim = range(X_iid), 
     ylab = 'User CPU time (seconds)', 
     xlab = expression('Number of markers,'~italic(m)), 
     las = 1, xaxt = 'n')
for(j in 1:(ncol(X_iid)-1)){ # minus one to remove function calls 
  lines(x = X, y = X_iid[,j], type = 'b', col = j)} 
legend('top', legend = colnames(X_iid), fill = 1:(ncol(X_iid)-1))

# Plot HMM for self
plot(NULL, xlim = range(X), ylim = range(X_hmm[,1:3]))
for(j in 1:3){lines(x = X, y = X_hmm[,j], type = 'b', col = j)} 
legend('top', legend = colnames(X_hmm[,1:3]), fill = 1:3)

#=============================================
# Plots for manuscript
#=============================================
if(PDF){pdf('../Plots/Plot_compute_time.pdf', height = 6, width = 5)}
par(family = 'serif')
plot(NULL,ylim = range(0, 0.005), 
     xlim = range(X), 
     ylab = 'User CPU time (seconds)', 
     xlab = expression('Number of markers,'~italic(m)), 
     las = 1, xaxt = 'n',bty = 'n', panel.first = grid(nx=NA,ny=NULL))
axis(side = 1, at = X)
abline(v = X, col = "lightgray", lty = "dotted")
lines(x = X, y = X_hmm[,'hmm_user.self'], type = 'b', pch = 16)
lines(x = X, y = X_iid[,'iid_optimize_user.self'], type = 'b', pch = 21, bg = 'white', lty = 'dashed', col = 'blue')
lines(x = X, y = X_iid[,'iid_optim_user.self'], type = 'b', pch = 21, bg = 'white', lty = 'dashed')
legend(x = 24, y = 0.005, bty = 'n', lty = 1:2, pch = c(16,21), col = 'black', pt.bg = 'white', 
       legend = c('HMM', 'Independence model'), inset = 0.2)
legend(x = 24, y = 0.004, bty = 'n', fill = c('black', 'blue'),  
       legend = c('optim(...)', 'optimize(...)'), inset = 0.2)


# Summaries relation for caption
signif(summary(lm(X_hmm[,'hmm_user.self'] ~ X))$coefficients[,1], 3)
signif(summary(lm(X_iid[,'iid_optimize_user.self'] ~ X))$coefficients[,1], 3)
signif(summary(lm(X_iid[,'iid_optim_user.self'] ~ X))$coefficients[,1], 3)

# Plot HMM function calls
plot(NULL, xlim = range(X), ylim = c(1,max(X_hmm[,'hmm_function'], X_iid[,'iid_optim_function'])), 
     ylab = 'Number of calls to the likelihood function', 
     xlab = expression('Number of markers,'~italic(m)), 
     las = 1, xaxt = 'n', type = 'b', bty = 'n', 
     panel.first = grid(nx=NA,ny=NULL))
abline(v = X, col = "lightgray", lty = "dotted")
axis(side = 1, at = X)
lines(x = X, y = X_hmm[,'hmm_function'], type = 'b', pch = 16)
lines(x = X, y = X_iid[,'iid_optim_function'], type = 'b', pch = 21, bg = 'white', lty = 'dashed')

legend(x = 192, y = 100, bty = 'n', lty = 1:2, pch = c(16,21), col = 'black', pt.bg = 'white', 
       legend = c('HMM', 'Independence model'), inset = 0.2)

if(PDF){dev.off()}

# Convergence codes (inspect manually)
convergence_iid = sapply(Comp_results, function(x) x[[2]])
convergence_hmm = sapply(Comp_results, function(x) x[[3]])
print(convergence_hmm) # All 0
