
rm(list = ls())
load(file = '../RData/tables_ibs_results.RData')
ms = names(tables_ibs)
rs = names(tables_ibs[[ms[1]]])
PDF = T
if(PDF){pdf('../Plots/Plot_ibs_results.pdf', 
            width = 10, height = 7)}

#================
# IBS
#================
par(mfcol = c(length(rs),length(ms)), family = 'serif',  mar = c(4,4,4,2), pty = 's')
for(m in ms){
  
  X = tables_ibs_dataset[[m]]
  frequencies = cbind(X$frequencies, 1-X$frequencies) 
  h_constant = mean(rowSums(frequencies^2)) # Expected heterozygosity for outbred diploid
  
  for(r in rs){
    Z = tables_ibs[[m]][[r]]$ibs 
    hist(Z, las = 1, panel.first = grid(), freq = F, 
         col = 'gray', xlim = c(0,1), 
         cex.axis = 0.75, xlab = bquote(widehat(IBS)[italic(m) == .(as.numeric(m))]),
         main = '', ylab = '')
    title(ylab = 'Density', line = 2)
    title(main = bquote(italic(r)==.(as.numeric(r))), line = 0)
    expectation = h_constant + (1-h_constant)*as.numeric(r)
    abline(v = expectation, lwd = 2, col = 'green')
    legend('topleft', fill = c('green'), cex = 0.75, inset = 0.075,
           legend = bquote(italic(bar(h)[m])+(1-italic(bar(h)[m]))*italic(r) == .(round(expectation,2))))
  }
}

#================
# IBS corrected 
#================
par(mfcol = c(length(rs),length(ms)), family = 'serif',  mar = c(4,4,2,2), pty = 's')
for(m in ms){
  for(r in rs){
    Z = tables_ibs[[m]][[r]]$ibsc 
    hist(Z, las = 1, panel.first = grid(), freq = F, 
         col = 'gray', xlim = c(0,1), cex.axis = 0.75,
         xlab = bquote(widehat(IBS)[italic(m) == .(as.numeric(m))]^(c)),
         main = '', ylab = '')
    title(ylab = 'Density', line = 2)
    title(main = bquote(italic(r)==.(as.numeric(r))), line = 0)
  }
}

#================
# IBS vs IBSc theoretical 
#================
par(family = 'serif', mar = c(5,5,1,1), mfrow = c(1,1))
hs = 1-seq(0.05, 0.45, length.out = 9)
cols = brewer.pal(length(hs), 'Spectral')
names(cols) = hs
X = seq(0,1,0.01)
plot(y = (X-hs[1])/(1-hs[1]), x = X, type = 'l', panel.first = grid(),  
     ylim = c(-.2,1), xaxt = 'n', 
     ylab = expression(widehat(IBS)[m]^(c)), bty = 'n', 
     xlab = expression(widehat(IBS)[m]), main = '', col  = cols[as.character(hs[1])])
for(h in hs[-1]){
  lines(y = (X-h)/(1-h), x = X, col = cols[as.character(h)])
}
legend('left', title = expression(italic(bar(h)[m])), 
       legend = format(hs, digits = 2, drop0trailing = F), col = cols, lwd = 2, bty = 'n')
axis(side = 1, at = seq(0,1,0.25), pos = c(0,0))



if(PDF){dev.off()}

  
  

