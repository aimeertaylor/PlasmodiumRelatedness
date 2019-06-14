
rm(list = ls())
require(RColorBrewer)
JPEG=T
par_df <- par(no.readonly = T)
par(par_df); par(mfrow = c(1,1), family = 'serif')
if(JPEG){jpeg('../Plots/Plots_m_Keff%d.jpeg', height = 18, width = 18, res = 500, units = 'cm')}

#====================================================
# K vs Keff
#====================================================
load('../RData/tables_many_repeats_bigKm.RData') # Load fixed r results and data sets 
X <- dimnames(tables_many_repeats_Km_iid)
alphas <- array(as.numeric(X[[3]]),dim = length(X[[3]]), dimnames = list(X[[3]]))

# Remove distances to create a list of frequency datasets
Ks_datasets <- tables_many_repeats_Km_dataset
Ks_dataset_names <- names(Ks_datasets)

# Calculate the effective cardinalities
Keff_ms <- sapply(Ks_datasets, function(x){
  Keff_ts = 1/rowSums(x^2)
  mean(Keff_ts)
})
names(Keff_ms) <- Ks_dataset_names

# Extract other variables 
Ks = as.numeric(dimnames(tables_many_repeats_Km_hmm)[[2]])
ms = as.numeric(dimnames(tables_many_repeats_Km_hmm)[[1]])
cols = rainbow(length(ms), end = 0.8); names(cols) = ms

par(par_df); par(mfrow = c(1,2), family = 'serif')
indalpha <- rbind(grepl(alphas[1], Ks_dataset_names), !grepl(alphas[1], Ks_dataset_names))

for(i in 1:2){
  ind <- indalpha[i,]
  plot(NULL, xlim = range(Ks), ylim = range(Ks), xlab = 'K', ylab = 'K effective', 
       xaxt = 'n', yaxt = 'n', main = ifelse(i == 1, 'dist. of K close to uniform', 'Skewed'))
  axis(side = 1, at = Ks); axis(side = 2, at = Ks)
  abline(h = Ks,v = Ks, lty = 'dotted')
  for(K in Ks){
    indK <- grepl(sprintf('K%s_',K), names(Keff_ms))
    Ys <- Keff_ms[indK & ind]
    points(y = Ys, rep(K, length(Ys)), col = cols)
  }
  legend('topleft', col = cols, pch = 1, legend = ms, 
         inset = 0.05, title = expression(italic(m)))
}



#====================================================
# m vs Keff
#====================================================
par(mfrow = c(1,1), mar = c(5,5,1,1))
# Since the number of markers has no bearing, we first take an average over m 
newcols <- array(NA, dim = c(length(Ks),2), dimnames = list(Ks,NULL))
for(K in as.character(Ks)){
  newcols[K,1] = mean(Keff_ms[paste(sprintf('K%s_alpha:%s_',K,names(alphas)[1]),ms,sep = '')])
  newcols[K,2] = mean(Keff_ms[paste(sprintf('K%s_alpha:%s_',K,names(alphas)[2]),ms,sep = '')])
}

# Separated uniform and skewed and re-name columns with mean Keff
table1 = tables_many_repeats_Km_hmm[,,1]
table2 = tables_many_repeats_Km_hmm[,,2]
colnames(table1) <- newcols[,1]
colnames(table2) <- newcols[,2]

# Bring together into a single table and order the columns. 
tableboth = cbind(table1,table2)
order = sort.int(as.numeric(colnames(tableboth)),index.return = T)
tablebothsorted <- tableboth[,order$ix]

Z <- order$x
cols <- rainbow(length(Z), end = 0.8)
names(cols) <- Z

# HMM
plot(NULL, bty = 'n', cex.lab = 1.5, 
     ylim = c(0.05, pmax(max(tablebothsorted), 0.25)),
     xlim = range(ms), xaxt = 'n', panel.first = grid(),
     ylab = expression('Root mean squared error given data generating'~italic(r)==0.5), 
     main = '', las = 2, cex.lab = 1.5, 
     xlab = expression('Number of markers,'~italic(m)))
axis(side = 1, at = ms)
abline(h = 0.1, lty = 'dashed')

legend(y = 0.222, x = 384, cex = 1.2, 
       inset = 0, legend = format(Z, digits = 2, drop0trailing = F), 
       col = cols,
       pch = 16, bty = 'n', title = expression(bar(italic(K))~"'"[italic(m)[cum]]))

for (z in as.character(Z)){
  lines(x = ms, y = tablebothsorted[,z], type = "b", pch = 16, col = cols[z])
  lines(x = ms, y = tablebothsorted[,z], type = 'b', pch = 16, col = cols[z])
}


# Unpackage ms and multiply by Keff_ms
m_vector = sapply(1:length(Keff_ms), function(x){
  as.numeric(strsplit(names(Keff_ms)[x], '_')[[1]][3])})
names(m_vector) <- names(Keff_ms)
Keff_ms_m <- m_vector * Keff_ms

# Match the unpacking from as.vector below
Z = dimnames(tables_many_repeats_Km_hmm)
order_vector = rbind(rep(Z[[1]], length(Z[[2]])), rep(Z[[2]], each = length(Z[[1]])), 
                     rep(Z[[3]], each = length(Z[[1]])*length(Z[[2]])))
names_order_vector = apply(order_vector, 2, function(x){sprintf('K%s_alpha:%s_%s',x[2],x[3],x[1])})

# Collate and sort
Y <- as.vector(tables_many_repeats_Km_hmm) # Unpacks column by column 
X <- Keff_ms_m[names_order_vector]
X_sort <- sort.int(X, index.return = T)
Y_sort <- Y[X_sort$ix]

# Colour by Keff 
cols = rainbow(length(Keff_ms),end = 0.8)
names(cols) = names(sort(Keff_ms))
cols = sapply(cols, adjustcolor, alpha.f = 0.65)

# Plots 
symbols(y = Y, x = X, las = 1, 
        circles = sqrt(m_vector[names_order_vector]/pi), # Area scales with m
        inches = 0.15, 
        ylim = c(0.05, pmax(max(tablebothsorted), 0.25)), 
        xlab = expression('Number of markers times average effective cardinality, '~italic(m)%*%bar(italic(K))~"'"[italic(m)[cum]]), 
        ylab = expression('Root mean squared error given data generating'~italic(r)==0.5), 
        cex.lab = 1.5, 
        bty ='n', panel.first = grid(), 
        bg = cols[names(Keff_ms[names_order_vector])])
fit = loess(Y_sort ~ X_sort$x, span=0.20)
#lines(y = predict(fit), x = X_sort$x)
legend(y = 0.222, cex = 1.2, 
       x = 6418,
       fill = cols[seq(1,length(Keff_ms),length.out = 10)], bty = 'n', 
       title = expression(bar(italic(K))~"'"[italic(m)[cum]]), 
       legend = format(sort(Keff_ms)[seq(1,length(Keff_ms),length.out = 10)], digits = 2, drop0trailing = F))
# Legend for symbol size
symbols(y = seq(0.206, 0.145, length.out = length(ms)),
        x = rep(max(X),length(ms)), bg = 'lightgray', 
        circles = sqrt(ms/pi), 
        add = T, inches = 0.155) # Area scales with m
text(y = seq(0.206, 0.145, length.out = length(ms)),cex = 1.2,
     x = rep(max(X),length(ms)),labels = ms, pos = 4, offset = 1.1)
text(y = 0.206+0.01, pos = 4, x = max(X), labels = expression(italic(m)), cex = 1.2)
abline(h = 0.1, lty = 'dashed')

if(JPEG){dev.off()}








