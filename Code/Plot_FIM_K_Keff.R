##########################################################
# Exploring precision and statdard deviation for different K and r
#
# FIM is a complicated function of K and r: 
# As r tends to 0, precision tends towards a linearly increasing function of K (FIM(K) = K-1 when r = 0)
# As r tends to 1, precision tends towards a logarithmically increasing function of K, with increasing precision 
# Within a biologically reasonable range of allelic richness (inc. microsat) this leads to the following trends 
# as function of increasing K (in relative terms):
# very, very rapid increases to high precision for r close to 1 (the FIM is not defined for r = 1)
# slow increases to low precision for r = 0.5
# rapid increase to moderate precision with K for r = 0

# Considering standard deviation, this translates to 
# rapid decrease to moderate stdev for r = 0 (biggest reduction in sd)
# slow decrease to high stdev for r = 0.5 (moderate reduction in sd)
# rapid decrease to low stdev for r > 0.5 (amounting to little reduction in sd when r close to 0)
##########################################################
rm(list = ls())
require(plotly)
PLOT3d <- F # Set to false to avoid plotting 3d 
PLOThm <- F # Set to false to avoid plotting heat maps
JPEG <- T # Set to false to print to screen 

if(JPEG){jpeg(file = '../Plots/Plot_FIM_K_eff%d.jpeg', res = 500, 
            height = 20, width = 20, units = 'cm')}

FIM_K <- function(x, r){ # Function to calculate FIM assuming equifrequent Fs
  1/(1-r) + (x-1)^2/(x*(1+(x-1)*r)) - 1/(x*(1-r))}

FIM_fs = function(fs,r){ # Function to calculate FIM for any fs
  t1 = 1/(1-r)
  t2 = (fs * (1 - fs)^2) / (r + fs*(1-r))
  t3 = fs^2 / (1-r)
  return(t1 + sum(t2 - t3))
}

# A vast number of rs (for surfaces etc.)
rs_many = seq(0.1,0.9,0.01) 
cols_rs_many = rainbow(n = length(rs_many), end = 0.8)
names(cols_rs_many) = rs_many

# A regular number of rs 
rs = seq(0.1,0.9,0.1) 
cols_rs = rainbow(n = length(rs), end = 0.8) # Same mapping 
names(cols_rs) = rs

# Descrete number of Ks
Ks = seq(2,15,1) 

# Calculate a FIM Surface using many rs
FIM_surface = sapply(rs_many, function(r){
  sapply(Ks, function(x) FIM_K(x, r))
})
dimnames(FIM_surface) = list(Ks,rs_many)


#=======================================
# FIM and stdev as func. of K and r
#=======================================
par(family = 'serif', mfrow = c(1,2))

plot(NULL, # FIM 
     xlim = range(Ks), ylim = c(0,100), bty = 'n',
     ylab = expression('FIM('~italic(K)~','~italic(r)~')'), 
     xlab = expression(italic(K)))
for(r in rev(rs_many)){
  lines(y = FIM_surface[,as.character(r)], x = Ks, col = cols_rs_many[as.character(r)])
}
plot(NULL, # Std dev. 
     xlim = c(2,20), ylim = c(0,1), bty = 'n',
     ylab = expression(sigma~'('~italic(K)~','~italic(r)~')'), 
     xlab = expression(italic(K)))
for(r in rev(rs_many)){
  lines(y = sqrt(1/FIM_surface[,as.character(r)]), x = Ks, col = cols_rs_many[as.character(r)])
}

legend('topright', col = cols_rs, legend = round(as.numeric(names(cols_rs)),2), 
       title = expression(italic(r)), bty = 'n', lwd = 2, cex = 0.7)


#============================================================
# Precision for equifrequent alleles against K for different r
# 3d plot (consider presenting in hap map)
#============================================================}
if(PLOT3d){
  plot_ly(z = ~FIM_surface, x = rs, y = Ks) %>% add_surface(colors = "PRGn") # 3d plot
  Surface0.5 = array(0.55, dim = dim(FIM_surface))
  plot_ly(z = ~sqrt(1/FIM_surface), x = ~rs, y = ~Ks) %>% add_surface(colors = "BuGn") %>% add_surface(z = ~Surface0.5, opacity = 0.8)
}


#==============================
# Gains plots many rs 
#==============================
FIM_2 = t(array(FIM_surface['2',], dim = rev(dim(FIM_surface))))
Gain = (sqrt(1/FIM_2)-sqrt(1/FIM_surface))/sqrt(1/FIM_2)*100
if(PLOT3d){
  plot_ly(z = ~Gain, x = ~rs, y = ~Ks) %>% add_surface(colors = "BuGn")
}

if(PLOThm){
  image(Gain, col = rev(sapply(rainbow(200, end = 0.6), adjustcolor, alpha.f = 0.5)), 
        xlab = expression(italic(K)), ylab = expression(italic(r)), xaxt = 'n', yaxt = 'n')
  contour(Gain, add = T)
  rs. = round(seq(min(rs), max(rs), by = 0.1),2)
  Ks. = round(seq(min(Ks), max(Ks), by = 1),2)
  axis(side = 1, at = (Ks. - min(Ks.))/(max(Ks.)-min(Ks.)), labels = Ks., cex.axis = 0.5)
  axis(side = 2, at = (rs. - min(rs.))/(max(rs.)-min(rs.)), labels = rs., las = 2)
}


#==============================
# Pct decrease in stdev vs Keff at single m 
#==============================
# Generate frequencies inc those where 
# Keff < K
p1 <- 1.75 # Multiplier given to large freq
All_fs <- lapply(2:max(Ks), function(K){
  feq = rep(1/K, K) # Equal frequencies
  f1 = p1/K # Large unequal frequency
  fr = rep((1-f1)/(K-1), K-1) # Remainding unequal frequencies
  cbind(feq = feq, fun = c(f1,fr))
})

# Generate hs
h1s = sapply(All_fs, function(X){
  colSums(X^2)
})

# Generate precisions
FIMs = lapply(rs, function(r){
  sapply(All_fs, function(X){
    apply(X, 2, FIM_fs, r = r)
  })
})
names(FIMs) = rs

# Calculate gain relatative to K = 2
GAIN = lapply(FIMs, function(X){
  X2 = X['feq',1] 
  (sqrt(1/X2)-sqrt(1/X))/sqrt(1/X2)*100
})

# Plots of pct decrease in standard deviation
XLIM = c(1, max(1/h1s))
YLIM = c(-50,80)
par(mfrow = c(1,2), mar = c(5.1,4.1,4.1,2.1)
)

# For K 
plot(NULL, ylim = YLIM, xlim = XLIM, las = 1, panel.first = grid(), 
     ylab = 'Percentage decrease in standard deviation', 
     xlab = expression(italic(K)[t]))
for(r in rs){
  Y = GAIN[[as.character(r)]]; COL = cols_rs[as.character(r)]
  points(y = Y['feq',], x = 1/h1s['feq',], pch = 16, col = COL, cex = 0.75)
  lines(y = Y['feq',], x = 1/h1s['feq',], col = COL, lwd = 0.75)
  
  # # Check agreement with FIM_K - yes
  # CHECK_Y = (sqrt(1/FIM_K(2,r)) - sqrt(1/sapply(Ks, FIM_K, r))) / sqrt(1/FIM_K(2,r)) * 100  
  # points(y = CHECK_Y, x = 1/h1s['feq',], pch = 1) # Check
}
legend('bottomright', inset = 0.1, col = cols_rs, title = expression(italic(r)), 
       legend = round(rs, 2), lwd = 2, bty = 'n')

# For K eff
plot(NULL, ylim = YLIM, xlim = XLIM, las = 1, panel.first = grid(), 
     ylab = 'Percentage decrease in standard deviation', 
     xlab = expression(italic(K)[t]~{"'"}))
for(r in rs){
  Y = GAIN[[as.character(r)]]; COL = cols_rs[as.character(r)]
  points(y = Y['feq',], x = 1/h1s['feq',], pch = 16, col = COL, cex = 0.75)
  points(y = Y['fun',], x = 1/h1s['fun',], pch = 17, col = COL, cex = 0.75)
  
  # Add lines 
  Y = as.vector(Y); X = rep(1/as.vector(h1s))
  X = sort.int(X, index.return = T)
  lines(y = Y[X$ix], x = X$x, col = COL, lwd = 0.75)
}
legend('bottomright', inset = 0.1, pch = c(16,17), 
       legend = c(expression(italic(K)[t]~{"'"}==italic(K)[t]), 
                  expression(italic(K)[t]~{"'"}<italic(K)[t])), cex = 1, lwd = 2, bty = 'n')


# # Average over rs 
# Y = rbind(rowMeans(sapply(FIMs, function(X)X['feq',])), 
#           rowMeans(sapply(FIMs, function(X)X['fun',])))
# points(y = Y[1,], x = 1/h1s['feq',], pch = 17, col = 'black')
# points(y = Y[2,], x = 1/h1s['fun',], pch = 20, col = 'black')
# # Add lines 
# Y = as.vector(Y); X = rep(1/as.vector(h1s))
# Y = sort.int(Y, index.return = T)
# lines(y = Y$x, x = X[Y$ix])


#==============================
# multiplicative increase in informativeness vs Keff at single m 
#==============================
# Calculate gain relatative to K = 2
GAIN = lapply(FIMs, function(X){
  X2 = X['feq',1] 
  #((X-X2)/X2) *100
  X/X2  
})

par(mfrow = c(1,2), mar = c(5,2,2,0.5), oma = c(0,2,0,0))
plot(NULL, ylim = c(0,7), xlim = c(1, max(1/h1s)), las = 1, panel.first = grid(), ylab = '', 
     xlab = expression('Effective marker cardinality,'~italic(K)[t]~"'"), 
     main = 'Equifrequent alleles', cex.lab = 1.2, cex.main = 1.5)

# Outer label
title(ylab = expression('Multiplicative increase in informativeness:'~FIM[italic(t)](italic(K)[t]~"'"~','~italic(r))%/%FIM[italic(t)](italic(K)[t]~"'" == 2~','~italic(r))), 
      cex.lab = 1.2, outer = T, line = 0, adj = 0.7)

for(r in rs){
  Y = GAIN[[as.character(r)]]; COL = cols_rs[as.character(r)]
  points(y = Y['feq',], x = 1/h1s['feq',], pch = 16, col = COL, cex = 0.75)
  lines(y = Y['feq',], x = 1/h1s['feq',], col = COL, lwd = 0.75)
}
legend('topleft', inset = 0.01, col = cols_rs, title = expression(italic(r)), 
       legend = round(rs, 2), lwd = 2, bty = 'n')
legend('bottomright', inset = 0.1, pch = 19, cex = 1.2, 
       legend = expression(italic(K)[t]~"'"== italic(K)[t]), lwd = 2, bty = 'n')

# For K eff
par(mar = c(5,0.5,2,2))
plot(NULL, ylim = c(0,7), xlim = c(1, max(1/h1s)), las = 1, panel.first = grid(), ylab = '', 
     xlab = expression('Effective marker cardinality,'~italic(K)[t]~"'"), yaxt = 'n', 
     main = 'Non-equifrequent alleles',cex.lab = 1.2, cex.main = 1.5)

for(r in rs){
  Y = GAIN[[as.character(r)]]; COL = cols_rs[as.character(r)]
  points(y = Y['feq',], x = 1/h1s['feq',], pch = 16, col = COL, cex = 0.75)
  points(y = Y['fun',], x = 1/h1s['fun',], pch = 17, col = COL, cex = 0.75)
  
  # Add lines 
  Y = as.vector(Y); X = rep(1/as.vector(h1s))
  X = sort.int(X, index.return = T)
  lines(y = Y[X$ix], x = X$x, col = COL, lwd = 0.75)
}
legend('bottomright', inset = 0.1, pch = c(19,17), cex = 1.2, 
       legend = c(expression(italic(K)[t]~"'"== italic(K)[t]),
                  expression(italic(K)[t]~"'"<italic(K)[t])), lwd = 2, bty = 'n')




#==============================
# Plot of absolute sd against Keff
#==============================
par(mfrow = c(3,3), mar = c(3,3,1,1))

for(r in rs){
  Y = sqrt(1/FIMs[[as.character(r)]])
  XLIM = range(1/h1s)
  COL = cols_rs[as.character(r)]
  plot(NULL, ylim = range(Y), xlim = XLIM, las = 1, panel.first = grid())
  points(y = Y['feq',], x = 1/h1s['feq',], pch = 17, col = COL)
  points(y = Y['fun',], x = 1/h1s['fun',], pch = 16, col = COL)
  
  # Add lines 
  Y = as.vector(Y); X = rep(1/as.vector(h1s))
  Y = sort.int(Y, index.return = T)
  lines(y = Y$x, x = X[Y$ix], col = COL)
}

# Average over rs
Y = rbind(rowMeans(sapply(FIMs, function(X)X['feq',])), 
          rowMeans(sapply(FIMs, function(X)X['fun',])))
plot(NULL, ylim = range(Y), xlim = XLIM, las = 1, panel.first = grid())
points(y = Y[1,], x = 1/h1s['feq',], pch = 20, col = 'black')
points(y = Y[2,], x = 1/h1s['fun',], pch = 17, col = 'black')
# Add lines 
Y = as.vector(Y); X = rep(1/as.vector(h1s))
Y = sort.int(Y, index.return = T)
lines(y = Y$x, x = X[Y$ix])



#===========================================
# Check that heuristic applies for numerically 
# drawn frequencies also: yes providing 
# r > 0.15 (ish)
#===========================================
load('../RData/tables_many_repeats_bigK.RData') # Load fixed m results and data sets 

# Get data set names and set points
names_datasets = names(tables_many_repeats_K_dataset)
pchs = rep(c(20,17),each = 10)
names(pchs) = names_datasets[-1]

# For the set of 96 loci, calculate h
h = sapply(names_datasets[-1],function(x){
  mean(rowSums(tables_many_repeats_K_dataset[[x]]^2))
})

Ylabs = rep(c('Average FIM (m = 96)', 'Average sd (m = 96)'), 2)
Xlabs = rep(c('h', 'K'), each = 2)

# Calculate avFIM for all rs
avFIMs = sapply(rs, function(r){
  avFIM = sapply(names_datasets[-1],function(x){
    mean(apply(tables_many_repeats_K_dataset[[x]], 1, FIM_fs, r))
  })})
colnames(avFIMs) = rs

par(mfrow = c(3,3), mar = c(3,3,1,1))
for(i in 1:length(Ylabs)){
  
  for(r in rs){
    
    # Sort if plotting precision of sd
    if(Ylabs[i] == 'Average FIM, m = 96'){
      Y = avFIMs[,as.character(r)]
    } else {
      Y = sqrt(1/avFIMs[,as.character(r)])}
    
    # Sort if plotting h or Keff
    if(Xlabs[i] == 'h'){
      X = h; Xlab = expression(italic(h)[m])
    } else {
      X = 1/h; expression('Average'~italic(K)~{"'"})}
    
    # Order ys 
    Y = sort.int(Y, index.return = T) 
    
    # Plot for specific r
    plot(y = Y$x, x = X[Y$ix], type = 'b', panel.first = grid(),
         pch = pchs[names(Y$x)], col = cols_rs[as.character(r)], 
         xlab = '', ylab = '', yaxt = 'n', xaxt = 'n')
    axis(side = 2, las = 2, line = -0.5, tick = F)
    axis(side = 1, line = -0.5, tick = F)
    title(ylab = Ylabs[i], xlab = Xlab, line = 1.85)
  }
  
  # Average over all rs
  if(Ylabs[i] == 'Average FIM, m = 96'){
    Y = avFIMs} else {Y = sqrt(1/avFIMs)}
  Y = rowMeans(Y)
  Y = sort.int(Y, index.return = T)
  
  plot(y = Y$x, x = X[Y$ix], type = 'b', pch = pchs[names(Y$x)], 
       panel.first = grid(), yaxt = 'n', xaxt = 'n')
  axis(side = 2, las = 2, line = -0.5, tick = F)
  axis(side = 1, line = -0.5, tick = F)
  legend('topright', col = c(cols_rs,"#000000FF"), title = expression(italic(r)), 
         legend = c(round(rs, 2),'Av.'), pch = pchs[names(Y$x)], bty = 'n')
  title(ylab = Ylabs[i], xlab = Xlab, line = 1.85)
}

if(JPEG){dev.off()}




