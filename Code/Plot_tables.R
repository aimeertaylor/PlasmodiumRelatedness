### This script plots tables using both the i.i.d model and the HMM 
### (including genotyping errors and allowing multiallelic loci)
### 
### The tables show the Root mean squared errors (MSEs)
### for various sample size, various rs, various fs, various Ks
###
### Terminology: f refers to a frequency of an allele at each position in the genome, 
### where there are K possible alleles at a given position 
### r refers to relatedness

#====================================================
## Set up
#====================================================
rm(list = ls())
dev.off()
PLOT_pdf <- T
library(RColorBrewer) # For color mapper functions
library(fields) # For image plot
par_df <- par(no.readonly = T)
green_mapper <- function(num){
  y <- colorRampPalette(brewer.pal(9, "Greens"))(num)
  z <- sapply(y, adjustcolor, alpha.f = 0.5)
  return(z)
}

# Load results
load(file = '../RData/tables_many_repeats_m.RData') 
load(file = '../RData/tables_asymptotics_m.RData')
load(file = '../RData/tables_many_repeats_littlek.RData') 
load(file = '../RData/tables_many_repeats_bigk.RData') 


# Set file for plots
if(PLOT_pdf){
  pdf(file = '../Plots/Plot_tables.pdf', 
      height = 8, width = 8)}


# ====================================================
# Visual check that rhat are normally distributed
# with mean 0 and variance 1/m * sigma_m^2
# ====================================================
names_thetas <- names(tables_many_repeats_m_thetas)
thetasmaf <- names_thetas[grepl('Proportional to MAF', names_thetas)] # Extract rs for prop to maf
thetas05maf <- lapply(thetasmaf, function(x){tables_many_repeats_m_thetas[[x]]$`0.5`})
rs05maf <- lapply(thetas05maf, function(x){x[,c(1,3)]})
par(mfrow = c(1,3), family = 'serif', pty = 's')

for(i in c(1,3,5)){
  X = hist(rs05maf[[i]][,2], xlim = c(0,1), freq = FALSE, col = 'gray', 
           ylab = 'Density', 
           xlab = expression(italic(widehat(r)[m])),
           main = bquote(italic(m) == .(gsub('m', '', strsplit(thetasmaf[i], split = '_')[[1]][1]))))
  
  Var = var(rs05maf[[i]][,2])
  curve(dnorm(x, 0.5, sd = sqrt(Var)), add = TRUE, col = 'blue', lwd = 2)
  legend('topleft', legend = c('Empirical', 'Normal'), fill = c('gray', 'blue'), bty = 'n')
}


#====================================================
#====================================================
Tables = list(hmm = tables_many_repeats_m_hmm,  
              iid = tables_many_repeats_m_iid)


# Plot and check distances are the same for the different strategies
par(mfcol = c(6,2), mar = c(1,1,1,1))
for(ft_strategy in c('Proportional to MAF','Uniformly at random')){
    
    # Get datasets for distances and allele frequencies 
    ind_data <- grepl(ft_strategy, names(tables_many_repeats_m_dataset))
    
    # Plot distances between markers
    for(i in rev((1:length(tables_many_repeats_m_dataset))[ind_data])){
      data =  tables_many_repeats_m_dataset[[i]]
      dt_finite <- data$dt[is.finite(data$dt)]
      hist(dt_finite, col = 'gray', main = '')
    }}


# Plot frequencies and generate summaries hm and average k eff.
for(ft_strategy in c('Proportional to MAF','Uniformly at random')){
  
  # Get datasets for distances and allele frequencies 
  ind_data <- grepl(ft_strategy, names(tables_many_repeats_m_dataset))
  
  # Plot distances between markers
  for(i in rev((1:length(tables_many_repeats_m_dataset))[ind_data])){
    data =  tables_many_repeats_m_dataset[[i]]
    hist(data$fs, col = 'gray', main = '')
  }}

# Generate summaries h_m and Keffs
ind_MAF = grepl('Proportional to MAF', names(tables_many_repeats_m_dataset))
ind_unif = grepl('Uniformly at random', names(tables_many_repeats_m_dataset))

h_ms = sapply(tables_many_repeats_m_dataset, function(x){mean(x$fs^2 + (1-x$fs)^2)})
K_effs = sapply(tables_many_repeats_m_dataset, function(x){mean(1/(x$fs^2 + (1-x$fs)^2))})
h_ms_MAF = h_ms[ind_MAF]
h_ms_unif = h_ms[ind_unif]
K_effs_MAF = K_effs[ind_MAF]
K_effs_unif = K_effs[ind_unif]

round(c(mean(h_ms_MAF), mean(K_effs_MAF)),2)
round(c(mean(h_ms_unif), mean(K_effs_unif)), 2)

#====================================================
# Plot MSE for different m and r
#====================================================
# Generate a plot per ft strategy
for(ft_strategy in c('Proportional to MAF','Uniformly at random')){
  
  # Generate plots for both hmm and idd
  for(Table_ind in 1:length(Tables)){
    
    # Get datasets for distances and allele frequencies 
    ind_data <- grepl(ft_strategy, names(tables_many_repeats_m_dataset))
    
    # Extract table with ft strategy
    Table <- Tables[[Table_ind]][,,ft_strategy]
    Dims <- dim(Table)
    Dimnames <- dimnames(Table)
    
    # Set plotting parameters
    par(mar = c(2,1,1,1), plt = c(0,1,0,1), pty = 's', family = 'serif', omi = c(1,0.2,0.2,0.2))
    graphics::layout(mat = cbind(1:6,7:12, array(13,dim = Dims)))
    cols <- green_mapper(100)
    
    # Plot distances between markers
    for(i in rev((1:length(tables_many_repeats_m_dataset))[ind_data])){
      data =  tables_many_repeats_m_dataset[[i]]
      dt_finite <- data$dt[is.finite(data$dt)]
      # Note that setting xlim to max global distance makes plot uniformative
      X <- hist(dt_finite, main = '', col = 'gray', breaks = 10, border = 'white',
                ylab = '', xlab = '', xaxt = 'n', yaxt = 'n') 
      axis(side = 1, at = c(0,max(dt_finite)), 
           labels = c(0, format(max(dt_finite), scientific = T, digits = 1)),
           tick = F, line = -1)
      axis(side = 2, at = c(0,max(X$counts)), tick = F, line = -0.8, las = 2)
      abline(h = max(X$counts))
      abline(v = mean(dt_finite), col = 'blue', lty = 'dotted')
      m_num <- strsplit(names(tables_many_repeats_m_dataset)[i], split= '_')[[1]][1]
      title(ylab = 'Distance count', xlab = 'Distance', line = 0, cex.lab = 0.75)
      title(main = bquote(italic(m)==.(strsplit(m_num, split = 'm')[[1]][2])), line = 0.5)
    }
    
    # Plot allele frequencies
    for(i in rev((1:length(tables_many_repeats_m_dataset))[ind_data])){
      data =  tables_many_repeats_m_dataset[[i]]
      X <- hist(pmin(data$fs, 1-data$fs), main = '', col = 'gray', breaks = 10, border = 'white',
                ylab = '', xlab = '', xaxt = 'n', yaxt = 'n')
      axis(side = 1, at = c(0,0.5), labels = c('0','0.5'),tick = F, line = -1)
      axis(side = 2, at = max(X$counts), tick = F, line = -0.8, las = 2)
      abline(h = max(X$counts))
      abline(v = mean(pmin(data$fs, 1-data$fs)), col = 'blue', lty = 'dotted')
      m_num <- strsplit(names(tables_many_repeats_m_dataset)[i], split= '_')[[1]][1]
      title(ylab = 'SNP count', xlab = 'MAF', line = 0, cex.lab = 0.75)
      title(main = bquote(italic(m)==.(strsplit(m_num, split = 'm')[[1]][2])), line = 0.5)
    }
    
    
    # Plot table
    par(plt = c(0,0.9,0,1))
    image(t(Table), col = cols,yaxt = 'n', xaxt = 'n', ylab = '', xlab = '')
    
    # Annotate horizontal axis
    axis(side = 1, at = seq(0,1,length.out = Dims[2]), labels = Dimnames[[2]], 
         tick = FALSE, cex.axis  = 1.5)
    mtext(side = 1, text = expression('Relatedness'~italic(r)), line = 3, cex = 1)
    mtext(side = 1, text =  paste(ft_strategy, 'under', names(Tables[Table_ind])), 
          line = 4, cex = 1)
    
    # Annotate cells
    for(i in 1:Dims[1]){
      for(j in 1:Dims[2]){
        text(x = (j-1)/(Dims[2]-1), y = (i-1)/(Dims[1]-1), 
             labels = round(Table,2)[i, j], cex = 1.5)    
      }}
    
    # Add image plot legend separately 
    # (otherwise mtext results in text for layout.show(n = 1), 
    #  which after much exploration of par(c('usr', 'mfg', 'plt')) 
    #  and graphics.reset = F I do not understand). 
    image.plot(t(Table), col = cols, legend.only =  T,  
               smallplot = c(x1 = 0.93, x2 = 0.96, y1 = 0, y2 = 1),
               legend.args = list(text = ''))
    
  }
}


#====================================================
# Plot average MSEs against m for different r
#====================================================
par(par_df); par(mfrow = c(2,2), family = 'serif')
for(Table_ind in 1:length(Tables)){
  
  # Generate a plot per ft strategy
  for(ft_strategy in c('Proportional to MAF','Uniformly at random')){
    
    Table <- Tables[[Table_ind]][,,ft_strategy] # Extract Table
    ms <- as.numeric(rownames(Table))
    rs <- as.numeric(colnames(Table))
    cols <- brewer.pal(length(rs), 'Spectral')
    max_m <- max(ms)
    plot(NULL, 
         ylim = c(0, pmax(max(Table), 0.35)), 
         xlim = range(ms), bty = 'n', xaxt = 'n', las = 2, 
         ylab = 'Root mean squared error', xlab = expression('Number of markers,'~italic(m)), 
         main = paste(ft_strategy, 'under', names(Tables[Table_ind])))
    axis(side = 1, at = ms)
    for(j in 1:ncol(Table)){
      lines(y = Table[,j], x = ms, pch = 20, panel.first = grid(), type = 'b', 
            col= cols[j])
    }
    legend('topright', bty = 'n', pch = 20, col = cols, 
           legend = format(rs, digits = 2, drop0trailing = F), 
           title = expression('Relatedness,'~italic(r)))
  }
}


### Tabulate SNP requirements for different CI widths and r 
## Problematic when non-monotonic for a given r (e.g. r = 0.01 under iid uniform)
## Problem may go away when run with more simulations??
## Need to sort so it takes maximum
Av_MSE = seq(0.05,0.25,0.05)
tables_snp_requirements <- list()

for(Table_ind in 1:length(Tables)){
  # Generate a plot per ft strategy
  for(ft_strategy in c('Proportional to MAF','Uniformly at random')){
    Table <- Tables[[Table_ind]][,,ft_strategy] # Extract Table
    # Tabulate intervals
    table_snp_requirements = t(sapply(Av_MSE, function(CI_width){
      apply(Table, 2, function(x){
        # Add >480 s.t. as many CIs as intervals and sort
        z = sort.int(c(x,'>480' = 0, '<24' = 1),index.return = T) 
        interval = findInterval(CI_width,z$x) # Find interval
        # Create vectors (inside loop since CI widths do not always decreasing with m)
        intervals = sapply(2:length(z$x),function(i){paste(c(names(z$x)[i-1],":", names(z$x)[i]), collapse = '')})
        intervals[interval] # Return interval
      })
    }))
    # Replace any with ">480:" and ":<24" with just ">480" and "<24"
    table_snp_requirements[grepl(">480:",table_snp_requirements)] <- ">480"
    table_snp_requirements[grepl(":<24",table_snp_requirements)] <- "<24"
    rownames(table_snp_requirements) = Av_MSE
    tables_snp_requirements[[paste(ft_strategy, 'under', names(Tables[Table_ind]))]] <- table_snp_requirements
  }
}

print(tables_snp_requirements$`Proportional to MAF under hmm`, quote = F)
print(tables_snp_requirements$`Uniformly at random under hmm`, quote = F)
#save(tables_snp_requirements)

#====================================================
# Comparing hmm to iid ('Proportional to MAF')
#====================================================
rfixed = 0.5
par(par_df)
ft_strategy <- 'Proportional to MAF'
hmm <- tables_many_repeats_k_hmm[,,ft_strategy]
iid <- tables_many_repeats_k_iid[,,ft_strategy]
ms <- as.numeric(rownames(iid))
ks <- colnames(hmm)
cols <- brewer.pal(length(ks), 'Spectral'); 
names(cols) <- ks
plot(NULL, 
     ylim = c(0, pmax(max(max(hmm, iid),0.35))), 
     xlim = range(ms), bty = 'n', xaxt = 'n', las = 2, 
     ylab = 'Root mean squared error', xlab = expression('Number of markers,'~italic(m)), 
     main = paste(ft_strategy),panel.first = grid())
axis(side = 1, at = ms)
for(k in ks){
  lines(y = hmm[,k], x = ms, pch = 20, type = 'b', col= cols[k])
  lines(y = iid[,k], x = ms, pch = 1, type = 'b', col = cols[k], lty = 'dashed')
}

legend('top', bty = 'n', pch = 20, col = cols[ks], legend = ks, lwd = 1, 
       title = expression('Number of generations,'~italic(k)))
legend('top', bty = 'n', lty = 1:2, pch = c(20,1), col = 'black', 
       legend = c('HMM', 'i.i.d model'), inset = 0.2)


# ### Asymptotics
# ms = as.numeric(rownames(table_asymptotics_m_iid))
# rs = as.numeric(colnames(table_asymptotics_m_iid))
# 
# # check visually for iid case
# par(par_df); par(mfrow = c(3,2), mar = c(4,4,1,1), family = 'serif')
# for (icol in c(1:5)){
#   plot(x = ms, y = table_asymptotics_m_iid[,icol], type = "b", pch = 16, bty = 'n',
#        ylim = c(0, 1), xaxt = 'n', panel.first = grid(),
#        ylab = 'Root mean squared error', xlab = expression('Number of markers,'~italic(m)))
#   title(main = bquote(italic(r)==.(rs[icol])), line = -1)
#   lines(x = ms, y = tables_many_repeats_m_iid[,icol,1], lty = 2, type = 'b')
#   axis(side = 1, at = ms)
#   legend('topright', inset = 0.1, legend = c('Asymptotic result', 'Result based on many repeats'), 
#          pch = c(16,21), bty = 'n', lty =1:2, pt.bg = 'white')
# }
# 
# # check visually for hmm
# par(par_df); par(mfrow = c(3,2), mar = c(4,4,1,1), family = 'serif')
# for (icol in c(1:5)){
#   plot(x = ms, y = table_asymptotics_m_hmm[,icol], type = "b", pch = 16, bty = 'n',
#        ylim = c(0, 1), xaxt = 'n', panel.first = grid(),
#        ylab = 'Root mean squared error', xlab = expression('Number of markers,'~italic(m)))
#   title(main = bquote(italic(r)==.(rs[icol])), line = -1)
#   lines(x = ms, y = tables_many_repeats_m_hmm[,icol,1], lty = 2, type = 'b')
#   axis(side = 1, at = ms)
#   legend('topright', inset = 0.1, legend = c('Asymptotic result', 'Result based on many repeats'), 
#          pch = c(16,21), bty = 'n', lty =1:2, pt.bg = 'white')
# }



if(PLOT_pdf){dev.off()}

