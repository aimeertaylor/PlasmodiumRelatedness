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
library(kableExtra) # For latex tables
par(family = 'serif')
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
if(PLOT_pdf){pdf(file = '../Plots/Plot_tables.pdf', height = 7, width = 7)}


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

kable(tables_snp_requirements$`Proportional to MAF under hmm`, format = 'latex')
kable(tables_snp_requirements$`Uniformly at random under hmm`, format = 'latex')
#save(tables_snp_requirements)


#====================================================
# Plot average MSEs against m for different r
#====================================================
par(par_df)

# Generate a plot per ft strategy
for(ft_strategy in c('Proportional to MAF','Uniformly at random')){
  
  Table <- Tables[['hmm']][,,ft_strategy] # Extract Table
  ms <- as.numeric(rownames(Table))
  rs <- as.numeric(colnames(Table))
  cols <- rainbow(length(rs), end = 0.8)
  max_m <- max(ms)
  plot(NULL, 
       ylim = c(0, pmax(max(Table), 0.35)), 
       xlim = range(ms), bty = 'n', xaxt = 'n', las = 2, 
       ylab = 'Root mean squared error', xlab = expression('Number of markers,'~italic(m)), 
       main = paste(ft_strategy, 'under', names(Tables['hmm'])))
  axis(side = 1, at = ms)
  for(j in 1:ncol(Table)){
    lines(y = Table[,j], x = ms, pch = 16, panel.first = grid(), type = 'b', 
          col= cols[j])
  }
  legend('topright', bty = 'n', pch = 16, col = cols, 
         legend = format(rs, digits = 2, drop0trailing = F), 
         title = expression('Relatedness,'~italic(r)))
}




#====================================================
# Comparing hmm to iid ('Proportional to MAF')
#====================================================
rfixed = 0.5
par(par_df)
ft_strategy <- 'Proportional to MAF'

# With variable r
hmm <- tables_many_repeats_m_hmm[,,ft_strategy]
iid <- tables_many_repeats_m_iid[,,ft_strategy]
ms <- as.numeric(rownames(iid))
rs <- colnames(hmm)
cols <- rainbow(length(rs), end = 0.8); 
names(cols) <- rs
plot(NULL, 
     ylim = c(0, pmax(max(max(hmm, iid),0.25))), 
     xlim = range(ms), bty = 'n', xaxt = 'n', las = 2, 
     ylab = 'Root mean squared error', xlab = expression('Number of markers,'~italic(m)), 
     main = paste(ft_strategy),panel.first = grid())
axis(side = 1, at = ms)
for(r in rs){
  lines(y = hmm[,r], x = ms, pch = 16, type = 'b', col= cols[r])
  lines(y = iid[,r], x = ms, pch = 21, type = 'b', col = cols[r], bg = 'white', lty = 'dashed')
}

legend(x = 96, y = 0.25, bty = 'n', pch = 16, col = cols[rs], legend = rs, lwd = 1, 
       title = expression('Relatedness,'~italic(r)))
legend(x = 288, y = 0.25, bty = 'n', lty = 1:2, pch = c(16,21), col = 'black', pt.bg = 'white', 
       legend = c('HMM', 'Independence model'), inset = 0.2)


# With variable k
hmm <- tables_many_repeats_k_hmm[,,ft_strategy]
iid <- tables_many_repeats_k_iid[,,ft_strategy]
ms <- as.numeric(rownames(iid))
ks <- colnames(hmm)
cols <- rainbow(length(ks), end = 0.8); 
names(cols) <- ks
plot(NULL, 
     ylim = c(0, pmax(max(max(hmm, iid),0.25))), 
     xlim = range(ms), bty = 'n', xaxt = 'n', las = 2, 
     ylab = 'Root mean squared error', xlab = expression('Number of markers,'~italic(m)), 
     main = paste(ft_strategy),panel.first = grid())
axis(side = 1, at = ms)
for(k in ks){
  lines(y = hmm[,k], x = ms, pch = 16, type = 'b', col= cols[k])
  lines(y = iid[,k], x = ms, pch = 21, type = 'b', col = cols[k], bg = 'white', lty = 'dashed')
}

legend(x = 96, y = 0.25, bty = 'n', pch = 16, col = cols[ks], legend = ks, lwd = 1, 
       title = expression('Switch rate parameter,'~italic(k)))
legend(x = 288, y = 0.25, bty = 'n', lty = 1:2, pch = c(16,21), col = 'black', pt.bg = 'white', 
       legend = c('HMM', 'Independence model'), inset = 0.2)



if(PLOT_pdf){dev.off()}

