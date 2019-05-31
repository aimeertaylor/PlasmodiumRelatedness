#====================================
# Extracting a typical k for r approx. 0.5 
# 8 based on Thai data set
#====================================
rm(list = ls())
load('../RData/TM_WGS_mles.RData')
rhats = mle_df$rhat_hmm
ind <- which(rhats < 0.525 & rhats > 0.475)
mean(mle_df$khat_hmm[ind]) # 8
range(mle_df$khat_hmm[ind])
plot(mle_df$khat_hmm[ind], mle_df$rhat_hmm[ind])
plot(y = mle_df$khat_hmm, x = mle_df$rhat_hmm, ylim = c(0,50))
