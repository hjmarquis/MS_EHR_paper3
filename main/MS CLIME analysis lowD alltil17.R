
library(ahaz)
library(glmnet)
library(doParallel)

version = "v2"
load(paste("MS CLIME analysis/data-analysis/dataLD_allyear_",version,".rda",sep=''))

source("core/ATEdata.R")

fitLD.alltil17 = out
rm(out)

save(fitLD.alltil17, file = paste("MS CLIME analysis/result/MS_LDresult_alltil17_",version,".rda",sep=''))
# save(fitLD.alltil17, file = "data analysis/SM-hrr-inter-5f.rda")
# load("data analysis/SM-hrr-inter-5f.rda")

(1-pnorm(abs(fitLD.alltil17$ate$hat/fitLD.alltil17$ate$hat.sd)))*2
(1-pnorm(abs(fitLD.alltil17$atecf$hat/fitLD.alltil17$atecf$hat.sd)))*2
(1-pnorm(abs(fitLD.alltil17$ate$checkhat/fitLD.alltil17$ate$checkhat.sd)))*2
(1-pnorm(abs(fitLD.alltil17$atecf$checkhat/fitLD.alltil17$atecf$checkhat.sd)))*2

hthetabeta =  fitLD.alltil17$fit$hthetabeta
hgr = fitLD.alltil17$fit$hgr

fitLD.alltil17$conf.coef = colnames(Z)[hthetabeta[-1]!=0 & (hgr[-1]!=0)]
fitLD.alltil17$Znames = colnames(Z)

save(fitLD.alltil17, file = paste("MS CLIME analysis/result/MS_LDresult_alltil17_",version,".rda",sep=''))
