
library(ahaz)
library(glmnet)
library(doParallel)

version = "v2"
load(paste("MS CLIME analysis/data-analysis/dataHD_allyear_",version,".rda",sep=''))

source("core/ATEdata.R")

fit.alltil17 = out
rm(out)

save(fit.alltil17, file = paste("MS CLIME analysis/result/MS_HDresult_alltil17_",version,".rda",sep=''))
# save(fit.alltil17, file = "data analysis/SM-hrr-inter-5f.rda")
# load("data analysis/SM-hrr-inter-5f.rda")

(1-pnorm(abs(fit.alltil17$ate$hat/fit.alltil17$ate$hat.sd)))*2
(1-pnorm(abs(fit.alltil17$atecf$hat/fit.alltil17$atecf$hat.sd)))*2
(1-pnorm(abs(fit.alltil17$ate$checkhat/fit.alltil17$ate$checkhat.sd)))*2
(1-pnorm(abs(fit.alltil17$atecf$checkhat/fit.alltil17$atecf$checkhat.sd)))*2

hthetabeta =  fit.alltil17$fit$hthetabeta
hgr = fit.alltil17$fit$hgr

fit.alltil17$conf.coef = colnames(Z)[hthetabeta[-1]!=0 & (hgr[-1]!=0)]
fit.alltil17$Znames = colnames(Z)

save(fit.alltil17, file = paste("MS CLIME analysis/result/MS_HDresult_alltil17_",version,".rda",sep=''))
