
library(ahaz)
library(glmnet)
library(doParallel)

version = "v1"
load(paste("MS CLIME analysis/data-analysis/dataLD_18m_",version,".rda",sep=''))

source("core/ATEdataLD.R")

fit.18m = out
rm(out)

save(fit.18m, file = paste("MS CLIME analysis/result/MS_LDresult_18m_",version,".rda",sep=''))
# save(fit.18m, file = "data analysis/SM-hrr-inter-5f.rda")
# load("data analysis/SM-hrr-inter-5f.rda")

(1-pnorm(abs(fit.18m$ate$hat/fit.18m$ate$hat.sd)))*2
(1-pnorm(abs(fit.18m$atecf$hat/fit.18m$atecf$hat.sd)))*2
(1-pnorm(abs(fit.18m$ate$checkhat/fit.18m$ate$checkhat.sd)))*2
(1-pnorm(abs(fit.18m$atecf$checkhat/fit.18m$atecf$checkhat.sd)))*2

hthetabeta =  fit.18m$fit$hthetabeta
hgr = fit.18m$fit$hgr

fit.18m$conf.coef = colnames(Z)[hthetabeta[-1]!=0 & (hgr[-1]!=0)]
fit.18m$Znames = colnames(Z)

save(fit.18m, file = paste("MS CLIME analysis/result/MS_LDresult_18m_",version,".rda",sep=''))
