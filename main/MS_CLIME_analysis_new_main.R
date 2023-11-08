
library(ahaz)
library(glmnet)

source("core/ATE_ada.R")

dat.name.list = paste(rep(c("low-d","high-d","comb"),4), "-data-",
                      rep(c("RN","DF"),each = 6), '_',
                      rep(c("MStrt","RxNorm"),2,each=3),sep = "")

# dat.name = dat.name.list[3]
for(dat.name in dat.name.list)
{
  load( paste("MS CLIME analysis/data-analysis/",dat.name,".rda",sep=''))
  
  res.file = paste("MS CLIME analysis/result/MS_result_new_",dat.name,".rda",sep='')
  
  # while(fit.alltil17$or<fit.alltil17$dr$est)
  # {
  # if(!file.exists(res.file))
    fit.alltil17 = ATE.ada(surv,D,Z)
    hthetabeta =  fit.alltil17$fit$hthetabeta
    hgr = fit.alltil17$fit$hgr
  
    fit.alltil17$conf.coef = colnames(Z)[hthetabeta[-1]!=0 & (hgr[-1]!=0)]
    fit.alltil17$Znames = colnames(Z)
  # }
    
    save(fit.alltil17, file =res.file)
    
    print(dat.name)
    print(exp(-c(fit.alltil17$or,fit.alltil17$ipw,
                 fit.alltil17$dr$est)*2))
    print(exp(-c(fit.alltil17$hdi0,fit.alltil17$hdi)*2))
    print(pnorm(-abs(fit.alltil17$dr$est/fit.alltil17$dr$sd))*2)
    print(c(sum(fit.alltil17$fit$hthetabeta!=0),
            sum(fit.alltil17$fit$hgr!=0)))
    # print(fit.alltil17$penalties)
  
}
