
library(ahaz)
library(glmnet)
library(doParallel)
library(doRNG)

source("core/ATE_ada.R")

dat.name.list = paste(rep(c("low-d","high-d","comb"),4), "-data-",
                      rep(c("RN","DF"),each = 6), '_',
                      rep(c("MStrt","RxNorm"),2,each=3),sep = "")

B = 100

dat.name = dat.name.list[3]

# 
# for(dat.name in dat.name.list)

  load( paste("MS CLIME analysis/data-analysis/",dat.name,"p10.rda",sep=''))
  
  # res.file = paste("MS CLIME analysis/result/MS_result_cvboot_",dat.name,".rda",sep='')
fit.alltil17 = ATE.ada(surv,D,Z)
    hthetabeta =  fit.alltil17$fit$hthetabeta
    hgr = fit.alltil17$fit$hgr
  
    fit.alltil17$conf.coef = colnames(Z)[hthetabeta[-1]!=0 & (hgr[-1]!=0)]
    fit.alltil17$Znames = colnames(Z)
    

  

  boot.cv <-foreach(b = 1:B, .packages=c("ahaz","glmnet","survival"),
                 .combine = rbind) %dopar%
  {
    bootid = sample(1:length(D),replace = T)
    
    bootfit = try(ATE.ada(Surv(jitter(surv[bootid,1], 1e-8),surv[bootid,2]),
                          D[bootid],Z[bootid,]),TRUE)
    
    if(inherits(bootfit,"try-error"))
    {
      return(NULL)
    }else
    {
      return(c(bootfit$or, bootfit$ipw,
               bootfit$hdi0, bootfit$dr$est,  bootfit$hdi,
               bootfit$dr.cf, bootfit$hdi.cf))
    }
  }
  
  ridge.pen = fit.alltil17$penalties
  ridge.pen$lambda.beta.nopen = NULL
  ridge.pen$lambda.hgr = NULL
  boot.cv2 <-foreach(b = 1:B, .packages=c("ahaz","glmnet","survival"),
                    .combine = rbind) %dopar%
    {
      bootid = sample(1:length(D),replace = T)
      
      bootfit = try(ATE.ada(Surv(jitter(surv[bootid,1], 1e-8),surv[bootid,2]),
                            D[bootid],Z[bootid,],
                            penalties =  ridge.pen),TRUE)
      
      if(inherits(bootfit,"try-error"))
      {
        return(NULL)
      }else
      {
        return(c(bootfit$or, bootfit$ipw,
                 bootfit$hdi0, bootfit$dr$est,  bootfit$hdi,
                 bootfit$dr.cf, bootfit$hdi.cf))
      }
    }
  
  boot <-foreach(b = 1:B, .packages=c("ahaz","glmnet","survival"),
                    .combine = rbind) %dopar%
    {
      bootid = sample(1:length(D),replace = T)
      
      bootfit = try(ATE.ada(Surv(jitter(surv[bootid,1], 1e-8),surv[bootid,2]),
                            D[bootid],Z[bootid,],
                            penalties =  fit.alltil17$penalties),TRUE)
      
      if(inherits(bootfit,"try-error"))
      {
        return(NULL)
      }else
      {
        return(c(bootfit$or, bootfit$ipw,
                  bootfit$hdi0,  bootfit$dr$est,bootfit$hdi,
                 bootfit$dr.cf, bootfit$hdi.cf))
      }
    }

  full.est = c(fit.alltil17$or, fit.alltil17$ipw,
               fit.alltil17$hdi0,   fit.alltil17$dr$est, fit.alltil17$hdi,
               fit.alltil17$dr.cf, fit.alltil17$hdi.cf)
  exp(full.est*2)
  
  format(apply(sign(boot)== -sign(matrix(rep(full.est,
                         each = nrow(boot)),nrow(boot))),2,mean)*2,
                                        digits=0, nsmall = 3)
  format(apply(sign(boot.cv)== -sign(matrix(rep(full.est,
         each = nrow(boot.cv)),nrow(boot.cv))),2,mean)*2,
         digits=0, nsmall = 3)
  
  format(apply(sign(boot.cv2)== -sign(matrix(rep(full.est,
                                                each = nrow(boot.cv2)),nrow(boot.cv2))),2,mean)*2,
         digits=0, nsmall = 3)
  
