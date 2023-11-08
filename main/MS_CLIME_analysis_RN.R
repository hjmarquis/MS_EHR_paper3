
.libPaths(c("~/R-3.6.1/library",.libPaths()))
library(ahaz)
library(glmnet)
library(doParallel)
library(doRNG)

source("core/ATE_ada.R")

dat.name.list = paste(rep(c("low-d","high-d","comb"),4), "-data-",
                      rep(c("RN","DF"),each = 6), '_',
                      rep(c("MStrt","RxNorm"),2,each=3),sep = "")
cl = makeCluster(getOption("cl.cores", 20))
registerDoParallel(cl)
registerDoRNG(seed = 531)
clusterEvalQ(cl, .libPaths(c("~/R-3.6.1/library",.libPaths())))

B = 10000

for(dat.name in dat.name.list[c(1,3)])
{
  load( paste("MS CLIME analysis/data-analysis/",dat.name,"p10.rda",sep=''))
  
  res.file = paste("MS CLIME analysis/result/MS_result_cvboot_",dat.name,".rda",sep='')
  # if(!file.exists(res.file))
  # {
    fit.alltil17 = ATE.ada(surv,D,Z)
    hthetabeta =  fit.alltil17$fit$hthetabeta
    hgr = fit.alltil17$fit$hgr
  
    fit.alltil17$conf.coef = colnames(Z)[hthetabeta[-1]!=0 & (hgr[-1]!=0)]
    fit.alltil17$Znames = colnames(Z)
    
    save(fit.alltil17, file =res.file)
  # }
  
  # load.list = load(res.file)
  # 
  # if(is.element("boot",load.list))
  # {
  #   next
  # }

    
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
  save(fit.alltil17, boot, boot.cv, boot.cv2, file =res.file)
}

stopCluster(cl)