
library(ahaz)
library(glmnet)
library(doParallel)
library(doRNG)

source("core/ATE_ada.R")
source("MS CLIME analysis/ate_dich.R")

dat.name.list = paste("reg", "-data-",
                      c("RN","DF"), "_MStrt",sep = "")
fit.name.list = paste("reg", "-data-",
                      c("RN","DF"), "_MStrt",sep = "")

# dat.name = dat.name.list[3]
np=10
cl = makeCluster(getOption("cl.cores", np))
registerDoParallel(cl)
registerDoRNG(seed = 531)

B = 10000
iset=2

for(iset in 1:length(dat.name.list))
{
  dat.name = dat.name.list[iset]
  fit.name = fit.name.list[iset]
  load( paste("MS CLIME analysis/data-analysis/",dat.name,".rda",sep=''))
  Z = Z[,c("FEMALE", "RACE", "AGE_AT_FIRSTMSICD", "FOLLOWUP_DURA", "DISEASE_DURA",
           "PRIOR_RELAPSE_12MONS", "PRIOR_RELAPSE_24MONS")]
  
  res.file = paste("MS CLIME analysis/result/MS_result_new_",fit.name,".rda",sep='')
  
  fit.alltil17 = ATE.ada(surv,D,Z)
  
  boot <-foreach(b = 1:B, .packages=c("ahaz","glmnet","survival")) %dopar%
    {
      bootid = sample(1:length(D),replace = T)
      
      bootfit = try(ATE.ada(Surv(jitter(surv[bootid,1], 1e-7),surv[bootid,2]),
                            D[bootid],Z[bootid,],
                            penalties = fit.alltil17$penalties),TRUE)
      
      if(inherits(bootfit,"try-error"))
      {
        return(NULL)
      }else
      {
        return(bootfit)
      }
    }
  
  fit.1yr = ate.dich.ada(as.numeric(surv[,1]<=1),D, Z)
  boot.1yr <-foreach(b = 1:B, .packages=c("glmnet")) %dopar%
    {
      bootid = sample(1:length(D),replace = T)
      
      bootfit = try(ate.dich.ada(as.numeric(surv[bootid,1]<=1),D[bootid], Z[bootid,],
                                 or1.pen = fit.1yr$or1.pen, or0.pen = fit.1yr$or0.pen,
                                 or1.lambda = fit.1yr$or1.lambda, or0.lambda = fit.1yr$or0.lambda,
                                 ps.pen = fit.1yr$ps.pen, ps.lambda = fit.1yr$ps.lambda),TRUE)
      
      if(inherits(bootfit,"try-error"))
      {
        return(NULL)
      }else
      {
        return(bootfit)
      }
    }
  
  fit.2yr = ate.dich.ada(as.numeric(surv[,1]<=2),D, Z)
  boot.2yr <-foreach(b = 1:B, .packages=c("glmnet")) %dopar%
    {
      bootid = sample(1:length(D),replace = T)
      
      bootfit = try(ate.dich.ada(as.numeric(surv[bootid,1]<=2),D[bootid], Z[bootid,],
                                 or1.pen = fit.2yr$or1.pen, or0.pen = fit.2yr$or0.pen,
                                 or1.lambda = fit.2yr$or1.lambda, or0.lambda = fit.2yr$or0.lambda,
                                 ps.pen = fit.2yr$ps.pen, ps.lambda = fit.2yr$ps.lambda),TRUE)
      
      if(inherits(bootfit,"try-error"))
      {
        return(NULL)
      }else
      {
        return(bootfit)
      }
    }
  
  save(fit.alltil17, boot, 
       fit.1yr, boot.1yr, 
       fit.2yr, boot.2yr,
       file =res.file)
}

stopCluster(cl)
