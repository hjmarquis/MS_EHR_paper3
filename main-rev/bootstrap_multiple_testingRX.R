
library(ahaz)
library(glmnet)

source("core/ATE_ada.R")
source("MS CLIME analysis/ate_dich.R")


dat.name.list = paste0(rep(c("low-d-data-","comb-data-"),each = 2),
                       rep(c("RN","DF"),2),"_RxNorm")

B = 1000


dat.name = dat.name.list[1]

for(dat.name in dat.name.list)
{
  load( paste("MS CLIME analysis/data-analysis/",dat.name,"p10.rda",sep=''))
  Y1 = (surv[,1]<= 1)*surv[,2]
  Y2 = (surv[,1]<= 2)*surv[,2]
  
  fit.alltil17 = ATE.ada(surv,D,Z)
  fit1 = ate.dich.ada(Y1,D,Z)
  fit2 = ate.dich.ada(Y2,D,Z)
  
  boot.est = matrix(NA, 3, B)
  
  for (b in 1:B)
  {
    print(b)
    start.time = Sys.time()
    cont = TRUE
    
    while(cont)
    {
      bootid = sample(1:length(D),replace = T)
      
      boot.surv = Surv(jitter(surv[bootid,1], 1e-8),surv[bootid,2])
      bootfit.ahaz = try(ATE.ada(boot.surv,
                                 D[bootid],Z[bootid,],
                                 penalties =  fit.alltil17$penalties),TRUE)
      if(inherits(bootfit.ahaz,"try-error"))
      {
        next
      }
      
      bootfit.yr1 = try(ate.dich.ada(Y1[bootid],
                                     D[bootid],Z[bootid,],
                                     fit1$or1.pen, fit1$or0.pen, fit1$ps.pen,
                                     fit1$or1.lambda, fit1$or0.lambda, fit1$ps.lambda),TRUE)
      if(inherits(bootfit.yr1,"try-error"))
      {
        next
      }
      
      bootfit.yr2 = try(ate.dich.ada(Y2[bootid],
                                     D[bootid],Z[bootid,],
                                     fit2$or1.pen, fit2$or0.pen, fit2$ps.pen,
                                     fit2$or1.lambda, fit2$or0.lambda, fit2$ps.lambda),TRUE)
      if(inherits(bootfit.yr2,"try-error"))
      {
        next
      }
      cont = FALSE
    }
    boot.est[1,b] = bootfit.yr1$dr
    boot.est[2,b] = bootfit.yr2$dr
    boot.est[3,b] = bootfit.ahaz$dr$est
    run.time = Sys.time() - start.time
    print(run.time)
  }
  save(
    fit1, fit2, fit.alltil17,
    boot.est, file = paste0("MS CLIME analysis/result/MS_result_new_",dat.name,"p10_boot_joint.rda"))

}

i = 1
i = 4

load(paste0("MS CLIME analysis/result/MS_result_new_",dat.name.list[i],"p10_boot_joint.rda"))

c(-fit1$dr, -fit2$dr ,exp( fit.alltil17$dr$est*2))
c(fit1$dr, fit2$dr ,exp( -fit.alltil17$dr$est*2))

ci = apply(boot.est,1,quantile, probs = c(0.025, 0.975))
ci = -apply(boot.est,1,quantile, probs = c(0.025, 0.975))
ci[,3] = exp(-ci[,3]*2)

est = c(fit1$dr, fit2$dr ,fit.alltil17$dr$est)
apply(sign(boot.est)!= sign(est),1, mean)*2

sig =apply(apply(boot.est,1,quantile, probs = c(0.025, 0.975)),
           2,diff)/(1.96*2)
est = apply(boot.est, 1, median)
Z.wald = (est/sig)
boot.std = (boot.est - est) / sig

apply(outer(apply(abs(boot.std), 2, max), abs(Z.wald), '>'),2,mean)
# apply(outer(apply(boot.std, 2, min), Z.wald, '<'),2,mean)*2

# 
# yo = mvrnorm(10000, rep(0,3), var(t(boot.std)))
# apply(outer(apply(abs(yo), 1, max), abs(Z.wald), '>'),2,mean) 
