
library(ahaz)
library(glmnet)
# 
# source("core/ATE_ada.R")
# source("MS CLIME analysis/ate_dich.R")


dat.name.list = paste0("comb-data-",c("RN","DF"),"_MStrt")

B = 1000

dat.name = dat.name.list[2]
load( paste("MS CLIME analysis/data-analysis/",dat.name,"p10_yr.rda",sep=''))

est = c(mean(Y1[D==1]) - mean(Y1[D==0]),
        mean(Y2[D==1]) - mean(Y2[D==0]),
        coef(ahaz(surv,D)))

boot.est = matrix(NA, 3, B)
for (b in 1:B)
{
  bootid = sample(1:length(D),replace = T)
  boot.surv = Surv(jitter(surv[bootid,1], 1e-7),surv[bootid,2])
  bootD = D[bootid]
  bootY1 = Y1[bootid]
  bootY2 = Y2[bootid]
  
  boot.est[,b] = c(mean(bootY1[bootD==1]) - mean(bootY1[bootD==0]),
                   mean(bootY2[bootD==1]) - mean(bootY2[bootD==0]),
                  coef(ahaz(boot.surv,bootD)))
}

apply(boot.est, 1, quantile, probs = c(.025, .975))
apply(sign(boot.est)!= sign(est), 1, mean)*2

sig =apply(apply(boot.est,1,quantile, probs = c(0.025, 0.975)),
           2,diff)/(1.96*2)
est = apply(boot.est, 1, median)
Z.wald = abs(est/sig)
boot.std = (boot.est - est) / sig

apply(outer(apply(abs(boot.std), 2, max), Z.wald, '>'),2,mean)


fit = ahaz(surv, D)
summary(fit)
exp(coef(fit)*2)
exp((coef(fit)+ summary(fit)$coef[,"Std. Error"]*1.96*c(-1,1))*2)

