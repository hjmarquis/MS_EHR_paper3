
rm(list = objects())
# .libPaths(c("~/R-3.6.1/library",.libPaths()))
# np = 20
library(ahaz)
library(glmnet)
# library(doParallel)
# library(doRNG)
# cl = makeCluster(getOption("cl.cores", np))
# registerDoParallel(cl)
# registerDoRNG(seed = 531)
# clusterEvalQ(cl, .libPaths(c("~/R-3.6.1/library",.libPaths())))
source("core/ATE_ada_offset.R")
source("core/ATE_ada.R")


dat.name.list = paste0("comb-data-",c("RN","DF"),"_MStrt")

B = 1000


dat.name = dat.name.list[2]


# for(dat.name in dat.name.list)
load( paste("MS CLIME analysis/data-analysis/",dat.name,"p10_yr.rda",sep=''))
load(paste("MS CLIME analysis/result/MS_result_new_",dat.name,"p10.rda",sep=''))

surv = surv[year>=2013,]
D = D[year>=2013]
Z = Z[year>=2013,]

or.var = colnames(Z)[fit.alltil17$fit$hthetabeta[-1]!=0]
ps.var = colnames(Z)[fit.alltil17$fit$hgr[-1]!=0]

or.min.nz = apply(Z[,or.var]!=0, 2, sum)
ps.min.nz = apply(Z[,ps.var]!=0, 2, sum)

offs = drop(Z %*% fit.alltil17$fit$hthetabeta[-1])
offD = drop(Z %*% fit.alltil17$fit$hgr[-1])
Zs = Z[,or.var]
ZD = Z[,ps.var]

if(file.exists(paste("MS CLIME analysis/result/MS_result_new_",dat.name,"p10_yr.rda",sep='')))
{
  load(paste("MS CLIME analysis/result/MS_result_new_",dat.name,"p10_yr.rda",sep=''))
}else
{
  fit.time.adj = ATE.ada.off(surv, D, ZD, Zs, offD, offs)
  
  exp(c(fit.time.adj$dr$est, fit.alltil17$dr$est)*2)
  
  save(fit.time.adj, file = paste("MS CLIME analysis/result/MS_result_new_",dat.name,"p10_yr.rda",sep=''))
}

boot.time.adj = rep(NA, B)
for (b in 1:B) 
{
  print(b)
  start.time = Sys.time()
  cont = TRUE
  
  while(cont)
  {
    bootid = sample(1:length(D),replace = T)
    
    boot.surv = Surv(jitter(surv[bootid,1], 1e-8),surv[bootid,2])
    bootfit = try(ATE.ada(boot.surv,
                          D[bootid],Z[bootid,],
                          penalties =  fit.alltil17$penalties),TRUE)
    if(inherits(bootfit,"try-error"))
    {
      next
    }
    
    offs = drop(Z[bootid,] %*% bootfit$fit$hthetabeta[-1])
    offD = drop(Z[bootid,] %*% bootfit$fit$hgr[-1]) + bootfit$fit$hgr[1]
    
    bootfit.time.adj = try(ATE.ada.off(boot.surv, 
                                       D[bootid],
                                       ZD[bootid,], Zs[bootid,], offD, offs,
                                       penalties = fit.time.adj$penalties),TRUE)
    cont = inherits(bootfit.time.adj,"try-error")
  }
  run.time = Sys.time() - start.time
  print(run.time)
  
  boot.time.adj[b] = bootfit.time.adj$dr$est
}

save(boot.time.adj, file = paste0("MS CLIME analysis/result/MS_result_new_",dat.name,"p10_yr_ahaz_boot.rda"))
# 
# rm(list = objects())
