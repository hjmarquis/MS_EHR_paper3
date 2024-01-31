
rm(list = objects())
# .libPaths(c("~/R-3.6.1/library",.libPaths()))
# np = 20
library(glmnet)
# library(doParallel)
# library(doRNG)
# cl = makeCluster(getOption("cl.cores", np))
# registerDoParallel(cl)
# registerDoRNG(seed = 531)
# clusterEvalQ(cl, .libPaths(c("~/R-3.6.1/library",.libPaths())))
source("MS CLIME analysis/ate_dich.R")


dat.name.list = paste0("comb-data-",c("RN","DF"),"_MStrt")

B = 1000


dat.name = dat.name.list[2]


# for(dat.name in dat.name.list)
load( paste("MS CLIME analysis/data-analysis/",dat.name,"p10_yr.rda",sep=''))
load(paste("MS CLIME analysis/result/MS_result_new_",dat.name,"p10.rda",sep=''))


D = D[year>=2013]
Z = Z[year>=2013,]
Y1 = Y1[year>=2013]

or1.var = colnames(Z)[fit1$or1[-1]!=0]
or0.var = colnames(Z)[fit1$or0[-1]!=0]
ps.var = colnames(Z)[fit1$ps[-1]!=0]


offY1 = fit1$or1[1] + drop( Z %*% fit1$or1[-1])
offY0 = fit1$or0[1] + drop( Z %*% fit1$or0[-1])
offD = fit1$ps[1] + drop( Z %*% fit1$ps[-1])

if(length(or1.var)>0)
{
  or1.min.nz = apply(Z[,or1.var,drop=F]!=0, 2, sum)
  Z1 = Z[,or1.var[or1.min.nz>20], drop = F]
}else{
  Z1 = Z[, NULL]
}

if(length(or0.var)>0)
{
  or0.min.nz = apply(Z[,or0.var,drop=F]!=0, 2, sum)
  Z0 = Z[,or0.var[or0.min.nz>20], drop = F]
}else{
  Z0 = Z[, NULL]
}

if(length(ps.var)>0)
{
  ps.min.nz = apply(Z[,ps.var,drop=F]!=0, 2, sum)
  ZD = Z[,ps.var[ps.min.nz>20], drop = F]
}else{
  ZD = Z[, NULL]
}

if(file.exists(paste("MS CLIME analysis/result/MS_result_new_",dat.name,"p10_yr_dich1.rda",sep='')))
{
  load(paste("MS CLIME analysis/result/MS_result_new_",dat.name,"p10_yr_dich1.rda",sep=''))
}else
{
  fit1.time.adj = ate.dich.ada.off(Y1, D, Z1, Z0, ZD, 
                                   offY1, offY0, offD)
  
  c(fit1$dr, fit1.time.adj$dr)
  
  save(fit1.time.adj, file = paste("MS CLIME analysis/result/MS_result_new_",dat.name,"p10_yr_dich1.rda",sep=''))
}


boot.time.adj = rep(NA , B)
for(b in 1:B)
{
  print(b)
  start.time = Sys.time()
  cont = TRUE
  
  while(cont)
  {
    bootid = sample(1:length(D),replace = T)
    
    bootfit = try(ate.dich.ada(Y1[bootid],
                               D[bootid],Z[bootid,],
                               fit1$or1.pen, fit1$or0.pen, fit1$ps.pen,
                               fit1$or1.lambda, fit1$or0.lambda, fit1$ps.lambda),TRUE)
    if(inherits(bootfit,"try-error"))
    {
      next
    }
    
    off1 = drop(Z[bootid,] %*% bootfit$or1[-1])+bootfit$or1[1]
    off0 = drop(Z[bootid,] %*% bootfit$or0[-1])+bootfit$or0[1]
    offD = drop(Z[bootid,] %*% bootfit$ps[-1]) + bootfit$ps[1]
    
    bootfit.time.adj = try(ate.dich.ada.off(Y1[bootid], D[bootid], 
                                            Z1[bootid,,drop=F], 
                                            Z0[bootid,,drop=F], 
                                            ZD[bootid,,drop=F], 
                                            offY1, offY0, offD,
                                            fit1.time.adj$or1.pen, 
                                            fit1.time.adj$or0.pen, 
                                            fit1.time.adj$ps.pen,
                                            fit1.time.adj$or1.lambda, 
                                            fit1.time.adj$or0.lambda, 
                                            fit1.time.adj$ps.lambda),TRUE)
    cont = inherits(bootfit.time.adj,"try-error")
  }
  run.time = Sys.time() - start.time
  print(run.time)
  
  boot.time.adj[b] = bootfit.time.adj$dr
}
save(boot.time.adj, file = paste0("MS CLIME analysis/result/MS_result_new_",dat.name,"p10_yr_dich1_boot.rda"))


