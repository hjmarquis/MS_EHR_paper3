
rm(list = objects())
library(glmnet)
source("MS CLIME analysis/ate_dich.R")

dat.name.list = paste0("comb-data-",c("RN","DF"),"_MStrt")

B = 1000

dat.name = dat.name.list[1]


# for(dat.name in dat.name.list)
load( paste("MS CLIME analysis/data-analysis/",dat.name,"p10_yr.rda",sep=''))
load(paste("MS CLIME analysis/result/MS_result_new_",dat.name,"p10.rda",sep=''))


# Two year

or1.var = colnames(Z)[fit2$or1[-1]!=0]
or0.var = colnames(Z)[fit2$or0[-1]!=0]
ps.var = colnames(Z)[fit2$ps[-1]!=0]

or1.min.nz = apply(Z[,or1.var], 2, function(x) min(table(year.grp, x!=0)[,"TRUE"]))
or0.min.nz = apply(Z[,or0.var], 2, function(x) min(table(year.grp, x!=0)[,"TRUE"]))
ps.min.nz = apply(Z[,ps.var], 2, function(x) min(table(year.grp, x!=0)[,"TRUE"]))

yr.col = model.matrix(~year.grp)[,-1]

offY1 = fit2$or1[1] + drop( Z %*% fit2$or1[-1])
Z1.0608 = Z[,or1.var[or1.min.nz>=20]] * (year.grp=="2006-2008")
colnames(Z1.0608) = paste0(colnames(Z1.0608), "_0608")
Z1.0911 = Z[,or1.var[or1.min.nz>=20]] * (year.grp=="2009-2011")
colnames(Z1.0911) = paste0(colnames(Z1.0911), "_0911")
Z1.1216 = Z[,or1.var[or1.min.nz>=20]] * (year.grp=="2012-2016")
colnames(Z1.1216) = paste0(colnames(Z1.1216), "_1216")
Z1 = cbind(yr.col, Z1.0608, Z1.0911, Z1.1216)

offY0 =fit2$or0[1] + drop(Z %*% fit2$or0[-1])
Z0.0608 = Z[,or0.var[or0.min.nz>=20]] * (year.grp=="2006-2008")
colnames(Z0.0608) = paste0(colnames(Z0.0608), "_0608")
Z0.0911 = Z[,or0.var[or0.min.nz>=20]] * (year.grp=="2009-2011")
colnames(Z0.0911) = paste0(colnames(Z0.0911), "_0911")
Z0.1216 = Z[,or0.var[or0.min.nz>=20]] * (year.grp=="2012-2016")
colnames(Z0.1216) = paste0(colnames(Z0.1216), "_1216")
Z0 = cbind(yr.col, Z0.0608, Z0.0911, Z0.1216)

offD = fit2$ps[1]+drop(Z %*% fit2$ps[-1])
ZD.0608 = Z[,ps.var[ps.min.nz>=20]] * (year.grp=="2006-2008")
colnames(ZD.0608) = paste0(colnames(ZD.0608), "_0608")
ZD.0911 = Z[,ps.var[ps.min.nz>=20]] * (year.grp=="2009-2011")
colnames(ZD.0911) = paste0(colnames(ZD.0911), "_0911")
ZD.1216 = Z[,ps.var[ps.min.nz>=20]] * (year.grp=="2012-2016")
colnames(ZD.1216) = paste0(colnames(ZD.1216), "_1216")
ZD = cbind(yr.col, ZD.0608, ZD.0911, ZD.1216)
  
if(file.exists(paste("MS CLIME analysis/result/MS_result_new_",dat.name,"p10_yr_dich2.rda",sep='')))
{
  load(paste("MS CLIME analysis/result/MS_result_new_",dat.name,"p10_yr_dich2.rda",sep=''))
}else
{
  fit2.time.adj = ate.dich.ada.off(Y2, D, Z1, Z0, ZD, 
                                   offY1, offY0, offD)
  
  c(fit2$dr, fit2.time.adj$dr)
  
  save(fit2.time.adj, file = paste("MS CLIME analysis/result/MS_result_new_",dat.name,"p10_yr_dich2.rda",sep=''))
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
    
    bootfit = try(ate.dich.ada(Y2[bootid],
                               D[bootid],Z[bootid,],
                               fit2$or1.pen, fit2$or0.pen, fit2$ps.pen,
                               fit2$or1.lambda, fit2$or0.lambda, fit2$ps.lambda),TRUE)
    if(inherits(bootfit,"try-error"))
    {
      next
    }
    
    off1 = drop(Z[bootid,] %*% bootfit$or1[-1])+bootfit$or1[1]
    off0 = drop(Z[bootid,] %*% bootfit$or0[-1])+bootfit$or0[1]
    offD = drop(Z[bootid,] %*% bootfit$ps[-1]) + bootfit$ps[1]
    
    bootfit.time.adj = try(ate.dich.ada.off(Y2[bootid], D[bootid], 
                                            Z1[bootid,,drop=F], 
                                            Z0[bootid,,drop=F], 
                                            ZD[bootid,,drop=F], 
                                            offY1, offY0, offD,
                                            fit2.time.adj$or1.pen, 
                                            fit2.time.adj$or0.pen, 
                                            fit2.time.adj$ps.pen,
                                            fit2.time.adj$or1.lambda, 
                                            fit2.time.adj$or0.lambda, 
                                            fit2.time.adj$ps.lambda),TRUE)
    cont = inherits(bootfit.time.adj,"try-error")
  }
  run.time = Sys.time() - start.time
  print(run.time)
  
  boot.time.adj[b] = bootfit.time.adj$dr
}
save(boot.time.adj, file = paste0("MS CLIME analysis/result/MS_result_new_",dat.name,"p10_yr_dich2_boot.rda"))


