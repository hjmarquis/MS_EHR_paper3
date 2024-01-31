
expit = function(x) 1/(1+exp(-x))
discretize.ave = function(x, nbin=5)
{
  if(length(unique(x))< nbin)
  {
    return(x)
  }
  grp = cut(x, unique(c(-Inf,quantile(x, probs = (1:(nbin-1))/nbin), Inf)))
  for (lev in levels(grp))
  {
    x[grp==lev] = mean(x[grp==lev])
  }
  return(x)
}

discretize = function(x, nbin=5)
{
  if(length(unique(x))< nbin)
  {
    return(as.numeric(factor(x)))
  }
  grp = cut(x, unique(c(-Inf,quantile(x, probs = (1:(nbin-1))/nbin), Inf)))
  return(as.numeric(grp))
}

dat.ehr.list = paste0("comb-data-",c("RN","DF"),"_MStrt")
dat.reg.list = paste("reg", "-data-",
                     c("RN","DF"), "_MStrt",sep = "")
pair.list = c("RN","DF")
# RN 
#==========================================================

for(k in 1:2)
{
  
  # Load CLIMB fit
  #---------------------------------------------------------
  
  load( paste("MS CLIME analysis/data-analysis/",dat.reg.list[k],".rda",sep=''))
  load(paste("MS CLIME analysis/result/MS_result_new_",dat.reg.list[k],".rda",sep=''))
  
  Z = Z[,c("FEMALE", "RACE", "AGE_AT_FIRSTMSICD", "FOLLOWUP_DURA", "DISEASE_DURA",
           "PRIOR_RELAPSE_12MONS", "PRIOR_RELAPSE_24MONS")]
  
  C.ps.ahaz = expit(drop(Z %*% fit.alltil17$fit$hgr[-1])+fit.alltil17$fit$hgr[1])
  C.or.ahaz = drop(Z %*% fit.alltil17$fit$hthetabeta[-1])
  
  C.ps.Y1  = expit(drop(Z %*% fit.1yr$or1[-1]) + fit.1yr$or1[1])
  C.or.Y1  = (expit(drop(Z %*% fit.1yr$or1[-1]) + fit.1yr$or1[1]) * C.ps.Y1 + 
                expit(drop(Z %*% fit.1yr$or0[-1]) + fit.1yr$or0[1]) * (1-C.ps.Y1))
  
  C.ps.Y2  = expit(drop(Z %*% fit.2yr$or1[-1]) + fit.2yr$or1[1])
  C.or.Y2  = (expit(drop(Z %*% fit.2yr$or1[-1]) + fit.2yr$or1[1]) * C.ps.Y2 + 
                expit(drop(Z %*% fit.2yr$or0[-1]) + fit.2yr$or0[1]) * (1-C.ps.Y2))
  
  # Load EHR fit
  #---------------------------------------------------------
  
  id.keep = id
  load( paste("MS CLIME analysis/data-analysis/",dat.ehr.list[k],"p10_yr.rda",sep=''))
  load(paste("MS CLIME analysis/result/MS_result_new_",dat.ehr.list[k],"p10.rda",sep=''))
  Z = Z[id %in% id.keep,]
  
  
  U.or.ahaz = exp(drop(Z %*% fit.alltil17$fit$hthetabeta[-1]))
  U.ps.ahaz = expit(drop(Z %*% fit.alltil17$fit$hgr[-1]) + fit.alltil17$fit$hgr[1])
  
  U.ps.Y1 = expit(drop(Z %*% fit1$ps[-1]) + fit1$ps[1])
  U.or1.Y1 = expit(drop(Z %*% fit1$or1[-1]) + fit1$or1[1])
  U.or0.Y1 = expit(drop(Z %*% fit1$or0[-1]) + fit1$or0[1])
  
  U.ps.Y2 =expit( drop(Z %*% fit2$ps[-1]) + fit2$ps[1])
  U.or1.Y2 = expit(drop(Z %*% fit2$or1[-1]) + fit2$or1[1])
  U.or0.Y2 = expit(drop(Z %*% fit2$or0[-1]) + fit2$or0[1])
  
  
  # Sensitivity
  #---------------------------------------------------------
  n = length(D)
  nbin = 10
  # ahaz
  C1 = discretize(C.ps.ahaz, nbin = nbin)
  C2 = discretize(C.or.ahaz, nbin = nbin)
  
  RReu = RRud = matrix(1, max(C1),max(C2))
  
  C.ps.ahaz = discretize.ave(C.ps.ahaz, nbin = nbin)
  C.or.ahaz = discretize.ave(C.or.ahaz, nbin = nbin)
  U.ps.ahaz = discretize.ave(U.ps.ahaz, nbin = nbin)
  U.or.ahaz = discretize.ave(U.or.ahaz, nbin = nbin)
  
  for (i in 1:max(C1)) 
  {
    for (j in 1:max(C2))
    {
      cell.pos = (C1==i) & (C2==j)
      if(all(!cell.pos))
      {
        next
      }
      
      RReu[i,j] = max(U.ps.ahaz[cell.pos]*(1-C.ps.ahaz[cell.pos])/(1-U.ps.ahaz[cell.pos])/C.ps.ahaz[cell.pos])
      RRud[i,j] = max(U.or.ahaz[cell.pos])/min(U.or.ahaz[cell.pos])
    }
  }
  
  BF.ahaz = RReu*RRud/(RReu + RRud - 1)
  colnames(BF.ahaz) = sort(unique(C.or.ahaz))
  row.names(BF.ahaz) = sort(unique(C.ps.ahaz))
  ave.BF.ahaz = sum(BF.ahaz * table(C1,C2))/n
  BF.ahaz[which(table(C1,C2)==0)] = NA
  
  # Y1
  
  C1 = discretize(C.ps.Y1, nbin = nbin)
  C2 = discretize(C.or.Y1, nbin = nbin)
  
  RReu = RRud = matrix(1, max(C1),max(C2))
  
  C.ps.Y1 = discretize.ave(C.ps.Y1, nbin = nbin)
  C.or.Y1 = discretize.ave(C.or.Y1, nbin = nbin)
  U.ps.Y1 = discretize.ave(U.ps.Y1, nbin = nbin)
  U.or1.Y1 = discretize.ave(U.or1.Y1, nbin = nbin)
  U.or0.Y1 = discretize.ave(U.or0.Y1, nbin = nbin)
  
  for (i in 1:max(C1)) 
  {
    for (j in 1:max(C2))
    {
      cell.pos = (C1==i) & (C2==j)
      if(all(!cell.pos))
      {
        next
      }
      RReu[i,j] = max(U.ps.Y1[cell.pos]*(1-C.ps.Y1[cell.pos])/(1-U.ps.Y1[cell.pos])/C.ps.Y1[cell.pos])
      RRud[i,j] = max(max(U.or1.Y1[cell.pos])/min(U.or1.Y1[cell.pos]),
                      max(U.or0.Y1[cell.pos])/min(U.or0.Y1[cell.pos]))
    }
  }
  
  BF.Y1 = RReu*RRud/(RReu + RRud - 1)
  colnames(BF.Y1) = sort(unique(C.or.Y1))
  row.names(BF.Y1) = sort(unique(C.ps.Y1))
  ave.BF.Y1 = sum(BF.Y1 * table(C1,C2))/n
  BF.Y1[which(table(C1,C2)==0)] = NA
  
  # Y2
  
  C1 = discretize(C.ps.Y2, nbin = nbin)
  C2 = discretize(C.or.Y2, nbin = nbin)
  
  RReu = RRud = matrix(1, max(C1),max(C2))
  
  C.ps.Y2 = discretize.ave(C.ps.Y2, nbin = nbin)
  C.or.Y2 = discretize.ave(C.or.Y2, nbin = nbin)
  U.ps.Y2 = discretize.ave(U.ps.Y2, nbin = nbin)
  U.or1.Y2 = discretize.ave(U.or1.Y2, nbin = nbin)
  U.or0.Y2 = discretize.ave(U.or0.Y2, nbin = nbin)
  
  for (i in 1:max(C1)) 
  {
    for (j in 1:max(C2))
    {
      cell.pos = (C1==i) & (C2==j)
      if(all(!cell.pos))
      {
        next
      }
      RReu[i,j] = max(U.ps.Y2[cell.pos]*(1-C.ps.Y2[cell.pos])/(1-U.ps.Y2[cell.pos])/C.ps.Y2[cell.pos])
      RRud[i,j] = max(max(U.or1.Y2[cell.pos])/min(U.or1.Y2[cell.pos]),
                      max(U.or0.Y2[cell.pos])/min(U.or0.Y2[cell.pos]))
    }
  }
  
  BF.Y2 = RReu*RRud/(RReu + RRud - 1)
  colnames(BF.Y2) = sort(unique(C.or.Y2))
  row.names(BF.Y2) = sort(unique(C.ps.Y2))
  ave.BF.Y2 = sum(BF.Y2 * table(C1,C2))/n
  BF.Y2[which(table(C1,C2)==0)] = NA
  
  save(BF.ahaz, BF.Y1, BF.Y2, ave.BF.ahaz, ave.BF.Y1, ave.BF.Y2, 
       file = paste0("MS CLIME analysis/result/heat_map_data_", pair.list[k], ".rda"))
  
}