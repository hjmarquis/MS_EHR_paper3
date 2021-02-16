dat.name.list = paste(rep(c("low-d"),2), "-data-",
                      rep(c("RN","DF"),each = 2), '_',
                      rep(c("MStrt","RxNorm"),2),sep = "")

library(weights)
library(sjstats)

bal.table = vector("list",length(dat.name.list))
nlist = ntrtlist = rep(0,length(dat.name.list))

for (i in 1:length(dat.name.list)) {
  rm(list = setdiff(objects(), 
                    c("PS.bal.table","nlist","ntrtlist","dat.name.list","i",
                      "bal.table")))
  load(paste("MS CLIME analysis/data-analysis/",dat.name.list[i],"p10.rda",sep=''))
  load( paste("MS CLIME analysis/result/MS_result_new_",dat.name.list[i],"p10.rda",sep=''))
  
  nlist[i] = n = length(D)
  p = ncol(Z)
  d = sum(surv[,2])
  ntrtlist[i]=ntrt = sum(D)
  
  
  # description table
  
  # Z = Z*rep(Zsd,each=n) + rep(Zm,each=n)
  
  
  varICD = grep("PheCode", colnames(Z))
  varCUI = grep("CUI\\.", colnames(Z))
  varCPT = grep("CPT", colnames(Z))
  varall = grep("OVERALL", colnames(Z))
  var3m = grep("3MONS", colnames(Z))
  varlowd = setdiff(1:p,c(varICD,varCUI,varCPT))
  
  
  var.binary = 1:2
  var.real = setdiff(varlowd, var.binary)
  
  descr.count = cbind(apply(Z[D==0,1:2],2,sum), 
                      apply(Z[D==1,1:2],2,sum),
                      apply(Z[,1:2],2,sum))
  descr.rate =apply(format(descr.count * rep(1/c(n-ntrt,ntrt,n)*100,each =2),digits = 0,nsmall=1, scientific=F,trim=TRUE),2,paste,"\\%")
  descr.count = format(descr.count,trim=TRUE)
  descr.count.test = format(c(chisq.test(table(D,Z[,1]))$p.value,
                              chisq.test(table(D,Z[,2]))$p.value),
                            digits = 0, nsmall = 2)
  
  descr.mean = format(cbind(apply(Z[D==0,var.real],2,mean), 
                            apply(Z[D==1,var.real],2,mean),
                            apply(Z[,var.real],2,mean)),
                      digits = 0, nsmall = 2,trim=TRUE)
  descr.sd = format(cbind(apply(Z[D==0,var.real],2,sd), 
                          apply(Z[D==1,var.real],2,sd), 
                          apply(Z[,var.real],2,sd)),
                    digits = 0, nsmall = 2,trim=TRUE)
  descr.test = format(apply(Z[,var.real],2,function(x)
    wilcox.test(x[D==0],x[D==1])$p.value),
    digits = 0, nsmall = 2)
  descr.test[descr.test == "0.00"] = "$<0.01$"
  
  descr.labels = c("Female", 
                   "Non-Hispanic white",  
                   "Age at first MS ICD", 
                   "Follow up duration", 
                   "Disease duration", 
                   "Overall healthcare utilizations", 
                   "Healthcare utilizations within 3 month", 
                   "Adjusted overall MS ICD", 
                   "Adjusted MS ICD within 3 month", 
                   "Adjusted MS CUI within 3 month", 
                   "Adjusted overall Steroid", 
                   "Adjusted Steroid within 3 month", 
                   "Adjusted overall MRI", 
                   "Adjusted MRI within 6 month", 
                   "Adjusted overall Hospitalization", 
                   "Adjusted Hospitalization within 3 month", 
                   "Adjusted overall Emergency", 
                   "Adjusted Emergency within 3 month", 
                   "Months on prior DMT",
                   "Num relapse prior 1 year",
                   "Num relapse prior 2 years"
  )
  clinical.var = c("FEMALE", "RACE", "AGE_AT_FIRSTMSICD","FOLLOWUP_DURA","DISEASE_DURA",
                   "HUTIL_OVERALL","HUTIL_3MONS","ADJ_MSICD_OVERALL","ADJ_MSICD_3MONS","ADJ_MSCUI_3MONS",
                   "ADJ_STEROID_OVERALL","ADJ_STEROID_6MONS","ADJ_MRI_OVERALL","ADJ_MRI_6MONS","ADJ_HOSP_OVERALL",
                   "ADJ_HOSP_3MONS","ADJ_ED_OVERALL","ADJ_ED_3MONS", "PRIORDMT_DURA", "PRIOR_RELAPSE_12MONS",
                   "PRIOR_RELAPSE_24MONS")
  descr.labels = descr.labels[match(colnames(Z)[varlowd],clinical.var)]
  
  # hline.pos = length(varlowd)
  
  bal.table[[i]] = cbind(descr.labels,matrix(paste(
    rbind(descr.count,descr.mean
          # ,ICDALL.mean,ICD3M.mean,
          # CPTALL.mean,CPT3M.mean, CUIALL.mean,CUI3M.mean
          )," (",
    rbind(descr.rate,descr.sd
          # ,ICDALL.sd,ICD3M.sd,
          # CPTALL.sd,CPT3M.sd, CUIALL.sd,CUI3M.sd
          ), ")",sep=''),ncol=3),
    c(descr.count.test,descr.test))
  # 
  # nlab = length(descr.labels)
  # descr.order = c(outer(c(0,nlab),1:nlab,"+"))
  # descr.labels = c(paste("\\multirow{2}{1.5in}{",
  #                        descr.labels,"}", sep="")
  #                  , rep("",nlab))[descr.order]
  # descr.stat = rbind(descr.count,descr.mean, 
  #                    matrix(paste("(",rbind(descr.rate,descr.sd), ")",sep=''),
  #                           ncol = 3))[descr.order, ]
  # descr.p = c(descr.count.test,descr.test,rep("",nlab))[descr.order]
  # 
  # 
  # bal.table[[i]] = cbind(descr.labels,descr.stat,descr.p)
}

PS.bal.table = vector("list",length(dat.name.list))
# nlist = ntrtlist = rep(0,length(dat.name.list))

for (i in 1:length(dat.name.list)) 
{
  rm(list = setdiff(objects(), 
                    c("PS.bal.table","nlist","ntrtlist","dat.name.list","i",
                      "bal.table")))
  load(paste("MS CLIME analysis/data-analysis/",dat.name.list[i],"p10.rda",sep=''))
  load( paste("MS CLIME analysis/result/MS_result_new_",dat.name.list[i],"p10.rda",sep=''))
  
  
  
  # nlist[i] = 
  n = length(D)
  p = ncol(Z)
  d = sum(surv[,2])
  # ntrtlist[i]=
  ntrt = sum(D)
  
  ps = 1/(1+exp(-fit.alltil17$fit$hgr[1]-drop(Z%*%fit.alltil17$fit$hgr[-1])))
  ipw = D/ps + (1-D)/(1-ps)
  ipw.stab = D*ipw*sum(D)/sum(D*ipw) + (1-D)*ipw*sum(1-D)/sum((1-D)*ipw) 
  ipw.stab = pmax(0.1,pmin(10,ipw.stab))
  ipw1 = ipw.stab[D==1] 
  ipw0 = ipw.stab[D==0]
  
  # description table
  
  # Z = Z*rep(Zsd,each=n) + rep(Zm,each=n)
  
  
  varICD = grep("PheCode", colnames(Z))
  varCUI = grep("CUI\\.", colnames(Z))
  varCPT = grep("CPT", colnames(Z))
  varall = grep("OVERALL", colnames(Z))
  var3m = grep("3MONS", colnames(Z))
  varlowd = setdiff(1:p,c(varICD,varCUI,varCPT))
  
  
  var.binary = 1:2
  var.real = setdiff(varlowd, var.binary)
  
  ipw.count = cbind(apply(ipw0*Z[D==0,1:2],2,sum), 
                    apply(ipw1*Z[D==1,1:2],2,sum),
                    apply(Z[,1:2],2,sum))
  ipw.rate =apply(format(ipw.count * rep(1/c(n-ntrt,ntrt,n)*100,each =2),digits = 0,nsmall=1, scientific=F,trim=TRUE),2,paste,"\\%")
  ipw.count = format(ipw.count,trim=TRUE,digits = 0,nsmall=0)
  ipw.count.test =  format(c(wtd.chi.sq(D,Z[,1],weight = ipw.stab)["p.value"],
                             wtd.chi.sq(D,Z[,2],weight = ipw.stab)["p.value"]),
                           digits = 0, nsmall = 2,trim=TRUE)
  
  ipw.mean = format(cbind(apply(ipw0*Z[D==0,var.real],2,mean), 
                          apply(ipw1*Z[D==1,var.real],2,mean),
                          apply(Z[,var.real],2,mean)),
                    digits = 0, nsmall = 2,trim=TRUE)
  ipw.sd = format(cbind(sqrt(apply(ipw0*Z[D==0,var.real]^2,2,mean)-
                               apply(ipw0*Z[D==0,var.real],2,mean)^2), 
                        sqrt(apply(ipw1*Z[D==1,var.real]^2,2,mean)-
                               apply(ipw1*Z[D==1,var.real],2,mean)^2), 
                        apply(Z[,var.real],2,sd)),
                  digits = 0, nsmall = 2,trim=TRUE)
  tmp.dat = data.frame(Z=Z[,3],D=factor(D), ipw.stab = ipw.stab)
  ipw.test = format( apply(Z[,var.real], 2, function(x)
    weighted_mannwhitney(Z ~ D + ipw.stab, data.frame(Z=x,D=factor(D), ipw.stab = ipw.stab))$p.value),
    digits = 0, nsmall = 2)
  ipw.test[ipw.test == "0.00"] = "$<0.01$"
  
  descr.labels = c("Female", 
                   "Non-Hispanic white",  
                   "Age at first MS ICD", 
                   "Follow up duration", 
                   "Disease duration", 
                   "Overall healthcare utilizations", 
                   "Healthcare utilizations within 3 month", 
                   "Adjusted overall MS ICD", 
                   "Adjusted MS ICD within 3 month", 
                   "Adjusted MS CUI within 3 month", 
                   "Adjusted overall Steroid", 
                   "Adjusted Steroid within 3 month", 
                   "Adjusted overall MRI", 
                   "Adjusted MRI within 6 month", 
                   "Adjusted overall Hospitalization", 
                   "Adjusted Hospitalization within 3 month", 
                   "Adjusted overall Emergency", 
                   "Adjusted Emergency within 3 month", 
                   "Months on prior DMT",
                   "Num relapse prior 1 year",
                   "Num relapse prior 2 years"
  )
  clinical.var = c("FEMALE", "RACE", "AGE_AT_FIRSTMSICD","FOLLOWUP_DURA","DISEASE_DURA",
                   "HUTIL_OVERALL","HUTIL_3MONS","ADJ_MSICD_OVERALL","ADJ_MSICD_3MONS","ADJ_MSCUI_3MONS",
                   "ADJ_STEROID_OVERALL","ADJ_STEROID_6MONS","ADJ_MRI_OVERALL","ADJ_MRI_6MONS","ADJ_HOSP_OVERALL",
                   "ADJ_HOSP_3MONS","ADJ_ED_OVERALL","ADJ_ED_3MONS", "PRIORDMT_DURA", "PRIOR_RELAPSE_12MONS",
                   "PRIOR_RELAPSE_24MONS")
  descr.labels = descr.labels[match(colnames(Z)[varlowd],clinical.var)]
  
  # hline.pos = length(varlowd)
  
  PS.bal.table[[i]] = cbind(descr.labels,matrix(paste(
    rbind(ipw.count,ipw.mean
          # ,ICDALL.mean,ICD3M.mean,
          # CPTALL.mean,CPT3M.mean, CUIALL.mean,CUI3M.mean
          )," (",
    rbind(ipw.rate,ipw.sd
          # ,ICDALL.sd,ICD3M.sd,
          # CPTALL.sd,CPT3M.sd, CUIALL.sd,CUI3M.sd
          ), ")",sep=''),ncol=3),
    c(ipw.count.test,ipw.test))
  
  # nlab = length(descr.labels)
  # descr.order = c(outer(c(0,nlab),1:nlab,"+"))
  # descr.labels = c(descr.labels, rep("",nlab))[descr.order]
  # descr.stat = rbind(ipw.count,ipw.mean, 
  #                    matrix(paste("(",rbind(ipw.rate,ipw.sd), ")",sep=''),
  #                           ncol = 3))[descr.order, ]
  # descr.p = c(ipw.count.test,ipw.test,rep("",nlab))[descr.order]
  # 
  # 
  # PS.bal.table[[i]] = cbind(descr.labels,descr.stat,descr.p)
  
}

for(i in 1:length(dat.name.list))
{
  bal.table[[i]] = cbind(bal.table[[i]][,c(1,4,2,3,5)],
                         PS.bal.table[[i]][,c(2,3,5)])
  
  tmp.tab = bal.table[[i]]
  tmp.tab = gsub("\\\\|\\$", "", tmp.tab)
  write.csv(tmp.tab, file = paste("MS CLIME analysis/report/Supp-Table-1-",dat.name.list[i],".csv", sep=''),
            row.names = F)
}