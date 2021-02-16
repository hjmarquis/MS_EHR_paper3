dat.name.list = paste(rep(c("comb"),2), "-data-",
                      rep(c("RN","DF"),each = 2), '_',
                      rep(c("MStrt","RxNorm"),2),sep = "")

library(weights)
library(sjstats)



bal.table = vector("list",length(dat.name.list))
nlist = ntrtlist = rep(0,length(dat.name.list))

for (i in 1:length(dat.name.list)) {
  # print(i)
  rm(list = setdiff(objects(), 
                    c("PS.bal.table","nlist","ntrtlist","dat.name.list","i",
                      "bal.table","max.char")))
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
                   "Ever Hospitalization within 3 month", 
                   "Ever overall Emergency", 
                   "Ever Emergency within 3 month", 
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
  
  
  select = which(fit.alltil17$fit$hgr[-1]!=0)
  Z = Z[,select]
  varICD = which(is.element(select,varICD))
  varCUI = which(is.element(select,varCUI))
  varCPT = which(is.element(select,varCPT))
  varall = which(is.element(select,varall))
  var3m = which(is.element(select,var3m))
  descr.labels = descr.labels[is.element(varlowd, select)]
  varlowd = which(is.element(select,varlowd))
  p = ncol(Z)
  
  var.binary = match(c("FEMALE","RACE","ADJ_HOSP_OVERALL","ADJ_HOSP_3MONS" ,     
                       "ADJ_ED_OVERALL" , "ADJ_ED_3MONS"),colnames(Z))
  var.binary = var.binary[!is.na(var.binary)]
  var.real = setdiff(1:p, var.binary)
  
  Z.binary = matrix(Z[,var.binary],n)
  descr.count = cbind(apply(Z.binary[D==0,],2,sum), 
                      apply(Z.binary[D==1,],2,sum),
                      apply(Z.binary,2,sum))
  descr.rate =apply(format(descr.count * rep(1/c(n-ntrt,ntrt,n)*100,
                                             each =length(var.binary)),digits = 0,nsmall=1, scientific=F,
                           trim = F),2,paste,"%")
  descr.count = format(descr.count,trim=F)
  descr.count.test = format(apply(Z.binary,2, function(x)
    chisq.test(table(D,x))$p.value),
    digits = 0, nsmall = 2)
  
  descr.mean = format(cbind(apply(Z[D==0,],2,mean), 
                            apply(Z[D==1,],2,mean),
                            apply(Z,2,mean)),
                      digits = 0, nsmall = 2,trim=F)
  descr.sd = format(cbind(apply(Z[D==0,],2,sd), 
                          apply(Z[D==1,],2,sd), 
                          apply(Z,2,sd)),
                    digits = 0, nsmall = 2,trim=F)
  descr.test = format(apply(Z,2,function(x)
    wilcox.test(x[D==0],x[D==1])$p.value),
    digits = 0, nsmall = 2)
  descr.test[descr.test == "0.00"] = "<0.01"
  
  varICDALL =intersect(varall,varICD)
  varICD3M = intersect(var3m,varICD)
  varCPTALL = intersect(varall,varCPT)
  varCPT3M = intersect(var3m,varCPT)
  varCUIALL = intersect(varall,varCUI)
  varCUI3M = intersect(var3m,varCUI)
  
  descr.labels.select = rep("",p)
  descr.labels.select[varlowd] = descr.labels
  descr.labels.select[varICDALL] = paste("Phecode",
                                         gsub("_", ".",gsub("_OVERALL", '', 
                                                            gsub("PheCode.", "",colnames(Z)[varICDALL]))),
                                         "overall")
  descr.labels.select[varICD3M] = paste("Phecode",
                                        gsub("_", ".",gsub("_3MONS", '', 
                                                           gsub("PheCode.", "",colnames(Z)[varICD3M]))),
                                        "within 3 months")
  descr.labels.select[varCPTALL] = paste("CPT",
                                         gsub("_", ".",gsub("_OVERALL", '', 
                                                            gsub("CPTGroup.", "",colnames(Z)[varCPTALL]))),
                                         "overall")
  descr.labels.select[varCPT3M] = paste("CPT",
                                        gsub("_", ".",gsub("_3MONS", '', 
                                                           gsub("CPTGroup.", "",colnames(Z)[varCPT3M]))),
                                        "within 3 months")
  descr.labels.select[varCUIALL] = paste("CUI",
                                         gsub("_", ".",gsub("_OVERALL", '', 
                                                            gsub("CUI.", "",colnames(Z)[varCUIALL]))),
                                         "overall")
  descr.labels.select[varCUI3M] = paste("CUI",
                                        gsub("_", ".",gsub("_3MONS", '', 
                                                           gsub("CUI.", "",colnames(Z)[varCUI3M]))),
                                        "within 3 months")
  
  
  
  bal.table[[i]] = data.frame(var = descr.labels.select, 
                              matrix(paste(descr.mean, ' (', descr.sd, ')',
                                    sep = ''),ncol = 3), p = descr.test)
  bal.table[[i]][var.binary,2:4] = matrix(paste(descr.count, ' (', descr.rate, ')',
                                               sep = ''),ncol = 3)
  bal.table[[i]][var.binary,'p'] = descr.count.test
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
  
  # ps = 1/(1+exp(-fit.alltil17$fit$hgr[1]-drop(Z%*%fit.alltil17$fit$hgr[-1])))
  # ipw = D/ps + (1-D)/(1-ps)
  # ipw.stab = D*ipw*sum(D)/sum(D*ipw) + (1-D)*ipw*sum(1-D)/sum((1-D)*ipw) 
  # ipw.stab = pmax(0.1,pmin(10,ipw.stab))
  # ipw1 = ipw.stab[D==1] 
  # ipw0 = ipw.stab[D==0]
  
  ps = 1/(1+exp(-fit.alltil17$fit$hgr[1]-drop(Z%*%fit.alltil17$fit$hgr[-1])))
  ps = pmin(0.9,pmax(0.1,ps))
  ipw = D/ps + (1-D)/(1-ps)
  ipw.stab = D*ipw*sum(D)/sum(D*ipw) + (1-D)*ipw*sum(1-D)/sum((1-D)*ipw) 
  ipw1 = ipw.stab[D==1] 
  ipw0 = ipw.stab[D==0]
  
  # description table
  
  # Z = Z*rep(Zsd,each=n) + rep(Zm,each=n)
  
  p = ncol(Z)
  
  varICD = grep("PheCode", colnames(Z))
  varCUI = grep("CUI\\.", colnames(Z))
  varCPT = grep("CPT", colnames(Z))
  varall = grep("OVERALL", colnames(Z))
  var3m = grep("3MONS", colnames(Z))
  varlowd = setdiff(1:p,c(varICD,varCUI,varCPT))
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
                   "Ever Hospitalization within 3 month", 
                   "Ever overall Emergency", 
                   "Ever Emergency within 3 month", 
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
  
  
  select = which(fit.alltil17$fit$hgr[-1]!=0)
  Z = Z[,select]
  varICD = which(is.element(select,varICD))
  varCUI = which(is.element(select,varCUI))
  varCPT = which(is.element(select,varCPT))
  varall = which(is.element(select,varall))
  var3m = which(is.element(select,var3m))
  descr.labels = descr.labels[is.element(varlowd, select)]
  varlowd = which(is.element(select,varlowd))
  p = ncol(Z)
  
  var.binary = match(c("FEMALE","RACE","ADJ_HOSP_OVERALL","ADJ_HOSP_3MONS" ,     
                       "ADJ_ED_OVERALL" , "ADJ_ED_3MONS"),colnames(Z))
  var.binary = var.binary[!is.na(var.binary)]
  var.real = setdiff(1:p, var.binary)
  
  Z.binary = matrix(Z[,var.binary],n)
  
  ipw.count = cbind(apply(ipw0*Z.binary[D==0,],2,sum), 
                    apply(ipw1*Z.binary[D==1,],2,sum),
                    apply(Z.binary,2,sum))
  ipw.rate =apply(format(ipw.count * rep(1/c(n-ntrt,ntrt,n)*100,
                                         each =length(var.binary)),digits = 0,nsmall=1, scientific=F, trim = FALSE),2,paste,"%")
  ipw.count = format(ipw.count,digits = 0, nsmall = 0, trim = FALSE)
  ipw.count.test =  format(
    apply(Z.binary,2, function(x)
      wtd.chi.sq(D,x,weight = ipw.stab)["p.value"]),
    digits = 0, nsmall = 2, trim = FALSE)
  
  ipw.mean = format(cbind(apply(ipw0*Z[D==0,],2,mean), 
                          apply(ipw1*Z[D==1,],2,mean),
                          apply(Z,2,mean)),
                    digits = 0, nsmall = 2, trim = FALSE)
  ipw.sd = format(cbind(sqrt(apply(ipw0*Z[D==0,]^2,2,mean)-
                               apply(ipw0*Z[D==0,],2,mean)^2), 
                        sqrt(apply(ipw1*Z[D==1,]^2,2,mean)-
                               apply(ipw1*Z[D==1,],2,mean)^2), 
                        apply(Z,2,sd)),
                  digits = 0, nsmall = 2, trim = FALSE)
  tmp.dat = data.frame(Z=Z[,3],D=factor(D), ipw.stab = ipw.stab)
  ipw.test = format( apply(Z, 2, function(x)
    weighted_mannwhitney(Z ~ D + ipw.stab, data.frame(Z=x,D=factor(D), ipw.stab = ipw.stab))$p.value),
    digits = 0, nsmall = 2)
  ipw.test[ipw.test == "0.00"] = "<0.01"
  
  varICDALL =intersect(varall,varICD)
  varICD3M = intersect(var3m,varICD)
  varCPTALL = intersect(varall,varCPT)
  varCPT3M = intersect(var3m,varCPT)
  varCUIALL = intersect(varall,varCUI)
  varCUI3M = intersect(var3m,varCUI)
  
  descr.labels.select = rep("",p)
  descr.labels.select[varlowd] = descr.labels
  descr.labels.select[varICDALL] = paste("Phecode",
                                         gsub("_", ".",gsub("_OVERALL", '', 
                                                            gsub("PheCode.", "",colnames(Z)[varICDALL]))),
                                         "overall")
  descr.labels.select[varICD3M] = paste("Phecode",
                                        gsub("_", ".",gsub("_3MONS", '', 
                                                           gsub("PheCode.", "",colnames(Z)[varICD3M]))),
                                        "within 3 months")
  descr.labels.select[varCPTALL] = paste("CPT",
                                         gsub("_", ".",gsub("_OVERALL", '', 
                                                            gsub("CPTGroup.", "",colnames(Z)[varCPTALL]))),
                                         "overall")
  descr.labels.select[varCPT3M] = paste("CPT",
                                        gsub("_", ".",gsub("_3MONS", '', 
                                                           gsub("CPTGroup.", "",colnames(Z)[varCPT3M]))),
                                        "within 3 months")
  descr.labels.select[varCUIALL] = paste("CUI",
                                         gsub("_", ".",gsub("_OVERALL", '', 
                                                            gsub("CUI.", "",colnames(Z)[varCUIALL]))),
                                         "overall")
  descr.labels.select[varCUI3M] = paste("CUI",
                                        gsub("_", ".",gsub("_3MONS", '', 
                                                           gsub("CUI.", "",colnames(Z)[varCUI3M]))),
                                        "within 3 months")
  # PS.bal.table[[i]] = cbind(descr.labels.select,matrix(paste(
  #   rbind(ipw.count,ipw.mean
  #         # ,ICDALL.mean,ICD3M.mean,
  #         # CPTALL.mean,CPT3M.mean, CUIALL.mean,CUI3M.mean
  #   )," (",
  #   rbind(ipw.rate,ipw.sd
  #         # ,ICDALL.sd,ICD3M.sd,
  #         # CPTALL.sd,CPT3M.sd, CUIALL.sd,CUI3M.sd
  #   ), ")",sep=''),ncol=3),
  #   c(ipw.count.test,ipw.test))
  PS.bal.table[[i]] = data.frame(var = descr.labels.select, 
                              matrix(paste(ipw.mean, ' (', ipw.sd, ')',
                                           sep = ''),ncol = 3), p = ipw.test)
  PS.bal.table[[i]][var.binary,2:4] = matrix(paste(ipw.count, ' (', ipw.rate, ')',
                                                sep = ''),ncol = 3)
  PS.bal.table[[i]][var.binary,'p'] = ipw.count.test

}

for(i in 1:length(dat.name.list))
{
  tmp.tab = cbind(bal.table[[i]][,c(1,4,2,3,5)],
                         PS.bal.table[[i]][,c(2,3,5)])
  # tmp.tab = gsub("\\\\|\\$", "", tmp.tab)
  write.csv(tmp.tab, file = paste("MS CLIME analysis/report/Supp-Table-2-",dat.name.list[i],".csv", sep=''),
            row.names = F)
}