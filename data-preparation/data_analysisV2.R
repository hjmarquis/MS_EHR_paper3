
library(survival)
library(ggplot2)

data.name.list = paste(rep(c("low-d","high-d","comb"),4), "-data-",
                                        rep(c("RN","DF"),each = 6), '_',
                                        rep(c("MStrt","RxNorm"),2,each=3),sep = "")
for(data.name in data.name.list)
{
  datdir = paste('MS CLIME analysis/data-processed/', data.name,'.rds', sep='')
  dat <- readRDS(datdir)


atrisk17 = dat$RELAPSE_DATE >= as.Date("2016-12-31")
dat$RELAPSE[atrisk17] = 0 
dat$RELAPSE_DATE = pmin(dat$RELAPSE_DATE, as.Date("2016-12-31"))
before17 = dat$start_date <= as.Date("2016-12-31")

dat = dat[before17,]


surv = Surv((dat$RELAPSE_DATE-dat$start_date + runif(nrow(dat),0.1))/365.2475,dat$RELAPSE)
D = as.numeric(factor(dat$medication_desc))-1

table(apply(sapply(dat,is.nan),2,any))

Z = as.matrix(dat[,-1:-10])
Z[is.nan(Z)] = 0

if(length(grep("RN",data.name))==1)
{
  CUI.trt = c("C1172734", "C1529600", "C0393022", "C4047978", "C0732355")
}else{
  CUI.trt = c("C1699926", "C0388087", "C2938762", "C0058218", "C3556178")
}

trt.col = unlist(sapply(CUI.trt, grep,x=colnames(Z)))

if(length(trt.col)>0)
{
  Z = Z[,-trt.col]
}

nonEHR.col = match(c("FEMALE","RACE", "AGE_AT_FIRSTMSICD","FOLLOWUP_DURA","DISEASE_DURA",
               "PRIORDMT_DURA", "PRIOR_RELAPSE_12MONS", "PRIOR_RELAPSE_24MONS"),
               colnames(Z))
dich.col = match(c("ADJ_HOSP_OVERALL","ADJ_HOSP_3MONS","ADJ_ED_OVERALL","ADJ_ED_3MONS"),
                 colnames(Z))

nonEHR.col = nonEHR.col[!is.na(nonEHR.col)]
dich.col = dich.col[!is.na(dich.col)]
EHR.col = setdiff(1:ncol(Z), c(nonEHR.col,dich.col))

var.overall = grep("OVERALL", colnames(Z))
var.3m = setdiff(grep("MONS", colnames(Z)), c(nonEHR.col,dich.col))

# Z[,var.overall] = log(1+Z[,var.overall])
# Z[,var.3m] = as.numeric(Z[,var.3m]!= 0) 

Z[,EHR.col] = log(1+Z[,EHR.col])
if(length(dich.col)>0)
{
  Z[,dich.col] = as.numeric(Z[,dich.col]>0)
}

clinical.col = setdiff(1:ncol(Z), c(grep("PheCode\\.",colnames(Z)),
                                    grep("CPTGroup\\.",colnames(Z)),
                                    grep("CUI\\.",colnames(Z))))

min.count = round(nrow(Z)*0.1)

keep.col = which(apply(Z!=0,2,sum)>=min.count)
# keep.col = sort(unique(c(which(apply(Z!=0,2,sum)>=min.count), clinical.col)))
Z = Z[,keep.col]

Zm = apply(Z, 2, mean)
Zsd = apply(Z, 2, sd)
Zmax = apply(Z, 2, max) - apply(Z, 2, min)


n = nrow(Z)
p = ncol(Z)

Z = (Z-matrix(rep(Zm,each=n),n))/matrix(rep(Zmax,each=n),n)

save(surv,D,Z, Zm,Zsd,Zmax,
     file = paste("MS CLIME analysis/data-analysis/",data.name,"p10_max.rda",sep=''))
}
