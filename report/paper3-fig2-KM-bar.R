library(survival)
library(ggplot2)
library(glmnet)
library(ahaz)
source("core/basehaz.R")

data.name.list = paste(c("low-d","comb"), "-data-",
                       "RN", '_',
                       "MStrt",sep = "")
title.list = c("Low-dimensional",
               "High-dimensional")

plist = list()
for(i in 1:length(data.name.list))
{
  data.name = data.name.list[i]
  load(paste("MS CLIME analysis/data-analysis/",data.name,"p10.rda",sep=''))

  # Crude plot
  #-----------------------------------------------------------------
  
  time.grid = unique(sort(pmin(surv[,1],3)))
  
  KMfit1 = summary(survfit(surv[D==1]~1),time.grid)
  KMfit0 = summary(survfit(surv[D==0]~1),time.grid)
  # KMfit1$surv = KMfit1$surv[KMfit1$time<=3]
  # KMfit0$surv = KMfit0$surv[KMfit0$time<=3]
  # KMfit1$time = KMfit1$time[KMfit1$time<=3]
  # KMfit0$time = KMfit0$time[KMfit0$time<=3]
  
  if(length(grep("data-RN",data.name))==1)
  {
    trt.level = c("Natalizumab","Rituximab")
  }else{
    trt.level = c("Dimethyl Fumarate", "Fingolimod")
  }
  
  ggdf = data.frame(time = c(0,KMfit0$time,
                             0, KMfit1$time,
                             1, 1, 2,2),
                    relapse1 = c(0,1-KMfit0$surv,
                                 0,1-KMfit1$surv,
                                 rep(NA,4)),
                    relapse2 = c(rep(NA,length(KMfit0$time)+length(KMfit1$time)+2),
                                 mean((surv[D==0,1] <= 1)*surv[D==0,2]),
                                 mean((surv[D==1,1] <= 1)*surv[D==1,2]),
                                 mean((surv[D==0,1] <= 2)*surv[D==0,2]),
                                 mean((surv[D==1,1] <= 2)*surv[D==1,2])),
                    trt = c(rep(trt.level, c(length(KMfit0$time),length(KMfit1$time))+1), 
                            rep(trt.level,2)), 
                    type = rep(c("curve","bar"),
                               c(length(KMfit0$time)+length(KMfit1$time)+2, 4)))
  
  
  p = ggplot(ggdf, aes(x = time, y = relapse1, group = trt,color = trt))
  p = p + geom_step()
  p = p + geom_bar( aes(x = time, y = relapse2, group = trt,fill = trt),
                    stat="identity",width = 0.2,position=position_dodge())
  p = p +  labs(title = paste("Crude"),
                fill = "Medication", color =  "Medication")
  p = p + ylab("Relapse Probability") + xlab("Year(s) since Initiation of Medication")
   
  
  plist = c(plist,list(p))
  
  
  # IPW plot
  #-----------------------------------------------------------------
  # hat gamma
  gr.ridge = cv.glmnet(Z,D,family = "binomial", alpha = 0)
  gr.adawgt = 1/abs(coef(gr.ridge, s = "lambda.min")[-1])
  
  cv.hgr = cv.glmnet(Z,D,family = "binomial", penalty.factor = gr.adawgt)
  hgr = as.numeric(
    predict(cv.hgr,type="coef",s = "lambda.min"))
  
  ps = 1/(1+exp(-hgr[1]-drop(Z %*% hgr[-1])))
  
  ipw = D*mean(D)/ps + (1-D)*(1-mean(D))/(1-ps)
  ipw = pmax(0.1,pmin(10,ipw))
  ipw1 = ipw[D==1]/sum(ipw[D==1])
  ipw0 = ipw[D==0]/sum(ipw[D==0])
  
  KMfit1 = summary(survfit(surv[D==1]~1, weights = ipw1),time.grid)
  KMfit0 = summary(survfit(surv[D==0]~1, weights = ipw0),time.grid)
  # KMfit1$surv = KMfit1$surv[KMfit1$time<=3]
  # KMfit0$surv = KMfit0$surv[KMfit0$time<=3]
  # KMfit1$time = KMfit1$time[KMfit1$time<=3]
  # KMfit0$time = KMfit0$time[KMfit0$time<=3]
  
  if(length(grep("data-RN",data.name))==1)
  {
    trt.level = c("Natalizumab","Rituximab")
  }else{
    trt.level = c("Dimethyl Fumarate", "Fingolimod")
  }
  
  
  
  ggdf = data.frame(time = c(0,KMfit0$time,
                             0, KMfit1$time,
                             1, 1, 2,2),
                    relapse1 = c(0,1-KMfit0$surv,
                                 0,1-KMfit1$surv,
                                 rep(NA,4)),
                    relapse2 = c(rep(NA,length(KMfit0$time)+length(KMfit1$time)+2),
                                 sum((surv[D==0,1] <= 1)*surv[D==0,2]*ipw0),
                                 sum((surv[D==1,1] <= 1)*surv[D==1,2]*ipw1),
                                 sum((surv[D==0,1] <= 2)*surv[D==0,2]*ipw0),
                                 sum((surv[D==1,1] <= 2)*surv[D==1,2]*ipw1)),
                    trt = c(rep(trt.level, c(length(KMfit0$time),length(KMfit1$time))+1), 
                            rep(trt.level,2)), 
                    type = rep(c("curve","bar"),
                               c(length(KMfit0$time)+length(KMfit1$time)+2, 4)))
  
  
  p = ggplot(ggdf, aes(x = time, y = relapse1, group = trt,color = trt))
  p = p + geom_step()
  p = p + geom_bar( aes(x = time, y = relapse2, group = trt,fill = trt),
                    stat="identity",width = 0.2,position=position_dodge())
  p = p +  labs(title = paste(title.list[i],": IPW"),
                fill = "Medication", color =  "Medication")
  p = p + ylab("Relapse Probability") + xlab("Year(s) since Initiation of Medication")
  p = p + theme(axis.title.x=element_text(size=rel(0.8)))
   
  
  plist = c(plist,list(p))
  
  # OR plot
  #-----------------------------------------------------------------
  
  # Ada Lasso
  beta.ridge = tune.ahazpen(surv,as.matrix(cbind(D,Z)),penalty=lasso.control(alpha=1e-7))
  beta.adawgt =  1/abs(as.numeric(
    predict(beta.ridge,type="coef",lambda = beta.ridge$lambda.min))[-1])
  
  # hat beta star
  cv.bhat = tune.ahazpen(surv,as.matrix(cbind(D,Z)),penalty="lasso",
                         penalty.wgt = c(0,beta.adawgt))
  hthetabeta =  as.numeric(
    predict(cv.bhat,type="coef",lambda = cv.bhat$lambda.min))
  
  lp = D*hthetabeta[1] + drop(Z %*% hthetabeta[-1])
  Haz0 = basehaz.ahaz(surv,lp)
  
  
  x.grid = sort(unique(c(seq(0,max(surv[,1]), length.out = 100), knots(Haz0))))
  x.grid = x.grid[x.grid<=3]
  
  lp0 = drop(Z %*% hthetabeta[-1])
  prob0 = apply(1-exp(-outer(x.grid,lp0,'*')-Haz0(x.grid)),1,mean)
  prob1 = apply(1-exp(-outer(x.grid,lp0+hthetabeta[1],'*')-Haz0(x.grid)),1,mean)
  
  # OR year 1
  Y = (surv[,1]<=1)*surv[,2]
  gr.ridge = cv.glmnet(Z[D==1,],Y[D==1],family = "binomial", alpha = 0)
  gr.adawgt = 1/abs(coef(gr.ridge, s = "lambda.min")[-1])
  
  cv.hgr = cv.glmnet(Z[D==1,],Y[D==1],family = "binomial", penalty.factor = gr.adawgt)
  or1 = as.numeric(
    predict(cv.hgr,type="coef",s = "lambda.min"))
  
  gr.ridge = cv.glmnet(Z[D==0,],Y[D==0],family = "binomial", alpha = 0)
  gr.adawgt = 1/abs(coef(gr.ridge, s = "lambda.min")[-1])
  
  cv.hgr = cv.glmnet(Z[D==0,],Y[D==0],family = "binomial", penalty.factor = gr.adawgt)
  or0 = as.numeric(
    predict(cv.hgr,type="coef",s = "lambda.min"))
  
  prob1.1yr = mean(1/(1+exp(-drop(Z%*%or1[-1])-or1[1])))
  prob0.1yr = mean(1/(1+exp(-drop(Z%*%or0[-1])-or0[1])))
  
  # OR year 2
  Y = (surv[,1]<=2)*surv[,2]
  gr.ridge = cv.glmnet(Z[D==1,],Y[D==1],family = "binomial", alpha = 0)
  gr.adawgt = 1/abs(coef(gr.ridge, s = "lambda.min")[-1])
  
  cv.hgr = cv.glmnet(Z[D==1,],Y[D==1],family = "binomial", penalty.factor = gr.adawgt)
  or1 = as.numeric(
    predict(cv.hgr,type="coef",s = "lambda.min"))
  
  gr.ridge = cv.glmnet(Z[D==0,],Y[D==0],family = "binomial", alpha = 0)
  gr.adawgt = 1/abs(coef(gr.ridge, s = "lambda.min")[-1])
  
  cv.hgr = cv.glmnet(Z[D==0,],Y[D==0],family = "binomial", penalty.factor = gr.adawgt)
  or0 = as.numeric(
    predict(cv.hgr,type="coef",s = "lambda.min"))
  
  prob1.2yr = mean(1/(1+exp(-drop(Z%*%or1[-1])-or1[1])))
  prob0.2yr = mean(1/(1+exp(-drop(Z%*%or0[-1])-or0[1])))
  
  
  if(length(grep("data-RN",data.name))==1)
  {
    trt.level = c("Natalizumab","Rituximab")
  }else{
    trt.level = c("Dimethyl Fumarate", "Fingolimod")
  }
  
  ggdf = data.frame(time = c(x.grid,
                             x.grid,
                             1, 1, 2,2),
                    relapse1 = c(prob0,
                                 prob1,
                                 rep(NA,4)),
                    relapse2 = c(rep(NA,length(x.grid)*2),
                                 prob0.1yr,
                                 prob1.1yr,
                                 prob0.2yr,
                                 prob1.2yr),
                    trt = c(rep(trt.level, each = length(x.grid)), 
                            rep(trt.level,2)), 
                    type = rep(c("curve","bar"),
                               c(length(x.grid)*2, 4)))
  
  
  p = ggplot(ggdf, aes(x = time, y = relapse1, group = trt,color = trt))
  p = p + geom_step()
  p = p + geom_bar( aes(x = time, y = relapse2, group = trt,fill = trt),
                    stat="identity",width = 0.2,position=position_dodge())
  p = p +  labs(title = paste(title.list[i],": OR"),
                fill = "Medication", color =  "Medication")
  p = p + ylab("Relapse Probability") + xlab("Year(s) since Initiation of Medication")
  p = p + theme(axis.title.x=element_text(size=rel(0.8)))
   
  
  
  
  plist = c(plist,list(p))
}

library(gridExtra)
library(grid)

plist[[4]] = NULL


g <- ggplotGrob(plist[[1]] + theme(legend.position="bottom"))$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
lheight <- sum(legend$height)

pdf("MS CLIME analysis/figure/paper3-figure2.pdf", width = 10, height = 6)
grid.arrange(arrangeGrob(plist[[1]]+ theme(legend.position="none"),
            plist[[2]]+ theme(legend.position="none"),
            plist[[3]]+ theme(legend.position="none"),
            plist[[4]]+ theme(legend.position="none"),
            plist[[5]]+ theme(legend.position="none"),
            layout_matrix = cbind(c(1,1),c(2,3),c(4,5)),
            widths = c(1.5,1,1)
            ),
            legend,
            ncol = 1,
            heights = unit.c(unit(1, "npc") - lheight, lheight)
)
dev.off()
