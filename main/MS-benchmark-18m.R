library(ahaz)
library(twang)

version = "v1"
load(paste("MS CLIME analysis/data-analysis/dataHD_18m_",version,".rda",sep=''))

p = ncol(Z)
varICD = grep("phecode", colnames(Z))
varCUI = grep("CUI_", colnames(Z))
varlowd = setdiff(1:p,c(varICD,varCUI))

naive = ahaz(surv, D)
adj.lowd = ahaz(surv, cbind(D,Z[,varlowd]))

ps.lowd = predict(glm(D~as.matrix(Z[,varlowd]),family = binomial),type = "response")

fitfX= formula(paste('D', '~',
                     paste(colnames(Z)[varlowd],collapse = '+')))
dat.lowd = cbind(data.frame(D=D),Z[,varlowd])


ps.twang <- ps(fitfX, data=dat.lowd, print.level = 0,verbose=F,
               stop.method = "ks.max")
ipw.twang = D/ps.twang$ps$ks.max.ATE+(1-D)/(1-ps.twang$ps$ks.max.ATE)
ipw.twang[D==1] = ipw.twang[D==1]/mean(ipw.twang[D==1])
ipw.twang[D==0] = ipw.twang[D==0]/mean(ipw.twang[D==0])
ipw.twang = pmax(pmin(ipw.twang,10),0.1)

ipw.lowd = ahaz(surv,D,weights = ipw.twang,
                robust = TRUE)

load(paste("MS CLIME analysis/result/MS_HDresult_18m_",version,".rda",sep=''))

sb.nohrr = sum(fit.18m$fit$hthetabeta[-1]!=0)
sg.nohrr = sum(fit.18m$fit$hgr[-1]!=0)
ps.lasso = 1/(1+exp(-drop(Z%*%fit.18m$fit$hgr[-1])-fit.18m$fit$hgr[1]))
ipw.lasso = D/ps.lasso + (1-D)/(1-ps.lasso)
ipw.lasso[D==1] = ipw.lasso[D==1]/mean(ipw.lasso[D==1])
ipw.lasso[D==0] = ipw.lasso[D==0]/mean(ipw.lasso[D==0])
ipw.lasso = pmax(pmin(ipw.lasso,10),0.1)
ipw.lasso.fit = ahaz(surv, D, weights = ipw.lasso,
                 robust = TRUE)

save(naive, adj.lowd,
     ipw.lowd,ipw.lasso.fit,
     ps.lasso,
     file = paste("MS CLIME analysis/result/MS_benchmark_18m_",version,".rda",sep=''))
