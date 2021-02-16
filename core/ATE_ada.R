source("core/phi.R")
source("core/psi.R")
source("core/phi-cf.R")
source("core/psi-cf.R")
source("core/pred-cf.R")
source("core/measure_v1.R")
source("core/utility.R")
# source("core/gshrinkv4.R")
# source("core/bshrinkv4.R")



ATE.ada = function(surv, D, Z,tol=1e-7,maxit=100,
                      maxATE = log(2)*(.Machine$double.max.exp-1),
                      nfold = 5,
                      penalties = NULL)
{  
  
  n = nrow(Z)
  p = ncol(Z)
  
  # Zm = apply(Z, 2, mean)
  # Zsd = apply(Z, 2, sd)
  # 
  # Z = (Z-matrix(rep(Zm,each=n),n,p))/matrix(rep(Zsd,each=n),n,p)
  
  # Nuisance parameters
  #======================
  #DZ = scale(as.matrix(cbind(D,Z)))    
  # Ada Lasso
  
  if(is.null(penalties))
  {
    beta.ridge = tune.ahazpen(surv,as.matrix(cbind(D,Z)),penalty=lasso.control(alpha=.Machine$double.eps),
                              lambda = exp(log(10)*seq(2,-2,length.out = 100)))
    beta.adawgt =  1/abs(as.numeric(
      predict(beta.ridge,type="coef",lambda = beta.ridge$lambda.min)))
    
    cv.bhat.nopen = tune.ahazpen(surv,as.matrix(cbind(D,Z)),penalty="lasso",
                                 penalty.wgt = c(0,beta.adawgt[-1]))
    theta.nopen =  as.numeric(
      predict(cv.bhat.nopen,type="coef",lambda = cv.bhat.nopen$lambda.min))[1]
   
    # hat beta star
    # cv.bhat = tune.ahazpen(surv,as.matrix(cbind(D,Z)),penalty="lasso",
    #                              penalty.wgt = beta.adawgt)
    # hthetabeta =  as.numeric(
    #   predict(cv.bhat,type="coef",lambda = cv.bhat$lambda.min))
    hthetabeta =   as.numeric(
      predict(cv.bhat.nopen,type="coef",lambda = cv.bhat.nopen$lambda.min))
    
    # sum(hthetabeta!=0)
    
    # hat gamma
    gr.ridge = cv.glmnet(Z,D,family = "binomial", alpha = 0)
    gr.adawgt = 1/abs(coef(gr.ridge, s = "lambda.min")[-1])
    
    cv.hgr = cv.glmnet(Z,D,family = "binomial", penalty.factor = gr.adawgt)
    hgr = as.numeric(
      predict(cv.hgr,type="coef",s = "lambda.min"))
  }else
  {
    cv.bhat.nopen = ahazpen(surv,as.matrix(cbind(D,Z)),penalty="lasso",
                                 penalty.wgt = c(0,penalties$beta.adawgt[-1]))
    theta.nopen =  as.numeric(
      predict(cv.bhat.nopen,type="coef",lambda = penalties$lambda.beta.nopen))[1]
    
    # cv.bhat = tune.ahazpen(surv,as.matrix(cbind(D,Z)),penalty="lasso",
    #                        penalty.wgt = penalties$beta.adawgt)
    # hthetabeta =  as.numeric(
    #   predict(cv.bhat,type="coef",lambda = penalties$lambda.bhat))
    hthetabeta =  as.numeric(
      predict(cv.bhat.nopen,type="coef",lambda = penalties$lambda.beta.nopen))
    
    cv.hgr = glmnet(Z,D,family = "binomial", penalty.factor = penalties$gr.adawgt)
    hgr = as.numeric(
      predict(cv.hgr,type="coef",s = penalties$lambda.hgr))
  }
  
  if(is.null(penalties))
  {
    penalties = list(
      lambda.betaridge = beta.ridge$lambda.min,
      beta.adawgt = beta.adawgt,
      # lambda.beta.nopen = cv.bhat.nopen$lambda.min,
  #    lambda.bhat = cv.bhat$lambda.min,
      lambda.grridge = gr.ridge$lambda.min,
      gr.adawgt = gr.adawgt
  # ,lambda.hgr = cv.hgr$lambda.min
  )
    }
  
  if( is.null(penalties$lambda.hgr))
  {
    cv.hgr = cv.glmnet(Z,D,family = "binomial", penalty.factor = penalties$gr.adawgt)
    hgr = as.numeric(
      predict(cv.hgr,type="coef",s = "lambda.min"))
    penalties$lambda.hgr = cv.hgr$lambda.min
  }else{
    cv.hgr = glmnet(Z,D,family = "binomial", penalty.factor = penalties$gr.adawgt)
    hgr = as.numeric(
      predict(cv.hgr,type="coef",s = penalties$lambda.hgr))
  }
  
  if( is.null(penalties$lambda.beta.nopen))
  {
    cv.bhat.nopen = tune.ahazpen(surv,as.matrix(cbind(D,Z)),penalty="lasso",
                                 penalty.wgt = c(0,penalties$beta.adawgt[-1]))
    theta.nopen =  as.numeric(
      predict(cv.bhat.nopen,type="coef",lambda = cv.bhat.nopen$lambda.min))[1]
    penalties$lambda.beta.nopen = cv.bhat.nopen$lambda.min
    
    # hat beta star
    # cv.bhat = tune.ahazpen(surv,as.matrix(cbind(D,Z)),penalty="lasso",
    #                              penalty.wgt = beta.adawgt)
    # hthetabeta =  as.numeric(
    #   predict(cv.bhat,type="coef",lambda = cv.bhat$lambda.min))
    hthetabeta =   as.numeric(
      predict(cv.bhat.nopen,type="coef",lambda = cv.bhat.nopen$lambda.min))
  }else{
    cv.bhat.nopen = ahazpen(surv,as.matrix(cbind(D,Z)),penalty="lasso",
                            penalty.wgt = c(0,penalties$beta.adawgt[-1]))
    theta.nopen =  as.numeric(
      predict(cv.bhat.nopen,type="coef",lambda = penalties$lambda.beta.nopen))[1]
    
    # cv.bhat = tune.ahazpen(surv,as.matrix(cbind(D,Z)),penalty="lasso",
    #                        penalty.wgt = penalties$beta.adawgt)
    # hthetabeta =  as.numeric(
    #   predict(cv.bhat,type="coef",lambda = penalties$lambda.bhat))
    hthetabeta =  as.numeric(
      predict(cv.bhat.nopen,type="coef",lambda = penalties$lambda.beta.nopen))
  }
  
  

  
  foldid = rep(1:nfold,ceiling(n/nfold))[1:n]
  
  hgr.cf = hthetabeta.cf = matrix(0,p+1,nfold)
  
  for(foldk in 1:nfold)
  {
    fitid = (foldid!=foldk) #& (foldid.tune!=k)
    Dk = D[fitid]
    Zk = Z[fitid,]
    survk = surv[fitid,]
    
    cv.bhat.nopen = ahazpen(survk,as.matrix(cbind(Dk,Zk)),penalty="lasso",
                            penalty.wgt = c(0,penalties$beta.adawgt[-1]))
    
    # cv.bhat = tune.ahazpen(survk,as.matrix(cbind(D,Z)),penalty="lasso",
    #                        penalty.wgt = penalties$beta.adawgt)
    # hthetabeta =  as.numeric(
    #   predict(cv.bhat,type="coef",lambda = penalties$lambda.bhat))
    hthetabeta.cf[,foldk] =  as.numeric(
      predict(cv.bhat.nopen,type="coef",lambda = penalties$lambda.beta.nopen))
    
    cv.hgr = cv.glmnet(Zk,Dk,family = "binomial", penalty.factor = penalties$gr.adawgt)
    hgr.cf[,foldk] = as.numeric(
      predict(cv.hgr,type="coef",s = penalties$lambda.hgr))
  }
  
  
  ps = 1/(1+exp(-hgr[1]-drop(Z %*% hgr[-1])))
  
  ipw = D/ps + (1-D)/(1-ps)
  ipw.stab = D*ipw*sum(D)/sum(D*ipw) + (1-D)*ipw*sum(1-D)/sum((1-D)*ipw) 
  ipw.stab = pmax(0.1,pmin(10,ipw.stab))
  
  ipw.fit = ahaz(surv,D,weights = ipw.stab)
  
  fit = list(hthetabeta = hthetabeta, hgr=hgr)
  # ATE 
  #===============
  
  # No cross-fitting --------------------------------
  
  
  dr = list()
  # hat : phi with lasso beta and lasso gamma
  pD = 1/(1+exp(-drop(Z%*%hgr[-1])-
                  hgr[1]))
  hazbZ = drop(Z%*%hthetabeta[-1])
  tmpfit = phi(hazbZ,hthetabeta[1],pD,
               D,surv,maxit,tol)
  dr$est = tmpfit$ate
  dr$sd = tmpfit$ate.sd
  
  
  # checkhat : psi with lasso beta and lasso gamma
  hdi0 =  psi(rep(0,n),pD,
               D,surv)$ate
  
  # checkhat : psi with lasso beta and lasso gamma
  hdi =  psi(hazbZ,pD,
             D,surv)$ate
  
  # Cross-fitting ----------------------------------
  pD = ps.cf(foldid, hgr.cf,Z)
  hazbZ = hazbZ.cf(foldid, 
                   hthetabeta.cf[-1,],Z)
  dr.cf = phi.cf(foldid,
                  hthetabeta.cf[1,],hazbZ, pD,
                  D,surv,maxit,tol)$ate
  
  # check : psi with regularized lasso beta and regularized lasso gamma
  hdi.cf = psi.cf(foldid,hazbZ, pD,
                  D,surv)$ate

  
  out = list(fit=fit,or = theta.nopen, ipw=coef(ipw.fit),
             dr = dr,
             hdi = hdi, hdi0 = hdi0, 
             dr.cf = dr.cf,
             hdi.cf = hdi.cf,
             penalties = penalties)
  
  return(out)
}