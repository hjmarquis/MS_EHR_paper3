
expit = function(x)
{
  return(1/(1+exp(-x)))
}

# A simple adaptive lasso
adalasso = function(..., adawgt, lambda.min)
{
  if(missing(penalty.factor))
  {
    ridgefit = cv.glmnet(alpha = 0, ...)
    adawgt = 1/abs(coef(ridgefit, s = "lambda.min")[-1])
  }
  if(missing(lambda.min))
  {
    adafit = cv.glmnet(penalty.factor = adawgt, ...)
    lambda.min = ada.fit$lambda.min
  }else
  {
    adafit = glmnet(penalty.factor = adawgt, ...)
  }
  return(drop(coef(adafit, s = lambda.min)))
}

tune.adalasso = function(...)
{
  ridgefit = cv.glmnet(alpha = 0, ...)
  adawgt = 1/abs(coef(ridgefit, s = "lambda.min")[-1])
  adafit = cv.glmnet(penalty.factor = adawgt, ...)
  lambda.min = ada.fit$lambda.min
  
  return(list(adawgt=adawgt,lambda.min=lambda.min))
}

tune.OR.PS = function(Yi,Ti,Xi)
{
  or1.penalty = tune.adalasso(x=Xi[Ti==1,], 
                              y=Yi[Ti==1], 
                              family = "binomial")
  
  or0.penalty = tune.adalasso(x=Xi[Ti==0,], 
                              y=Yi[Ti==0], 
                              family = "binomial")
  ps.penalty = tune.adalasso(x=Xi, 
                             y=Ti, 
                             family = "binomial")
  
  return(list(or1 = or1.penalty,
              or0 = or0.penalty,
              ps = ps.penalty))
}

# The AIPW
AIPW = function(Yi,Ti,Xi, reg = TRUE, 
                stablize = TRUE, min.ipw = 0.1, max.ipw = 10,
                ps.adawgt, ps.lambda.min,
                or1.adawgt, or1.lambda.min,
                or0.adawgt, or0.lambda.min)
{
  if(reg)
  {
    OR1.coef = adalasso(x=Xi[Ti==1,], 
                        y=Yi[Ti==1], 
                        family = "binomial",
                        or1.adawgt, or1.lambda.min)
    OR0.coef = adalasso(x=Xi[Ti==0,], 
                        y=Yi[Ti==0], 
                        family = "binomial",
                        or0.adawgt, or0.lambda.min)
    PS.coef = adalasso(x=Xi, 
                       y=Ti, 
                       family = "binomial",
                       ps.adawgt, ps.lambda.min)
    
  }else
  {
    OR1.coef = coef(glm(Yi~Xi, family = binomial, subset = (Ti==1)))
    OR0.coef = coef(glm(Yi~Xi, family = binomial, subset = (Ti==0)))
    
    PS.coef = coef(glm(Ti~Xi, family = binomial))
  }
  
  pred.OR1 = expit(drop(Xi %*% OR1.coef[-1]) + OR1.coef[1])
  pred.OR0 = expit(drop(Xi %*% OR0.coef[-1]) + OR0.coef[1])
  
  pred.PS = expit(drop(Xi %*% PS.coef[-1]) + PS.coef[1])
  pred.IPW = Ti/pred.PS + (1-Ti)/(1-pred.PS)
  
  if(stablize)
  {
    pred.IPW[Ti == 1] = pred.IPW[Ti == 1]*sum(Ti==1)/sum(pred.IPW[Ti == 1])
    pred.IPW[Ti == 0] = pred.IPW[Ti == 0]*sum(Ti==0)/sum(pred.IPW[Ti == 0])
  }
  pred.IPW = pmin(max.ipw, pmax(min.ipw, pred.IPW))
  
  return(
    mean(pred.OR1 - pred.OR0 
         + Ti*(Yi-pred.OR1)*pred.IPW
         - (1-Ti)*(Yi-pred.OR0)*pred.IPW)
  )
}


# OR

OR = function(Yi,Ti,Xi, reg = TRUE,
              or1.adawgt, or1.lambda.min,
              or0.adawgt, or0.lambda.min)
{
  if(reg)
  {
    OR1.coef = adalasso(x=Xi[Ti==1,], 
                        y=Yi[Ti==1], 
                        family = "binomial",
                        or1.adawgt, or1.lambda.min)
    OR0.coef = adalasso(x=Xi[Ti==0,], 
                        y=Yi[Ti==0], 
                        family = "binomial",
                        or0.adawgt, or0.lambda.min)
      
  }else
  {
    OR1.coef = coef(glm(Yi~Xi, family = binomial, subset = (Ti==1)))
    OR0.coef = coef(glm(Yi~Xi, family = binomial, subset = (Ti==0)))
  }
  pred.OR1 = expit(drop(Xi %*% OR1.coef[-1]) + OR1.coef[1])
  pred.OR0 = expit(drop(Xi %*% OR0.coef[-1]) + OR0.coef[1])
  
  return(mean(pred.OR1-pred.OR0))
}

# IPW

IPW =  function(Yi,Ti,Xi, reg = TRUE, 
                stablize = TRUE, min.ipw = 0.1, max.ipw = 10,
                adawgt, lambda.min)
{
  if(reg)
  {
    PS.coef = adalasso(x=Xi, 
                       y=Ti, 
                       family = "binomial",
                       adawgt, lambda.min)
  }else
  {
    PS.coef = coef(glm(Ti~Xi, family = binomial))
  }
  pred.PS = expit(drop(Xi %*% PS.coef[-1]) + PS.coef[1])
  if(stablize)
  {
    pred.IPW[Ti == 1] = pred.IPW[Ti == 1]*sum(Ti==1)/sum(pred.IPW[Ti == 1])
    pred.IPW[Ti == 0] = pred.IPW[Ti == 0]*sum(Ti==0)/sum(pred.IPW[Ti == 0])
  }
  pred.IPW = pmin(max.ipw, pmax(min.ipw, pred.IPW))
  
  return(mean(Yi*Ti*pred.IPW - Yi*(1-Ti)*pred.IPW))
}

# The OR for joint population
OR.joint = function(Yi,Ti,Xi, Ri, reg = TRUE)
{
  if(reg)
  {
    OR1.coef = adalasso(x=Xi[Ti*Ri==1,], 
                        y=Yi[Ti*Ri==1], 
                        family = "binomial")
    OR0.coef = adalasso(x=Xi[(1-Ti)*Ri==1,], 
                        y=Yi[(1-Ti)*Ri==1], 
                        family = "binomial")
  }else
  {
    OR1.coef = coef(glm(Yi~Xi, family = binomial, subset = (Ti*Ri==1)))
    OR0.coef = coef(glm(Yi~Xi, family = binomial, subset = ((1-Ti)*Ri==1)))
  }
  
  pred.OR1 = expit(drop(Xi %*% OR1.coef[-1]) + OR1.coef[1])
  pred.OR0 = expit(drop(Xi %*% OR0.coef[-1]) + OR0.coef[1])
  
  return(mean(pred.OR1-pred.OR0))
}


# The IPW for the joint population
IPW.joint =  function(Yi,Ti,Xi, Ri, reg = TRUE)
{
  if(reg)
  {
    PST.coef = adalasso(x=Xi, 
                        y=Ti, 
                        family = "binomial")
    PSR1.coef = adalasso(x=Xi[Ti==1,], 
                         y=Ri[Ti==1], 
                         family = "binomial")
    PSR0.coef = adalasso(x=Xi[Ti==0,], 
                         y=Ri[Ti==0], 
                         family = "binomial")
  }else
  {
    PST.coef = coef(glm(Ti~Xi, family = binomial))
    PSR1.coef = coef(glm(Ri~Xi, family = binomial, subset = (Ti==1)))
    PSR0.coef = coef(glm(Ri~Xi, family = binomial, subset = (Ti==0)))
  }
  pred.PST = expit(drop(Xi %*% PST.coef[-1]) + PST.coef[1])
  pred.R1 = expit(drop(Xi %*% PSR1.coef[-1]) + PSR1.coef[1])
  pred.R0 = expit(drop(Xi %*% PSR0.coef[-1]) + PSR0.coef[1])
  
  return(mean((Yi*Ti/(pred.PST*pred.R1) - Yi*(1-Ti)/((1-pred.PST)*pred.R0))[Ri==1]))
}

# The AIPW for joint population
AIPW.joint =  function(Yi,Ti,Xi, Ri, reg = TRUE)
{
  if(reg)
  {
    OR1.coef = adalasso(x=Xi[Ti*Ri==1,], 
                        y=Yi[Ti*Ri==1], 
                        family = "binomial")
    OR0.coef = adalasso(x=Xi[(1-Ti)*Ri==1,], 
                        y=Yi[(1-Ti)*Ri==1], 
                        family = "binomial")
    PST.coef = adalasso(x=Xi, 
                       y=Ti, 
                       family = "binomial")
    PSR1.coef = adalasso(x=Xi[Ti==1,], 
                        y=Ri[Ti==1], 
                        family = "binomial")
    PSR0.coef = adalasso(x=Xi[Ti==0,], 
                         y=Ri[Ti==0], 
                         family = "binomial")
  }else
  {
    OR1.coef = coef(glm(Yi~Xi, family = binomial, subset = (Ti*Ri==1)))
    OR0.coef = coef(glm(Yi~Xi, family = binomial, subset = ((1-Ti)*Ri==1)))
    PST.coef = coef(glm(Ti~Xi, family = binomial))
    PSR1.coef = coef(glm(Ri~Xi, family = binomial, subset = (Ti==1)))
    PSR0.coef = coef(glm(Ri~Xi, family = binomial, subset = (Ti==0)))
  }
  pred.OR1 = expit(drop(Xi %*% OR1.coef[-1]) + OR1.coef[1])
  pred.OR0 = expit(drop(Xi %*% OR0.coef[-1]) + OR0.coef[1])
  pred.PST = expit(drop(Xi %*% PST.coef[-1]) + PST.coef[1])
  pred.R1 = expit(drop(Xi %*% PSR1.coef[-1]) + PSR1.coef[1])
  pred.R0 = expit(drop(Xi %*% PSR0.coef[-1]) + PSR0.coef[1])
  
  return( mean(pred.OR1-pred.OR0+ 
    (Yi-pred.OR1)*Ti*Ri/(pred.PST*pred.R1) - 
            (Yi-pred.OR0)*(1-Ti)*Ri/((1-pred.PST)*pred.R0)))
}


# The OR estimate for SSL
OR.SSL = function(Yi,Ti,Xi, Ri, Si, reg = TRUE)
{
  if(reg)
  {
    Y1.coef = adalasso(x = cbind(Xi[Ti*Ri==1,], Si[Ti*Ri==1,]),
                       y=Yi[Ti*Ri==1], 
                       family = "binomial")
    Yi[Ti*(1-Ri)==1] = expit(drop(cbind(Xi[Ti*(1-Ri)==1,],Si[Ti*(1-Ri)==1,]) 
                            %*% Y1.coef[-1]) + Y1.coef[1])
    Y0.coef = adalasso(x=cbind(Xi[(1-Ti)*Ri==1,], Si[(1-Ti)*Ri==1,]),
                       y=Yi[(1-Ti)*Ri==1], 
                       family = "binomial")
    Yi[(1-Ti)*(1-Ri)==1] = expit(drop(cbind(Xi[(1-Ti)*(1-Ri)==1,],Si[(1-Ti)*(1-Ri)==1,]) 
                                  %*% Y0.coef[-1]) + Y0.coef[1])
    Xtmp = rbind(Xi[Ti*Ri==1,],Xi[Ti*(1-Ri)==1,], Xi[Ti*(1-Ri)==1,])
    Ytmp = c(Yi[Ti*Ri==1], rep(1:0, each = sum(Ti*(1-Ri))))
    wgt = c(rep(1,sum(Ti*Ri)), Yi[Ti*(1-Ri)==1], 1-Yi[Ti*(1-Ri)==1])
    OR1.coef = adalasso(x=Xtmp, 
                        y=Ytmp, 
                        family = "binomial",
                        weights = wgt)
    Xtmp = rbind(Xi[(1-Ti)*Ri==1,],Xi[(1-Ti)*(1-Ri)==1,], Xi[(1-Ti)*(1-Ri)==1,])
    Ytmp = c(Yi[(1-Ti)*Ri==1], rep(1:0, each = sum((1-Ti)*(1-Ri))))
    wgt = c(rep(1,sum((1-Ti)*Ri)), Yi[(1-Ti)*(1-Ri)==1], 1-Yi[(1-Ti)*(1-Ri)==1])
    OR0.coef = adalasso(x=Xtmp, 
                        y=Ytmp, 
                        family = "binomial",
                        weights = wgt)
  }else
  {
    Y1.coef = coef(glm(Yi~Xi+Si, family = binomial, subset = (Ti*Ri==1)))
    Yi[Ti*(1-Ri)==1] = expit(drop(cbind(Xi[Ti*(1-Ri)==1,],Si[Ti*(1-Ri)==1,]) 
                                  %*% Y1.coef[-1]) + Y1.coef[1])
    Y0.coef = coef(glm(Yi~Xi+Si, family = binomial, subset = ((1-Ti)*Ri==1)))
    Yi[(1-Ti)*(1-Ri)==1] = expit(drop(cbind(Xi[(1-Ti)*(1-Ri)==1,],Si[(1-Ti)*(1-Ri)==1,]) 
                                      %*% Y0.coef[-1]) + Y0.coef[1])
    OR1.coef = coef(glm(Yi~Xi, family = binomial, subset = (Ti==1)))
    OR0.coef = coef(glm(Yi~Xi, family = binomial, subset = (Ti==0)))
  }
  pred.OR1 = expit(drop(Xi %*% OR1.coef[-1]) + OR1.coef[1])
  pred.OR0 = expit(drop(Xi %*% OR0.coef[-1]) + OR0.coef[1])
  
  return(mean(pred.OR1-pred.OR0))
}

# The IPW estimate for SSL (Similar to David's DiPS)
IPW.SSL = function(Yi,Ti,Xi, Ri, Si, reg = TRUE)
{
  if(reg)
  {
    PS.coef = adalasso(x=Xi, 
                        y=Ti, 
                        family = "binomial")
    pred.PS = expit(drop(Xi %*% PS.coef[-1]) + PS.coef[1])
    Y1.coef = adalasso(x = cbind(Xi, Si, 1/pred.PS)[Ti*Ri==1,],
                       y=Yi[Ti*Ri==1], 
                       family = "binomial")
    Y0.coef = adalasso(x = cbind(Xi, Si, 1/(1-pred.PS))[(1-Ti)*Ri==1,],
                       y=Yi[(1-Ti)*Ri==1], 
                       family = "binomial")
  }else
  {
    PS.coef = coef(glm(Ti~Xi, family = binomial))
    pred.PS = expit(drop(Xi %*% PS.coef[-1]) + PS.coef[1])
    U1 = 1/pred.PS
    Y1.coef = coef(glm(Yi~Xi+Si+U1, family = binomial, subset = (Ti*Ri==1)))
    U0 = 1/(1-pred.PS)
    Y0.coef = coef(glm(Yi~Xi+Si+U0, family = binomial, subset = ((1-Ti)*Ri==1)))
  }
  Yi[Ti*(1-Ri)==1] = expit(drop( cbind(Xi, Si, 1/pred.PS)[Ti*(1-Ri)==1,] 
                                %*% Y1.coef[-1]) + Y1.coef[1])
  Yi[(1-Ti)*(1-Ri)==1] = expit(drop(cbind(Xi, Si, 1/(1-pred.PS))[(1-Ti)*(1-Ri)==1,]
                                    %*% Y0.coef[-1]) + Y0.coef[1])
  
  return(mean(Yi*Ti/pred.PS - Yi*(1-Ti)/(1-pred.PS)))
}

# The AIPW estimate for SSL
AIPW.SSL = function(Yi,Ti,Xi, Ri, Si, reg = TRUE)
{
  if(reg)
  {
    Y1.coef = adalasso(x = cbind(Xi[Ti*Ri==1,], Si[Ti*Ri==1,]),
                       y=Yi[Ti*Ri==1], 
                       family = "binomial")
    Yi[Ti*(1-Ri)==1] = expit(drop(cbind(Xi[Ti*(1-Ri)==1,],Si[Ti*(1-Ri)==1,]) 
                                  %*% Y1.coef[-1]) + Y1.coef[1])
    Y0.coef = adalasso(x=cbind(Xi[(1-Ti)*Ri==1,], Si[(1-Ti)*Ri==1,]),
                       y=Yi[(1-Ti)*Ri==1], 
                       family = "binomial")
    Yi[(1-Ti)*(1-Ri)==1] = expit(drop(cbind(Xi[(1-Ti)*(1-Ri)==1,],Si[(1-Ti)*(1-Ri)==1,]) 
                                      %*% Y0.coef[-1]) + Y0.coef[1])
    Xtmp = rbind(Xi[Ti*Ri==1,],Xi[Ti*(1-Ri)==1,], Xi[Ti*(1-Ri)==1,])
    Ytmp = c(Yi[Ti*Ri==1], rep(1:0, each = sum(Ti*(1-Ri))))
    wgt = c(rep(1,sum(Ti*Ri)), Yi[Ti*(1-Ri)==1], 1-Yi[Ti*(1-Ri)==1])
    OR1.coef = adalasso(x=Xtmp, 
                        y=Ytmp, 
                        family = "binomial",
                        weights = wgt)
    Xtmp = rbind(Xi[(1-Ti)*Ri==1,],Xi[(1-Ti)*(1-Ri)==1,], Xi[(1-Ti)*(1-Ri)==1,])
    Ytmp = c(Yi[(1-Ti)*Ri==1], rep(1:0, each = sum((1-Ti)*(1-Ri))))
    wgt = c(rep(1,sum((1-Ti)*Ri)), Yi[(1-Ti)*(1-Ri)==1], 1-Yi[(1-Ti)*(1-Ri)==1])
    OR0.coef = adalasso(x=Xtmp, 
                        y=Ytmp, 
                        family = "binomial",
                        weights = wgt)
    PS.coef = adalasso(x=Xi, 
                       y=Ti, 
                       family = "binomial")
    pred.PS = expit(drop(Xi %*% PS.coef[-1]) + PS.coef[1])
    Y1.coef = adalasso(x = cbind(Xi, Si, 1/pred.PS)[Ti*Ri==1,],
                       y=Yi[Ti*Ri==1], 
                       family = "binomial")
    Y0.coef = adalasso(x = cbind(Xi, Si, 1/(1-pred.PS))[(1-Ti)*Ri==1,],
                       y=Yi[(1-Ti)*Ri==1], 
                       family = "binomial")
  }else
  {
    Y1.coef = coef(glm(Yi~Xi+Si, family = binomial, subset = (Ti*Ri==1)))
    Yi[Ti*(1-Ri)==1] = expit(drop(cbind(Xi[Ti*(1-Ri)==1,],Si[Ti*(1-Ri)==1,]) 
                                  %*% Y1.coef[-1]) + Y1.coef[1])
    Y0.coef = coef(glm(Yi~Xi+Si, family = binomial, subset = ((1-Ti)*Ri==1)))
    Yi[(1-Ti)*(1-Ri)==1] = expit(drop(cbind(Xi[(1-Ti)*(1-Ri)==1,],Si[(1-Ti)*(1-Ri)==1,]) 
                                      %*% Y0.coef[-1]) + Y0.coef[1])
    OR1.coef = coef(glm(Yi~Xi, family = binomial, subset = (Ti==1)))
    OR0.coef = coef(glm(Yi~Xi, family = binomial, subset = (Ti==0)))
    
    PS.coef = coef(glm(Ti~Xi, family = binomial))
    pred.PS = expit(drop(Xi %*% PS.coef[-1]) + PS.coef[1])
    U1 = 1/pred.PS
    Y1.coef = coef(glm(Yi~Xi+Si+U1, family = binomial, subset = (Ti*Ri==1)))
    U0 = 1/(1-pred.PS)
    Y0.coef = coef(glm(Yi~Xi+Si+U0, family = binomial, subset = ((1-Ti)*Ri==1)))
  }
  pred.OR1 = expit(drop(Xi %*% OR1.coef[-1]) + OR1.coef[1])
  pred.OR0 = expit(drop(Xi %*% OR0.coef[-1]) + OR0.coef[1])
  
  Yi[Ti*(1-Ri)==1] = expit(drop( cbind(Xi, Si, 1/pred.PS)[Ti*(1-Ri)==1,] 
                                 %*% Y1.coef[-1]) + Y1.coef[1])
  Yi[(1-Ti)*(1-Ri)==1] = expit(drop(cbind(Xi, Si, 1/(1-pred.PS))[(1-Ti)*(1-Ri)==1,]
                                    %*% Y0.coef[-1]) + Y0.coef[1])
  
  return(mean(pred.OR1-pred.OR0 + 
                (Yi-pred.OR1)*Ti/pred.PS-(Yi-pred.OR0)*(1-Ti)/(1-pred.PS)))
}