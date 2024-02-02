# yr.dich = 5
# strat = year.grp


require(glmnet)
ate.dich.ada = function(Y,trt,Z,
                              or1.pen, or0.pen, ps.pen, 
                              or1.lambda, or0.lambda, ps.lambda)
{
  n = length(Y)
  out = list()
  # OR
  if(missing(or1.pen))
  {
    tmp.ridge = cv.glmnet(Z[trt==1,], Y[trt==1],
                         family = "binomial", 
                         alpha = 0)
    
    out$or1.pen = or1.pen =  1/abs(coef(tmp.ridge, s = "lambda.min"))
  }
  
  if(missing(or1.lambda))
  {
    tmp.ada = cv.glmnet(Z[trt==1,], Y[trt==1],
                        family = "binomial", 
                       alpha = 1, penalty.factor = or1.pen)
    or1.coef = as.numeric(coef(tmp.ada, s= "lambda.min"))
    out$or1.lambda = or1.lambda = tmp.ada$lambda.min
  }else
  {
    tmp.ada = glmnet(Z[trt==1,], Y[trt==1],
                     family = "binomial", 
                    alpha = 1, penalty.factor = or1.pen)
    or1.coef = as.numeric(coef(tmp.ada, s= or1.lambda))
  }
  
  if(missing(or0.pen))
  {
    tmp.ridge = cv.glmnet(Z[trt==0,], Y[trt==0],
                          family = "binomial", 
                          alpha = 0)
    
    out$or0.pen = or0.pen =  1/abs(coef(tmp.ridge, s = "lambda.min"))
  }
  
  if(missing(or0.lambda))
  {
    tmp.ada = cv.glmnet(Z[trt==0,], Y[trt==0],
                        family = "binomial", 
                        alpha = 1, penalty.factor = or0.pen)
    or0.coef = as.numeric(coef(tmp.ada, s= "lambda.min"))
    out$or0.lambda = or0.lambda = tmp.ada$lambda.min
  }else
  {
    tmp.ada = glmnet(Z[trt==0,], Y[trt==0],
                     family = "binomial", 
                     alpha = 1, penalty.factor = or0.pen)
    or0.coef = as.numeric(coef(tmp.ada, s= or0.lambda))
  }

  # PS
  if(missing(ps.pen))
  {
    tmp.ridge = cv.glmnet(Z, trt, family = "binomial",
                         alpha = 0)
    out$ps.pen = ps.pen = 1/abs(coef(tmp.ridge, s = "lambda.min")[-1])
  }
  if(missing(ps.lambda))
  {
    tmp.ada = cv.glmnet(Z, trt, family = "binomial",
                       alpha = 1, penalty.factor = ps.pen)
    ps.coef = as.numeric(coef(tmp.ada, s= "lambda.min"))
    out$ps.lambda = ps.lambda = tmp.ada$lambda.min
  }else
  {
    tmp.ada = glmnet(Z, trt, family = "binomial",
                    alpha = 1, penalty.factor = ps.pen)
    ps.coef = as.numeric(coef(tmp.ada, s= ps.lambda))
  }
  
  or1 = 1/(1+exp(-or1.coef[1]-drop( Z %*% or1.coef[-1])))
  or0 = 1/(1+exp(-or0.coef[1]-drop( Z %*% or0.coef[-1])))
  
  ps = 1/(1+exp(-ps.coef[1]-drop( Z %*% ps.coef[-1])))
  ps = pmin(0.9, pmax(0.1, ps))
  ipw1 = trt/ps
  ipw0 = (1-trt)/(1-ps)
  ipw1 = ipw1/mean(ipw1)
  ipw0 = ipw0/mean(ipw0)
  
  ate.dich.dr = mean(or1 + ipw1*(Y-or1) - or0 - ipw0*(Y-or0))
  return(c(list(dr = ate.dich.dr, 
                or1 = or1.coef,or0 = or0.coef,
                ps = ps.coef),
           out))
}

ate.dich.ada.off = function(Y,trt,Z1, Z0, ZD, 
                            off1, off0, offD,
                        or1.pen, or0.pen, ps.pen, 
                        or1.lambda, or0.lambda, ps.lambda)
{
  n = length(Y)
  out = list()
  # OR
  if(ncol(Z1)>1)
  {
    if(missing(or1.pen))
    {
      tmp.ridge = cv.glmnet(Z1[trt==1,], Y[trt==1], offset = off1[trt==1], 
                            family = "binomial", 
                            alpha = 0)
      
      out$or1.pen = or1.pen =  1/abs(coef(tmp.ridge, s = "lambda.min"))
    }
    
    if(missing(or1.lambda))
    {
      tmp.ada = cv.glmnet(Z1[trt==1,], Y[trt==1], offset = off1[trt==1], 
                          family = "binomial", 
                          alpha = 1, penalty.factor = or1.pen)
      or1.coef = as.numeric(coef(tmp.ada, s= "lambda.min"))
      out$or1.lambda = or1.lambda = tmp.ada$lambda.min
    }else
    {
      tmp.ada = glmnet(Z1[trt==1,], Y[trt==1], offset = off1[trt==1], 
                       family = "binomial", 
                       alpha = 1, penalty.factor = or1.pen)
      or1.coef = as.numeric(coef(tmp.ada, s= or1.lambda))
    }
    or1 = 1/(1+exp(-or1.coef[1]-drop( Z1 %*% or1.coef[-1])-off1))
  }else if(ncol(Z1)==1){
    or1.coef = as.numeric(coef(glm(Y~Z1, family = binomial,
                                   offset = off1, 
                                   subset = trt == 1)))
    or1 = 1/(1+exp(-or1.coef[1]- drop(Z1* or1.coef[-1])-off1))
    out$or1.pen = out$or1.lambda = NA
  }else 
  {
    or1.coef = as.numeric(coef(glm(Y~1, family = binomial,
                                   offset = off1, 
                                   subset = trt == 1)))
    or1 = 1/(1+exp(-or1.coef-off1))
    out$or1.pen = out$or1.lambda = NA
  }
  
  if(ncol(Z0)>1)
  {
    if(missing(or0.pen))
    {
      tmp.ridge = cv.glmnet(Z0[trt==0,], Y[trt==0], offset = off0[trt==0], 
                            family = "binomial", 
                            alpha = 0)
      
      out$or0.pen = or0.pen =  1/abs(coef(tmp.ridge, s = "lambda.min"))
    }
    
    if(missing(or0.lambda))
    {
      tmp.ada = cv.glmnet(Z0[trt==0,], Y[trt==0], offset = off0[trt==0], 
                          family = "binomial", 
                          alpha = 1, penalty.factor = or0.pen)
      or0.coef = as.numeric(coef(tmp.ada, s= "lambda.min"))
      out$or0.lambda = or0.lambda = tmp.ada$lambda.min
    }else
    {
      tmp.ada = glmnet(Z0[trt==0,], Y[trt==0], offset = off0[trt==0], 
                       family = "binomial", 
                       alpha = 1, penalty.factor = or0.pen)
      or0.coef = as.numeric(coef(tmp.ada, s= or0.lambda))
    }
    or0 = 1/(1+exp(-or0.coef[1]-drop( Z0 %*% or0.coef[-1])-off0))
  }else if(ncol(Z0)==1){
    or0.coef = as.numeric(coef(glm(Y~Z0, family = binomial,
                                   offset = off0, 
                                   subset = trt == 0)))
    or0 = 1/(1+exp(-or0.coef[1]- drop(Z0* or0.coef[-1])-off0))
    out$or0.pen = out$or0.lambda = NA
  }else 
  {
    or0.coef = as.numeric(coef(glm(Y~1, family = binomial,
                                   offset = off0, 
                                   subset = trt == 0)))
    or0 = 1/(1+exp(-or0.coef-off0))
    out$or0.pen = out$or0.lambda = NA
  }
  
  # PS
  if(ncol(ZD) > 1)
  {
    if(missing(ps.pen))
    {
      tmp.ridge = cv.glmnet(ZD, trt,  offset = offD, 
                            family = "binomial",
                            alpha = 0)
      out$ps.pen = ps.pen = 1/abs(coef(tmp.ridge, s = "lambda.min")[-1])
    }
    if(missing(ps.lambda))
    {
      tmp.ada = cv.glmnet(ZD, trt,  offset = offD, 
                          family = "binomial",
                          alpha = 1, penalty.factor = ps.pen)
      ps.coef = as.numeric(coef(tmp.ada, s= "lambda.min"))
      out$ps.lambda = ps.lambda = tmp.ada$lambda.min
    }else
    {
      tmp.ada = glmnet(ZD, trt,  offset = offD, 
                       family = "binomial",
                       alpha = 1, penalty.factor = ps.pen)
      ps.coef = as.numeric(coef(tmp.ada, s= ps.lambda))
    }
    ps = 1/(1+exp(-ps.coef[1]-drop( ZD %*% ps.coef[-1])-offD))
  }else if(ncol(ZD)==1)
  {
    ps.coef = as.numeric(coef(glm(trt~ZD, family = binomial,
                                   offset = offD)))
    ps = 1/(1+exp(-ps.coef[1]-drop( ZD * ps.coef[-1])-offD))
    out$ps.pen = out$ps.lambda = NA
  }else
  {
    ps.coef =  as.numeric(coef(glm(trt~1, family = binomial,
                                    offset = offD)))
    ps = 1/(1+exp(-ps.coef-offD))
    out$ps.pen = out$ps.lambda = NA
  }
    
    
    
  ps = pmin(0.9, pmax(0.1, ps))
  ipw1 = trt/ps
  ipw0 = (1-trt)/(1-ps)
  ipw1 = ipw1/mean(ipw1)
  ipw0 = ipw0/mean(ipw0)
  
  ate.dich.dr = mean(or1 + ipw1*(Y-or1) - or0 - ipw0*(Y-or0))
  return(c(list(dr = ate.dich.dr, 
                or1 = or1.coef,or0 = or0.coef,
                ps = ps.coef),
           out))
}