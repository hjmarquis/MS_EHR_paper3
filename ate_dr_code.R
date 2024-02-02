
# Inputs
#-----------------------------------------------
# Must supply Y, trt, Z. 
# Y: vector of outcomes
# trt: vector of treatment assignments
# Z: matrix of confounding adjustments
# family: indicator for binary or continuous Y
# psmax and psmin: maximal and minimal ps allowed. 
#                  Obeservations with extreme ps will be dicarded
# or1.pen, or0pen, ps.pen, or1.lambda, or0.lambda, ps.lambda:
#   user specified parameters for adpative lasso

# Outputs
#-----------------------------------------------
#
# dr: doubly robust estimator for ATE
# or1, or0, ps: model coefficients for OR and PS models
# or1.pen, or0pen, ps.pen, or1.lambda, or0.lambda, ps.lambda: 
#   parameters used for adaptive lasso

require(glmnet)
ate.dr.adalasso = function(Y,trt,Z, 
                           family = ifelse(all(Y %in% 0:1), 
                                           "binomial",
                                           "gaussian"),
                           psmax = 1, psmin=0,
                              or1.pen, or0.pen, ps.pen, 
                              or1.lambda, or0.lambda, ps.lambda, 
                           ...)
{
  n = length(Y)
  out = list()
  # OR
  if(missing(or1.pen))
  {
    tmp.ridge = cv.glmnet(Z[trt==1,], Y[trt==1],
                         family = family, 
                         alpha = 0,...)
    
    out$or1.pen = or1.pen =  1/abs(coef(tmp.ridge, s = "lambda.min"))
  }
  
  if(missing(or1.lambda))
  {
    tmp.ada = cv.glmnet(Z[trt==1,], Y[trt==1],
                        family = family, 
                       alpha = 1, penalty.factor = or1.pen, 
                       ...)
    or1.coef = as.numeric(coef(tmp.ada, s= "lambda.min"))
    out$or1.lambda = or1.lambda = tmp.ada$lambda.min
  }else
  {
    tmp.ada = glmnet(Z[trt==1,], Y[trt==1],
                     family = family, 
                    alpha = 1, penalty.factor = or1.pen, 
                    ...)
    or1.coef = as.numeric(coef(tmp.ada, s= or1.lambda))
  }
  
  if(missing(or0.pen))
  {
    tmp.ridge = cv.glmnet(Z[trt==0,], Y[trt==0],
                          family = family, 
                          alpha = 0, 
                          ...)
    
    out$or0.pen = or0.pen =  1/abs(coef(tmp.ridge, s = "lambda.min"))
  }
  
  if(missing(or0.lambda))
  {
    tmp.ada = cv.glmnet(Z[trt==0,], Y[trt==0],
                        family = family, 
                        alpha = 1, penalty.factor = or0.pen, 
                        ...)
    or0.coef = as.numeric(coef(tmp.ada, s= "lambda.min"))
    out$or0.lambda = or0.lambda = tmp.ada$lambda.min
  }else
  {
    tmp.ada = glmnet(Z[trt==0,], Y[trt==0],
                     family = family, 
                     alpha = 1, penalty.factor = or0.pen, 
                     ...)
    or0.coef = as.numeric(coef(tmp.ada, s= or0.lambda))
  }

  # PS
  if(missing(ps.pen))
  {
    tmp.ridge = cv.glmnet(Z, trt, family = family,
                         alpha = 0, 
                         ...)
    out$ps.pen = ps.pen = 1/abs(coef(tmp.ridge, s = "lambda.min")[-1])
  }
  if(missing(ps.lambda))
  {
    tmp.ada = cv.glmnet(Z, trt, family = family,
                       alpha = 1, penalty.factor = ps.pen, 
                       ...)
    ps.coef = as.numeric(coef(tmp.ada, s= "lambda.min"))
    out$ps.lambda = ps.lambda = tmp.ada$lambda.min
  }else
  {
    tmp.ada = glmnet(Z, trt, family = family,
                    alpha = 1, penalty.factor = ps.pen, 
                    ...)
    ps.coef = as.numeric(coef(tmp.ada, s= ps.lambda))
  }
  
  if(family == "binomial")
  {
    link = function(x) 1/(1+exp(-x))
  }else{
    link = identity
  }
  or1 = link(or1.coef[1]+drop( Z %*% or1.coef[-1]))
  or0 = link(or0.coef[1]+drop( Z %*% or0.coef[-1]))
  
  ps = 1/(1+exp(-ps.coef[1]-drop( Z %*% ps.coef[-1])))
  ps = pmin(psmax, pmax(psmin, ps))
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
