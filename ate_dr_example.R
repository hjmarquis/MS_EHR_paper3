source("ate_dr_code.R")

n = 200
p = 2

Z = matrix(rnorm(n*p),n,p)
trt = rbinom(n,1, 1/(1+exp(-X[,1]/2)))
Y = rbinom(n,1, 1/(1+exp(-X[,1]/2 + trt*X[,1])))

fit = ate.dr.adalasso(Y, trt, Z, 
                      or1.pen = rep(1,p),
                      or0.pen = rep(1,p),
                      ps.pen = rep(1,p),
                      or1.lambda = 0, 
                      or0.lambda = 0, 
                      ps.lambda = 0, 
                      lambda = 0)

