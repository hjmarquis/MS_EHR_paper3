
dat.name.list = paste0("comb-data-",c("RN","DF"),"_MStrt")

est = rep(0,3)
boot.est = matrix(0,3,1000)

dat.name = dat.name.list[2]
load(paste("MS CLIME analysis/result/MS_result_new_",dat.name,"p10_yr.rda",sep=''))
load(paste0("MS CLIME analysis/result/MS_result_new_",dat.name,"p10_yr_ahaz_boot.rda"))
exp(-fit.time.adj$dr$est*2)

est[3] = fit.time.adj$dr$est
boot.est[3,] = boot.time.adj

quantile(exp(-boot.time.adj*2), c(0.025, 0.975))
mean(sign(fit.time.adj$dr$est)!=sign(boot.time.adj))*2

load(paste("MS CLIME analysis/result/MS_result_new_",dat.name,"p10_yr_dich2.rda",sep=''))
load(paste0("MS CLIME analysis/result/MS_result_new_",dat.name,"p10_yr_dich2_boot.rda"))

est[2] = fit2.time.adj$dr
boot.est[2,] = boot.time.adj

fit2.time.adj$dr
quantile(boot.time.adj, c(0.025, 0.975))
mean(sign(fit2.time.adj$dr)!=sign(boot.time.adj))*2

load(paste("MS CLIME analysis/result/MS_result_new_",dat.name,"p10_yr_dich1.rda",sep=''))
load(paste0("MS CLIME analysis/result/MS_result_new_",dat.name,"p10_yr_dich1_boot.rda"))

est[1] = fit1.time.adj$dr
boot.est[1,] = boot.time.adj
fit1.time.adj$dr
quantile(boot.time.adj, c(0.025, 0.975))
mean(sign(fit1.time.adj$dr)!=sign(boot.time.adj))*2

sig =apply(apply(boot.est,1,quantile, probs = c(0.025, 0.975)),
           2,diff)/(1.96*2)
est = apply(boot.est, 1, median)
Z.wald = abs(est/sig)
boot.std = (boot.est - est) / sig

apply(sign(boot.est)!= sign(est),1, mean)*2

apply(boot.std > Z.wald,1,mean) *2 
apply(outer(apply(boot.std, 2, max), Z.wald, '>'),2,mean) *2 
apply(outer(apply(abs(boot.std), 2, max), Z.wald, '>'),2,mean)


apply(boot.std < Z.wald,1,mean) *2 
apply(outer(apply(boot.std, 2, min), Z.wald, '<'),2,mean) *2 


yo = mvrnorm(10000, rep(0,3), var(t(boot.std)))

apply(outer(apply(abs(yo), 1, max), Z.wald, '>'),2,mean) 
