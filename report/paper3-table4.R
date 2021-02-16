
# Table for RN

dat.name.list = paste(c("low-d","comb"), "-data-",
                      "RN", '_',
                      "MStrt",sep = "")

out.table = data.frame(outcome = rep(c("1Y", "2Y", "T2R"),each = 3),
                       method = rep(c("OR","IPW","DR"),3),
                       est.ld = "",  CI.ld = "", p.ld = "",
                       est.hd = "",  CI.hd = "", p.hd = "",
                       stringsAsFactors = F)

load( paste("MS CLIME analysis/result/MS_result_new_",dat.name.list[1],"p10.rda",sep=''))
boot = boot[,1:3]
out.table$est.ld[7] = format(exp(fit.alltil17$or * 2),digits=0, nsmall = 3)
out.table$est.ld[8] = format(exp(fit.alltil17$ipw * 2),digits=0, nsmall = 3)
out.table$est.ld[9] = format(exp(fit.alltil17$dr$est * 2),digits=0, nsmall = 3)

out.table$CI.ld[7:9] = paste("(",apply(format(apply(exp(2*boot),2,quantile,probs = c(.025,.975)),
                                                    digits=0, nsmall = 3),2,paste,collapse = ","),
                                   ")", sep='')
out.table$p.ld[7:9] = format(apply(sign(boot)== -sign(matrix(rep(c(fit.alltil17$or,
                                                                            fit.alltil17$ipw,
                                                                            fit.alltil17$dr$est),
                                                                          each = nrow(boot)),nrow(boot),3)),2,mean)*2,
                                      digits=0, nsmall = 3)

load( paste("MS CLIME analysis/result/MS_result_new_",dat.name.list[2],"p10.rda",sep=''))
boot = boot[,1:3]
out.table$est.hd[7] = format(exp(fit.alltil17$or * 2),digits=0, nsmall = 3)
out.table$est.hd[8] = format(exp(fit.alltil17$ipw * 2),digits=0, nsmall = 3)
out.table$est.hd[9] = format(exp(fit.alltil17$dr$est * 2),digits=0, nsmall = 3)

out.table$CI.hd[7:9] = paste("(",apply(format(apply(exp(2*boot),2,quantile,probs = c(.025,.975)),
                                              digits=0, nsmall = 3),2,paste,collapse = ","),
                             ")", sep='')
out.table$p.hd[7:9] = format(apply(sign(boot)== -sign(matrix(rep(c(fit.alltil17$or,
                                                                   fit.alltil17$ipw,
                                                                   fit.alltil17$dr$est),
                                                                 each = nrow(boot)),nrow(boot),3)),2,mean)*2,
                             digits=0, nsmall = 3)

write.csv(out.table,"MS CLIME analysis/result/table4RN.csv", row.names = F)



# Table for DF

dat.name.list = paste(c("low-d","comb"), "-data-",
                      "DF", '_',
                      "MStrt",sep = "")

out.table = data.frame(outcome = rep(c("1Y", "2Y", "T2R"),each = 3),
                       method = rep(c("OR","IPW","DR"),3),
                       est.ld = "",  CI.ld = "", p.ld = "",
                       est.hd = "",  CI.hd = "", p.hd = "",
                       stringsAsFactors = F)

load( paste("MS CLIME analysis/result/MS_result_new_",dat.name.list[1],"p10.rda",sep=''))
boot = boot[,1:3]
out.table$est.ld[7] = format(exp(-fit.alltil17$or * 2),digits=0, nsmall = 3)
out.table$est.ld[8] = format(exp(-fit.alltil17$ipw * 2),digits=0, nsmall = 3)
out.table$est.ld[9] = format(exp(-fit.alltil17$dr$est * 2),digits=0, nsmall = 3)

out.table$CI.ld[7:9] = paste("(",apply(format(apply(exp(-2*boot),2,quantile,probs = c(.025,.975)),
                                              digits=0, nsmall = 3),2,paste,collapse = ","),
                             ")", sep='')
out.table$p.ld[7:9] = format(apply(sign(boot)== -sign(matrix(rep(c(fit.alltil17$or,
                                                                   fit.alltil17$ipw,
                                                                   fit.alltil17$dr$est),
                                                                 each = nrow(boot)),nrow(boot),3)),2,mean)*2,
                             digits=0, nsmall = 3)

load( paste("MS CLIME analysis/result/MS_result_new_",dat.name.list[2],"p10.rda",sep=''))
boot = boot[,1:3]
out.table$est.hd[7] = format(exp(-fit.alltil17$or * 2),digits=0, nsmall = 3)
out.table$est.hd[8] = format(exp(-fit.alltil17$ipw * 2),digits=0, nsmall = 3)
out.table$est.hd[9] = format(exp(-fit.alltil17$dr$est * 2),digits=0, nsmall = 3)

out.table$CI.hd[7:9] = paste("(",apply(format(apply(exp(-2*boot),2,quantile,probs = c(.025,.975)),
                                              digits=0, nsmall = 3),2,paste,collapse = ","),
                             ")", sep='')
out.table$p.hd[7:9] = format(apply(sign(boot)== -sign(matrix(rep(c(fit.alltil17$or,
                                                                   fit.alltil17$ipw,
                                                                   fit.alltil17$dr$est),
                                                                 each = nrow(boot)),nrow(boot),3)),2,mean)*2,
                             digits=0, nsmall = 3)

write.csv(out.table,"MS CLIME analysis/result/table4DF.csv", row.names = F)
