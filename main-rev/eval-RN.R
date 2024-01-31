
# E-value for ahaz in RR
#-------------------------------------------------------------

load( "MS CLIME analysis/data-analysis/comb-data-RN_MStrtp10_yr.rda")
load("MS CLIME analysis/result/MS_result_new_comb-data-RN_MStrtp10.rda")      

eval.ahaz.est = 1/0.903 # 1.107
eval.ahaz.ci = 1/0.944 # 1.059

# E-value for 2 year
#-------------------------------------------------------------

mu1 = 1/(1+exp(-fit2$or1[1]-drop(Z %*% fit2$or1[-1])))
mu0 = 1/(1+exp(-fit2$or0[1]-drop(Z %*% fit2$or0[-1])))
ps = 1/(1+exp(-fit2$ps[1]-drop(Z %*% fit2$ps[-1])))

nbin = 2
cut.bin = (1:(nbin-1))/nbin

mu1.bin = cut(mu1, c(0,quantile(mu1,cut.bin),1))
mu0.bin = cut(mu0, c(0,quantile(mu0,cut.bin),1))
ps.bin = cut(ps, c(0,quantile(ps,cut.bin),1))


table(mu1.bin,mu0.bin,ps.bin)

eval.tab = data.frame(mu1 = rep(levels(mu1.bin), each = nbin^2),
                      mu0 = rep(levels(mu0.bin), nbin, each = nbin),
                      ps = rep(levels(ps.bin), nbin^2),
                      n = NA, f = NA, p1 = NA, p0 = NA)
ps = pmin(0.9, pmax(0.1, ps))
ipw1 = D/ps
ipw0 = (1-D)/(1-ps)
ipw1 = ipw1/mean(ipw1)
ipw0 = ipw0/mean(ipw0)

for(i.mu1 in 1:nbin)
{
  for(i.mu0 in 1:nbin)
  {
    for(i.ps in 1:nbin)
    {
      tab.pos = nbin^2*(i.mu1-1) + nbin*(i.mu0-1) + i.ps
      bin.pos = which((mu1.bin==levels(mu1.bin)[i.mu1]) & 
                        (mu0.bin==levels(mu0.bin)[i.mu0]) &
                        (ps.bin==levels(ps.bin)[i.ps]))
      eval.tab$n[tab.pos] = length(bin.pos)
      eval.tab$f[tab.pos] = sum(D[bin.pos])
      eval.tab$p1[tab.pos] = mean(mu1[bin.pos] + D[bin.pos]*ipw1[bin.pos]*(Y2[bin.pos]-mu1[bin.pos]))
      eval.tab$p0[tab.pos] = mean(mu0[bin.pos] + (1-D)[bin.pos]*ipw0[bin.pos]*(Y2[bin.pos]-mu0[bin.pos]))
    }
  }
}

ate = ((eval.tab$p1 - eval.tab$p0) %*% eval.tab$n) / sum(eval.tab$n)

eval.rd = function(BF, eval.tab, lower = TRUE)
{
  if(lower)
  {
    return(((eval.tab$p1/BF-eval.tab$p0)*(eval.tab$f*BF+1-eval.tab$f))%*% eval.tab$n / sum(eval.tab$n))
  }else{
    return(((eval.tab$p1-eval.tab$p0/BF)*(eval.tab$f*BF+1-eval.tab$f))%*% eval.tab$n / sum(eval.tab$n))
  }
}

BF.grid = seq(1,3,0.05)

plot(BF.grid, sapply(BF.grid, eval.rd, eval.tab = eval.tab, lower = FALSE) , type ="l" )

eval.Y2.est = uniroot(function(x) eval.rd(x, eval.tab, lower = FALSE), 
        interval = c(1,3))$root # 2.26

B = 1000

BF.grid = seq(1,3,0.01)
rd.up = matrix(NA, B, length(BF.grid))
boot.tab = eval.tab
source("MS CLIME analysis/ate_dich.R")
for(b in 1:1000)
{
  bootfit = try(stop(),TRUE)
  
  while(inherits(bootfit,"try-error"))
  {
    bootid = sample(1:length(D),replace = T)
    bootfit = try(ate.dich.ada(Y2[bootid],
                             D[bootid],Z[bootid,],
                             fit2$or1.pen, fit2$or0.pen, fit2$ps.pen,
                             fit2$or1.lambda, fit2$or0.lambda, fit2$ps.lambda),TRUE)
  }
  mu1 = 1/(1+exp(-bootfit$or1[1]-drop(Z[bootid,] %*% bootfit$or1[-1])))
  mu0 = 1/(1+exp(-bootfit$or0[1]-drop(Z[bootid,] %*% bootfit$or0[-1])))
  ps = 1/(1+exp(-bootfit$ps[1]-drop(Z[bootid,] %*% bootfit$ps[-1])))
  
  ps = pmin(0.9, pmax(0.1, ps))
  ipw1 = D/ps
  ipw0 = (1-D)/(1-ps)
  ipw1 = ipw1/mean(ipw1)
  ipw0 = ipw0/mean(ipw0)
  for(i.mu1 in 1:nbin)
  {
    for(i.mu0 in 1:nbin)
    {
      for(i.ps in 1:nbin)
      {
        tab.pos = nbin^2*(i.mu1-1) + nbin*(i.mu0-1) + i.ps
        bin.pos = which((mu1.bin[bootid]==levels(mu1.bin)[i.mu1]) & 
                          (mu0.bin[bootid]==levels(mu0.bin)[i.mu0]) &
                          (ps.bin[bootid]==levels(ps.bin)[i.ps]))
        if(length(bin.pos)==0)
        {
          boot.tab$n[tab.pos]=boot.tab$f[tab.pos]=boot.tab$p1[tab.pos]=boot.tab$p0[tab.pos]=0
        }else{
          boot.tab$n[tab.pos] = length(bin.pos)
          boot.tab$f[tab.pos] = sum(D[bootid[bin.pos]])
          boot.tab$p1[tab.pos] = mean(mu1[bin.pos] + D[bootid[bin.pos]]*ipw1[bin.pos]*(Y2[bootid[bin.pos]]-mu1[bin.pos]))
          boot.tab$p0[tab.pos] = mean(mu0[bin.pos] + (1-D)[bootid[bin.pos]]*ipw0[bin.pos]*(Y2[bootid[bin.pos]]-mu0[bin.pos]))
        }
      }
    }
  }

  rd.up[b, ] = sapply(BF.grid, eval.rd, eval.tab = boot.tab, lower = FALSE)

}

rd.up.975 = apply(rd.up, 2, quantile, prob = 0.975)
plot(BF.grid, rd.up.975 , type ="l" ) # 1.31




# E-value for 1 year
#-------------------------------------------------------------

mu1 = 1/(1+exp(-fit1$or1[1]-drop(Z %*% fit1$or1[-1])))
mu0 = 1/(1+exp(-fit1$or0[1]-drop(Z %*% fit1$or0[-1])))
ps = 1/(1+exp(-fit1$ps[1]-drop(Z %*% fit1$ps[-1])))

nbin = 2
cut.bin = (1:(nbin-1))/nbin

# mu1.bin = cut(mu1, c(0,quantile(mu1,cut.bin),1))
mu0.bin = factor(mu0)
ps.bin = cut(ps, c(0,quantile(ps,cut.bin),1))


table(mu0.bin,ps.bin)

eval.tab = data.frame(mu0 = rep(levels(mu0.bin), each = nbin),
                      ps = rep(levels(ps.bin), nbin),
                      n = NA, f = NA, p1 = NA, p0 = NA)
ps = pmin(0.9, pmax(0.1, ps))
ipw1 = D/ps
ipw0 = (1-D)/(1-ps)
ipw1 = ipw1/mean(ipw1)
ipw0 = ipw0/mean(ipw0)


for(i.mu0 in 1:nbin)
{
  for(i.ps in 1:nbin)
  {
    tab.pos = nbin*(i.mu0-1) + i.ps
    bin.pos = which((mu0.bin==levels(mu0.bin)[i.mu0]) &
                      (ps.bin==levels(ps.bin)[i.ps]))
    eval.tab$n[tab.pos] = length(bin.pos)
    eval.tab$f[tab.pos] = sum(D[bin.pos])
    eval.tab$p1[tab.pos] = mean(mu1[bin.pos] + D[bin.pos]*ipw1[bin.pos]*(Y1[bin.pos]-mu1[bin.pos]))
    eval.tab$p0[tab.pos] = mean(mu0[bin.pos] + (1-D)[bin.pos]*ipw0[bin.pos]*(Y1[bin.pos]-mu0[bin.pos]))
  }
}

ate = ((eval.tab$p1 - eval.tab$p0) %*% eval.tab$n) / sum(eval.tab$n)

eval.rd = function(BF, eval.tab, lower = TRUE)
{
  if(lower)
  {
    return(((eval.tab$p1/BF-eval.tab$p0)*(eval.tab$f*BF+1-eval.tab$f))%*% eval.tab$n / sum(eval.tab$n))
  }else{
    return(((eval.tab$p1-eval.tab$p0/BF)*(eval.tab$f*BF+1-eval.tab$f))%*% eval.tab$n / sum(eval.tab$n))
  }
}

BF.grid = seq(1,5,0.05)

plot(BF.grid, sapply(BF.grid, eval.rd, eval.tab = eval.tab, lower = FALSE) , type ="l" )

eval.Y1.est = uniroot(function(x) eval.rd(x, eval.tab, lower = FALSE), 
                      interval = c(1,3))$root # 1.5

B = 1000

BF.grid = seq(1,2,0.01)
rd.up = matrix(NA, B, length(BF.grid))
boot.tab = eval.tab

for(b in 1:1000)
{
  bootfit = try(stop(),TRUE)
  
  while(inherits(bootfit,"try-error"))
  {
    bootid = sample(1:length(D),replace = T)
    bootfit = try(ate.dich.ada(Y1[bootid],
                               D[bootid],Z[bootid,],
                               fit1$or1.pen, fit1$or0.pen, fit1$ps.pen,
                               fit1$or1.lambda, fit1$or0.lambda, fit1$ps.lambda),TRUE)
  }
  mu1 = 1/(1+exp(-bootfit$or1[1]-drop(Z[bootid,] %*% bootfit$or1[-1])))
  mu0 = 1/(1+exp(-bootfit$or0[1]-drop(Z[bootid,] %*% bootfit$or0[-1])))
  ps = 1/(1+exp(-bootfit$ps[1]-drop(Z[bootid,] %*% bootfit$ps[-1])))
  
  ps = pmin(0.9, pmax(0.1, ps))
  ipw1 = D/ps
  ipw0 = (1-D)/(1-ps)
  ipw1 = ipw1/mean(ipw1)
  ipw0 = ipw0/mean(ipw0)

  for(i.mu0 in 1:nbin)
  {
    for(i.ps in 1:nbin)
    {
      tab.pos = nbin*(i.mu0-1) + i.ps
      bin.pos = which((mu0.bin[bootid]==levels(mu0.bin)[i.mu0]) &
                        (ps.bin[bootid]==levels(ps.bin)[i.ps]))
      if(length(bin.pos)==0)
      {
        boot.tab$n[tab.pos]=boot.tab$f[tab.pos]=boot.tab$p1[tab.pos]=boot.tab$p0[tab.pos]=0
      }else{
        boot.tab$n[tab.pos] = length(bin.pos)
        boot.tab$f[tab.pos] = sum(D[bin.pos])
        boot.tab$p1[tab.pos] = mean(mu1[bin.pos] + D[bootid[bin.pos]]*ipw1[bin.pos]*(Y1[bootid[bin.pos]]-mu1[bin.pos]))
        boot.tab$p0[tab.pos] = mean(mu0[bin.pos] + (1-D)[bootid[bin.pos]]*ipw0[bin.pos]*(Y1[bootid[bin.pos]]-mu0[bin.pos]))
      }
    }
  }
  
  rd.up[b, ] = sapply(BF.grid, eval.rd, eval.tab = boot.tab, lower = FALSE)
}


rd.up.975 = apply(rd.up, 2, quantile, prob = 0.975)
plot(BF.grid, rd.up.975 , type ="l" ) # 1.13
cbind(BF.grid, rd.up.975)
