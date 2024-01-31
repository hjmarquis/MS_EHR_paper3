
library(ggplot2)


pair.list = c("RN","DF")
pair.name = c("RTX vs NTZ", "DMF vs FGL")
arm1.name = c("RTX", "FGL")

heat.lim = c(-3.2,3.2)

for(k in 1:2)
{
  
  load(paste0("MS CLIME analysis/result/heat_map_data_", pair.list[k], ".rda"))
  
  ggdf = data.frame(
    Risk = rep(format(exp(as.numeric(colnames(BF.ahaz))*2),digits = 0, nsmall = 2), each = nrow(BF.ahaz)), 
    Propensity = rep(format(as.numeric(rownames(BF.ahaz)),digits = 0, nsmall = 2), ncol(BF.ahaz)),
    Confounding = log(c(BF.ahaz))
  )
  
  # heat.lim = c(-1,1)*ceiling(max(abs(ggdf$Confounding))*10)/10
  
  p <- ggplot(ggdf, aes(Propensity, Risk, fill = Confounding)) + geom_tile()
  p <- p + scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                                midpoint = 0, limit = heat.lim, space = "Lab", 
                                name="Log Confounding")
  p <- p + theme(legend.position="bottom")
  p <- p + labs(title=paste0(pair.name[k],": Time-to-relapse (Average Confounding ",
                             format(ave.BF.ahaz,digits=0,nsmall=2),")"))
  p <- p + ylab("Relative Risk for 2-year Relapse") + xlab(paste0("Propensity for ", arm1.name[k]))
  assign(paste0("p.",pair.list[k],".ahaz"),p)
  
  
  ggdf = data.frame(
    Risk = rep(format(as.numeric(colnames(BF.Y1)),digits = 0, nsmall = 2), each = nrow(BF.Y1)), 
    Propensity = rep(format(as.numeric(rownames(BF.Y1)),digits = 0, nsmall = 2), ncol(BF.Y1)),
    Confounding = log(c(BF.Y1))
  )
  # heat.lim = c(-1,1)*ceiling(max(abs(ggdf$Confounding))*10)/10
  
  p <- ggplot(ggdf, aes(Propensity, Risk, fill = Confounding)) + geom_tile()
  p <- p + scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                                midpoint = 0, limit = heat.lim, space = "Lab", 
                                name="Log Confounding")
  p <- p + theme(legend.position="bottom")
  p <- p + labs(title=paste0(pair.name[k],": 1 Year Relapse (Average Confounding ",
                             format(ave.BF.Y1,digits=0,nsmall=2),")"))
  p <- p + ylab("Odds Ratio for 1-year Relapse") + xlab(paste0("Propensity for ", arm1.name[k]))
  assign(paste0("p.",pair.list[k],".Y1"),p)
  
  
  ggdf = data.frame(
    Risk = rep(format(as.numeric(colnames(BF.Y2)),digits = 0, nsmall = 2), each = nrow(BF.Y2)), 
    Propensity = rep(format(as.numeric(rownames(BF.Y2)),digits = 0, nsmall = 2), ncol(BF.Y2)),
    Confounding = log(c(BF.Y2))
  )
  # heat.lim = c(-1,1)*ceiling(max(abs(ggdf$Confounding))*10)/10
  
  p <- ggplot(ggdf, aes(Propensity, Risk, fill = Confounding)) + geom_tile()
  p <- p + scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                                midpoint = 0, limit = heat.lim, space = "Lab", 
                                name="Log Confounding")
  p <- p + theme(legend.position="bottom")
  p <- p + labs(title=paste0(pair.name[k],": 2 Year Relapse (Average Confounding ",
                             format(ave.BF.Y2,digits=0,nsmall=2),")"))
  p <- p + ylab("Odds Ratio for 2-year Relapse") + xlab(paste0("Propensity for ", arm1.name[k]))
  assign(paste0("p.",pair.list[k],".Y2"),p)
}

p.RN.Y1
p.RN.Y2
p.RN.ahaz


p.DF.Y1
p.DF.Y2
p.DF.ahaz

library(gridExtra)
library(grid)

# plist[[4]] = NULL

tiff("MS CLIME analysis/figure/paper3-CLIMB_vs_EHR_heat_map.tiff", width = 5000, height = 3500,
     res = 300)
grid.arrange(arrangeGrob(p.RN.Y1, p.RN.Y2, p.RN.ahaz, 
                         p.DF.Y1, p.DF.Y2, p.DF.ahaz, 
                         layout_matrix = t(matrix(1:6,3,2)),
                         widths = c(1,1,1)
)
)
dev.off()
