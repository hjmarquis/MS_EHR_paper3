
load("MS CLIME analysis/data-processed/common.rda")
MS_trt_new = MS_trt_new[!is.na(MS_trt_new$PatientNum),]

data.name.list = paste0("comb-data-",c("RN","DF"),"_MStrt")

datdir = paste('MS CLIME analysis/data-processed/', data.name.list[1],'.rds', sep='')
dat <- readRDS(datdir)


start.pos = match(paste0(dat$PatientNum,'_',dat$start_date),
                  paste0(MS_trt_new$PatientNum,'_',MS_trt_new$start_date))

MS_trt_RN = dat[,c("PatientNum", "medication_desc", "start_date","RELAPSE_DATE")]
MS_trt_RN$last_DMT_date = as.Date("2016-12-31")
MS_trt_RN$next_DMT_date = as.Date(NA)


for(i in 1:nrow(dat))
{
  for(j in start.pos[i]:nrow(MS_trt_new))
  {
    if(MS_trt_new$PatientNum[j] != MS_trt_RN$PatientNum[i])
    {
      break
    }
    if(MS_trt_new$medication_desc[j] != MS_trt_new$medication_desc[start.pos[i]])
    {
      MS_trt_RN$next_DMT_date[i] = MS_trt_new$start_date[j]
      break
    }
    MS_trt_RN$last_DMT_date[i] = MS_trt_new$stop_date[j]
  }
}



datdir = paste('MS CLIME analysis/data-processed/', data.name.list[2],'.rds', sep='')
dat <- readRDS(datdir)

start.pos = match(paste0(dat$PatientNum,'_',dat$start_date),
                  paste0(MS_trt_new$PatientNum,'_',MS_trt_new$start_date))

MS_trt_DF = dat[,c("PatientNum", "medication_desc", "start_date","RELAPSE_DATE")]
MS_trt_DF$last_DMT_date = as.Date("2016-12-31")
MS_trt_DF$next_DMT_date = as.Date(NA)

which(is.na(MS_trt_DF$last_DMT_date) & is.na(MS_trt_DF$next_DMT_date))
MS_trt_new[which(MS_trt_new$PatientNum == 21103),]

for(i in 1:nrow(dat))
{
  for(j in start.pos[i]:nrow(MS_trt_new))
  {
    if(MS_trt_new$PatientNum[j] != MS_trt_DF$PatientNum[i])
    {
      break
    }
    if(MS_trt_new$medication_desc[j] != MS_trt_new$medication_desc[start.pos[i]])
    {
      MS_trt_DF$next_DMT_date[i] = MS_trt_new$start_date[j]
      break
    }
    if(!is.na(MS_trt_new$stop_date[j]))
    {
      MS_trt_DF$last_DMT_date[i] = MS_trt_new$stop_date[j]
    }
  }
}


RTX.duration = as.numeric(pmax(MS_trt_RN$last_DMT_date, MS_trt_RN$next_DMT_date, na.rm = TRUE)-
                            MS_trt_RN$start_date,
                          units = "days")[MS_trt_RN$medication_desc=="Rituximab"]/365.2425
NTZ.duration = as.numeric(pmax(MS_trt_RN$last_DMT_date, MS_trt_RN$next_DMT_date, na.rm = TRUE)-
                            MS_trt_RN$start_date,
                          units = "days")[MS_trt_RN$medication_desc=="Natalizumab"]/365.2425

FGL.duration = as.numeric(pmax(MS_trt_DF$last_DMT_date, MS_trt_DF$next_DMT_date, na.rm = TRUE)-
                            MS_trt_DF$start_date,
                          units = "days")[MS_trt_DF$medication_desc=="Fingolimod"]/365.2425
DMF.duration = as.numeric(pmax(MS_trt_DF$last_DMT_date, MS_trt_DF$next_DMT_date, na.rm = TRUE)-
                            MS_trt_DF$start_date,
                          units = "days")[MS_trt_DF$medication_desc=="Dimethyl Fumarate"]/365.2425

library(ggplot2)
time.grid = sort(unique(c(0,pmin(c(RTX.duration,NTZ.duration),3))))

ggdf = data.frame( time = rep(time.grid,2), 
                   adherence = c(1-ecdf(RTX.duration)(time.grid),
                                 1-ecdf(NTZ.duration)(time.grid)),
                   trt = rep(c("RTX","NTZ"), each = length(time.grid)))

p = ggplot(ggdf, aes(x = time, y = adherence, group = trt,color = trt))
p = p + geom_step() + ylim(0,1)
p = p +  labs(title = "Adherence for RTX vs NTZ", color =  "Medication")
p = p + ylab("Adherence Rate") + xlab("Year(s) since Initiation of Medication")
p.RN <- p

time.grid = sort(unique(c(0,pmin(c(DMF.duration,FGL.duration),3))))

ggdf = data.frame( time = rep(time.grid,2), 
                   adherence = c(1-ecdf(DMF.duration)(time.grid),
                                 1-ecdf(FGL.duration)(time.grid)),
                   trt = rep(c("DMF","FGL"), each = length(time.grid)))

p = ggplot(ggdf, aes(x = time, y = adherence, group = trt,color = trt))
p = p + geom_step() + ylim(0,1)
p = p +  labs(title = "Adherence for DMF vs FGL", color =  "Medication")
p = p + ylab("Adherence Rate") + xlab("Year(s) since Initiation of Medication")
p.DF <- p

library(gridExtra)
library(grid)

# plist[[4]] = NULL

tiff("MS CLIME analysis/figure/paper3-adherence.tiff", width = 3500, height = 1500,
     res = 300)
grid.arrange(arrangeGrob(p.RN, p.DF, 
                         widths = c(1,1)
)
)
dev.off()

