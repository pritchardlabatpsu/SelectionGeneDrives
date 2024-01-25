
## Analyze tumor volumes for overall survival experiment

## Set up

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(NatParksPalettes)

## Analyze D0 tumor volumes
# Draw 3 mice to sack from 10% and 20% GD arms to analyze tumors

# Import D0 measurements
# df = read.table("PC9OS_TumorVolumes_092823.csv",header=T,sep=",",stringsAsFactors=F)
# 
# # Calculate volumes
# df$volume1 = pi*df$length1*df$width1*df$height1/6
# df$volume2 = pi*df$length2*df$width2*df$height2/6
# df$volume2 = ifelse(is.na(df$volume2),0,df$volume2)
# df$volume = rowSums(df[,c("volume1","volume2")])
# 
# # Draw 3 mice from each group with roughly same mean and standard deviation as population
# 
# set.seed(10)
# cndtns = c(10,20)
# to.sack = vector(length=length(cndtns))
# names(to.sack) = cndtns
# for (condition in cndtns) {
#   
#   mean = mean(df$volume[df$condition==condition])
#   sd = sd(df$volume[df$condition==condition])
#   
#   # Randomly subsample
#   subsamples = as.data.frame(matrix(nrow=10000,ncol=3))
#   colnames(subsamples) = c("mice","mean","sd")
#   for (i in 1:nrow(subsamples)) {
#     
#     # Randomly select 3 mice
#     sub.mice = sample(df$mouse[df$condition==condition],3)
#     subsamples$mice[i] = paste(sub.mice,collapse=",")
#     
#     # Calculate mean and sd of subset
#     subsamples$mean[i] = mean(df$volume[df$mouse%in%sub.mice])
#     subsamples$sd[i] = sd(df$volume[df$mouse%in%sub.mice])
#     
#   }
#   
#   subsamples$mean.diff = abs(subsamples$mean-mean)
#   subsamples$sd.diff = abs(subsamples$sd-sd)
#   
#   subsamples$tot.diff = subsamples$mean.diff+subsamples$sd.diff
#   
#   if (condition==10) {subsamples = subsamples[grepl("22",subsamples$mice),]} # for condition 10GD, include mouse 22 (outlier) to be sacked
# 
#   to.sack[as.character(condition)] = subsamples$mice[which.min(subsamples$tot.diff)] # record mice to sack
#   
# }
# 
# # Plot populations before and after sacking
# plot(df$condition,df$volume)
# 
# to.sack.vec = as.numeric(unlist(strsplit(to.sack,",")))
# plot(df$condition[!df$mouse%in%to.sack.vec],df$volume[!df$mouse%in%to.sack.vec])
# 
# # Summarize
# 
# sub.df = df[!df$mouse%in%to.sack.vec,]
# 
# smmry.df = sub.df %>%
#   group_by(condition) %>%
#   summarize(vol.mean = mean(volume,na.rm=T),
#             vol.sd = sd(volume,na.rm=T),
#             vol.sem = sd(volume,na.rm=T)/length(volume))
# 
# ggplot(smmry.df,aes(x=condition,y=vol.mean))+theme_bw()+
#   geom_errorbar(aes(ymin=vol.mean-vol.sd,ymax=vol.mean+vol.sd))+
#   geom_point()

## Analyze D22 tumor volumes
# Draw 3 mice to sack from gene drive arm to analyze tumors

# Import D22 measurements
df = read.table("PC9OS_TumorVolumes_102023.csv",header=T,sep=",",stringsAsFactors=F)

# Calculate volumes
df$volume1 = pi*df$length1*df$width1*df$height1/6
df$volume2 = pi*df$length2*df$width2*df$height2/6
df$volume2 = ifelse(is.na(df$volume2),0,df$volume2)
df$volume = rowSums(df[,c("volume1","volume2")])

# Draw 3 mice from gene drive group with roughly same mean and standard deviation as population on D22

set.seed(10)
cndtns = c(20)
to.sack = vector(length=length(cndtns))
names(to.sack) = cndtns
for (condition in cndtns) {

  mean = mean(df$volume[df$condition==condition])
  sd = sd(df$volume[df$condition==condition])

  # Randomly subsample
  subsamples = as.data.frame(matrix(nrow=10000,ncol=3))
  colnames(subsamples) = c("mice","mean","sd")
  for (i in 1:nrow(subsamples)) {

    # Randomly select 3 mice
    sub.mice = sample(df$mouse[df$condition==condition],3)
    subsamples$mice[i] = paste(sub.mice,collapse=",")

    # Calculate mean and sd of subset
    subsamples$mean[i] = mean(df$volume[df$mouse%in%sub.mice])
    subsamples$sd[i] = sd(df$volume[df$mouse%in%sub.mice])

  }

  subsamples$mean.diff = abs(subsamples$mean-mean)
  subsamples$sd.diff = abs(subsamples$sd-sd)

  subsamples$tot.diff = subsamples$mean.diff+subsamples$sd.diff

  if (condition==10) {subsamples = subsamples[grepl("22",subsamples$mice),]} # for condition 10GD, include mouse 22 (outlier) to be sacked

  to.sack[as.character(condition)] = subsamples$mice[which.min(subsamples$tot.diff)] # record mice to sack

}

# Plot populations before and after sacking
# plot(df$condition,df$volume)

to.sack.vec = as.numeric(unlist(strsplit(to.sack,",")))
# plot(df$condition[!df$mouse%in%to.sack.vec],df$volume[!df$mouse%in%to.sack.vec])

# Summarize

sub.df = df[!df$mouse%in%to.sack.vec,]

smmry.df = sub.df %>%
  group_by(condition) %>%
  summarize(vol.mean = mean(volume,na.rm=T),
            vol.sd = sd(volume,na.rm=T),
            vol.sem = sd(volume,na.rm=T)/length(volume))

# ggplot(smmry.df,aes(x=condition,y=vol.mean))+theme_bw()+
#   geom_errorbar(aes(ymin=vol.mean-vol.sd,ymax=vol.mean+vol.sd))+
#   geom_point()

## Plot tumor volumes over time

# Import and parse data

csv.files = list.files(pattern="PC9OS_TumorVolumes_.*csv$")
cmpld.df = data.frame()
start.date = as.Date("2023-09-28")

for (i in 1:length(csv.files)) {
  
  file.i = csv.files[i]
  date.i = as.Date(strsplit(file.i,"_|\\.")[[1]][3],format = "%m%d%y")
  
  df.i = read.table(file.i,header=T,stringsAsFactors=F,sep=",")
  
  df.i = df.i %>%
    mutate(day = as.numeric(date.i-start.date),.before=mouse)
  
  cmpld.df = rbind(cmpld.df,df.i)
  
}

# Filter for mice that weren't sacked
mice = cmpld.df$mouse[cmpld.df$day==4]
cmpld.df = cmpld.df %>%
  filter(mouse%in%mice)

cmpld.df$mouse = as.factor(cmpld.df$mouse)
cmpld.df$condition[cmpld.df$condition==20] = 10

# Intermediate calculations
cmpld.df$volume1 = pi*cmpld.df$length1*cmpld.df$width1*cmpld.df$height1/6
cmpld.df$length2 = as.numeric(cmpld.df$length2)
cmpld.df$volume2 = pi*cmpld.df$length2*cmpld.df$width2*cmpld.df$height2/6
cmpld.df$volume2 = ifelse(is.na(cmpld.df$volume2),0,cmpld.df$volume2)
cmpld.df$volume = rowSums(cmpld.df[,c("volume1","volume2")])

# Find relative change in volume
cmpld.df$rel.volume = NA
for (mouse.i in mice) {
  
  init.vol = cmpld.df$volume[cmpld.df$day==0&cmpld.df$mouse==mouse.i]
  
  cmpld.df$rel.volume[cmpld.df$mouse==mouse.i] = cmpld.df$volume[cmpld.df$mouse==mouse.i]/init.vol
  
}

# Plot individual trajectories
to.sack = c(4,33,40,47) # note: also including mouse 4 here, because it was inadvertently lost by ARP on 11/29
cndtns = c(0,10)
clrs = rev(natparks.pals("Acadia"))[c(1,9)]
for (i in 1:2) {
  
  cndtn = cndtns[i]
  plot.cndtn = cmpld.df$condition==cndtn&!cmpld.df$mouse%in%to.sack # remove S2 sacked mice
  
  # # Plot absolute volumes
  # g=ggplot(cmpld.df[plot.cndtn,],aes(x=day,y=volume,color=mouse.flank))+theme_bw()+
  #   geom_point()+
  #   geom_line()+
  #   ggtitle(paste("Condition:",cndtn))+
  #   theme(plot.title = element_text(hjust=0.5))
  
  # Plot relative volumes
  n.mice = length(unique(cmpld.df$mouse[plot.cndtn]))
  clr = clrs[i]
  
  g=ggplot(cmpld.df[plot.cndtn,],aes(x=day,y=rel.volume,color=mouse))+theme_bw()+
    geom_point(size=1.5,alpha=0.8)+
    geom_line(size=1.2,alpha=0.8)+
    # geom_vline(xintercept=57)+
    ggtitle(paste0(cndtn,"% Gene Drive"))+
    xlab("Days")+ylab("Relative Tumor Volume")+
    scale_color_manual(values=rep(clr,n.mice))+
    scale_y_continuous(limits=c(0,1.6))+
    scale_x_continuous(breaks=c(0,20,40,60,80),limits=c(0,85))+
    theme(plot.title = element_text(size=20,hjust=0.5,face="bold",color=clr),
          axis.title = element_text(size=18,face="bold"),
          axis.text = element_text(size=16,color="black",face="bold"))+
    guides(color="none")
  print(g)
  # ggsave(paste0(cndtn,"GeneDriveRelTumorVolume.png"),width=5,height=5)

}


# Summarize and plot data

smmry.df = cmpld.df %>%
  group_by(condition,day) %>%
  summarize(vol.mean = mean(volume),
            vol.sem = sd(volume)/sqrt(length(volume)),
            rel.vol.mean = mean(rel.volume),
            rel.vol.sem = sd(rel.volume)/sqrt(length(rel.volume)))


## Plot survival curve

# Import and parse data
sack.df = read.table("PC9OS_TerminalEndpoints.csv",sep=",",header=T,stringsAsFactors=F)

sack.df$sack.date = as.Date(sack.df$sack.date,format = "%m/%d/%y")
sack.df$day = sack.df$sack.date-start.date

sack.df = sack.df[order(sack.df$day),]

cndtns = c(0,10)
n.mice.0 = c(sum(sack.df$condition==0),sum(sack.df$condition==10))
n.mice = n.mice.0

surv.plt = data.frame(day = c(0,0),
                      condition = cndtns,
                      n.mice = n.mice)

for (i in 1:nrow(sack.df)) {
  
  if (!is.na(sack.df$day[i])) {
    
    sack.day.i = sack.df$day[i]
    cndtn.i = sack.df$condition[i]
    n.mice.i = n.mice[cndtns==cndtn.i]
    
    n.mice[cndtns==cndtn.i] = n.mice.i-1
    
    sub.df = data.frame(day = rep(sack.day.i,2),
                        condition = rep(cndtn.i,2),
                        n.mice = c(n.mice.i,n.mice.i-1))
    
    surv.plt = rbind(surv.plt,
                     sub.df)
    
  }
  
}

# Add current mice counts
end.date = Sys.Date()-start.date
end.date = as.Date("2023-12-22")-start.date

sub.df = data.frame(day = rep(end.date,2),
                    condition = cndtns,
                    n.mice = n.mice)
surv.plt = rbind(surv.plt,
                 sub.df)

# Calculate percent survival
for (cndtn.i in cndtns) {
  surv.plt$perc.surv[surv.plt$condition==cndtn.i] = surv.plt$n.mice[surv.plt$condition==cndtn.i]/n.mice.0[cndtns==cndtn.i]
}

surv.plt$condition = as.factor(surv.plt$condition)
clrs = rev(natparks.pals("Acadia"))[c(1,9)]

ggplot(surv.plt,aes(x=day,y=perc.surv*100,color=condition))+theme_bw()+
  geom_line(size=2)+
  scale_color_manual(values=clrs)+
  xlab("Days")+ylab("Percent Survival")+
  scale_y_continuous(limits=c(0,100))+
  scale_x_continuous(breaks=c(0,20,40,60,80))+
  ggtitle("Survival Curve")+
  theme(plot.title = element_text(size=20,hjust=0.5,face="bold"),
        axis.title = element_text(size=18,face="bold"),
        axis.text = element_text(size=16,color="black",face="bold"))+
  guides(color="none")
# ggsave("SurvivalCurve.png",width=5,height=4.5)


## Log-rank test

library(survival)

sack.df$day[is.na(sack.df$day)] = end.date

survdiff(Surv(day)~condition,data=sack.df)

## Plot absolute tumor volumes at day 0

ggplot(cmpld.df[cmpld.df$day==0,],aes(x=as.factor(condition),y=volume,color=as.factor(condition)))+theme_bw()+
  geom_boxplot(size=1.5)+
  scale_color_manual(values=clrs)+
  xlab("Group")+
  ylab("D0 Tumor Volume (mm3)")+
  ggtitle("Initial Tumor Volume")+
  scale_x_discrete(labels=c("0% Gene Drive","10% Gene Drive"))+
  scale_y_continuous(limits=c(50,450))+
  theme(plot.title = element_text(size=20,hjust=0.5,face="bold"),
        axis.title = element_text(size=18,face="bold"),
        axis.text = element_text(size=16,color="black",face="bold"),
        axis.title.x = element_blank())+
  guides(color="none")
# ggsave("D0_TumorVolumeBoxplot.png",width=5,height=4.5)

## Alternatively, plot just gene drive tumor volumes at day 0

library(ggbeeswarm)

cmpld.df$m.35 = cmpld.df$mouse==35

ggplot(cmpld.df[cmpld.df$day==0&cmpld.df$condition==10,],aes(x=as.factor(condition),y=volume))+theme_bw()+
  geom_boxplot(width=0.2)+
  geom_dotplot(aes(fill=m.35),dotsize=1.5,alpha=0.85,binaxis = "y",stackdir = "center")+
  scale_fill_manual(values=c("#94A5B3","#D65C83"))+
  xlab("Group")+
  ylab("Day 0 Tumor Volume (mm3)")+
  ggtitle("Initial\nTumor Volume")+
  scale_x_discrete(labels=c("0% Gene Drive","10% Gene Drive"))+
  scale_y_continuous(limits=c(50,450))+
  theme(plot.title = element_text(size=20,hjust=0.5,face="bold"),
        axis.title = element_text(size=18,face="bold"),
        axis.text = element_text(size=16,color="black",face="bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  guides(fill="none")
# ggsave("D0_TumorVolume.png",width=3,height=4.5)

## Plot gene drive relative tumor volumes at day 22

ggplot(cmpld.df[cmpld.df$day==22&cmpld.df$condition==10,],aes(x=as.factor(condition),y=rel.volume))+theme_bw()+
  geom_boxplot(width=0.2)+
  geom_dotplot(aes(fill=m.35),dotsize=1.5,alpha=0.85,binaxis = "y",stackdir = "center")+
  scale_fill_manual(values=c("#94A5B3","#D65C83"))+
  xlab("Group")+
  ylab("Day 22 Tumor Volume\nRelative to Day 0 Tumor Volume")+
  ggtitle("Rebound\nTumor Volume")+
  scale_x_discrete(labels=c("0% Gene Drive","10% Gene Drive"))+
  scale_y_continuous(limits=c(0.6,1.8))+
  theme(plot.title = element_text(size=20,hjust=0.5,face="bold"),
        axis.title = element_text(size=18,face="bold"),
        axis.text = element_text(size=16,color="black",face="bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  guides(fill="none")
# ggsave("D22_RelTumorVolume.png",width=3,height=4.5)

## Plot %GD at D0 and D22

GD.df = data.frame(
  day = c(0,0,0,22,22,22),
  perc.GFP = c(7.3,8.7,12.4,94.61,90.97,86.03)
)

ggplot(GD.df,aes(x=as.factor(day),y=perc.GFP))+theme_bw()+
  geom_point(size=5,alpha=0.5,fill="#1AD31A",color="black",shape=21,stroke=1.5)+
  xlab("Treatment Day")+
  ylab("Percent Gene Drive Cells")+
  ggtitle("Subpopulation Analysis")+
  scale_y_continuous(limits=c(0,100),breaks=c(0,20,40,60,80,100))+
  theme(plot.title = element_text(size=20,hjust=0.5,face="bold"),
        axis.title = element_text(size=18,face="bold"),
        axis.text = element_text(size=16,color="black",face="bold"))

# ggsave("SubpopulationAnalysis.png",width=4,height=4.5)

## Plot weight data

# Calculate change in mouse weight
cmpld.df$rel.weight = NA
for (mouse.i in mice) {
  
  init.weight = cmpld.df$weight[cmpld.df$day==0&cmpld.df$mouse==mouse.i]
  
  cmpld.df$rel.weight[cmpld.df$mouse==mouse.i] = cmpld.df$weight[cmpld.df$mouse==mouse.i]/init.weight
  
}

# Summarize

wght.smmry = cmpld.df %>%
  group_by(day,condition) %>%
  summarize(mean.rel.weight = mean(rel.weight,na.rm=T),
            sem.rel.weight = sd(rel.weight,na.rm=T)/sqrt(length(rel.weight))) %>%
  filter(!is.nan(mean.rel.weight))

ggplot(wght.smmry[wght.smmry$day<=60,],aes(x=day,y=mean.rel.weight,color=as.factor(condition)))+theme_bw()+
  geom_errorbar(aes(ymin=mean.rel.weight-sem.rel.weight,ymax=mean.rel.weight+sem.rel.weight),size=1.25)+
  geom_point(size=2.5)+
  geom_line(size=2)+
  scale_color_manual("Condition",labels=c("0% Gene Drive","10% Gene Drive"),values=clrs)+
  scale_y_continuous(limits=c(0.9,1.15),breaks=c(0.9,1,1.1))+
  ggtitle("Weight Tracking")+
  xlab("Day")+ylab("Weight Relative to Day 0")+
  theme(plot.title = element_text(size=20,hjust=0.5,face="bold"),
        axis.title = element_text(size=18,face="bold"),
        axis.text = element_text(size=16,color="black",face="bold"),
        legend.title = element_text(size=16,color="black",face="bold"),
        legend.text = element_text(size=14,color="black",face="bold"))
# ggsave("WeightTracking.png",width=6.5,height=4.5)
