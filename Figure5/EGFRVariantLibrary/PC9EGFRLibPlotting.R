
### Plot Dynamics of PC9s with EGFR library C797S GD experiment

## Set up
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)

## Import data
df = read.csv("EGFRLibraryPC9Summary_031322.csv",header=T,stringsAsFactors=F)

## Plot data

GD.col = "#1AD31A"
lib.col = "#EE1D23"

df$population = factor(df$population,levels=c("lib","GD"))

# Plot stdev rather than sem so error bars are visible
df$sd = df$sem*sqrt(3)

# Condition 1
ggplot(df[df$condition==1,],aes(x=day,y=mean/1e6,color=population))+theme_bw()+
  geom_errorbar(aes(ymin=mean/1e6-sd/1e6,ymax=mean/1e6+sd/1e6),width=0.3,size=1.25)+
  geom_point(size=4)+
  geom_line(size=2)+
  xlab("Days")+ylab("Population Size (M)")+
  ggtitle("Without Gene Drive")+
  scale_y_continuous(breaks=c(0,1,2),limits=c(0,2.5))+
  scale_x_continuous(breaks=c(0,6,12,18))+
  scale_color_manual("",values=lib.col)+
  theme(
    plot.title = element_text(size=22,face="bold",hjust=0.5),
    axis.title = element_text(size=20,face="bold"),
    axis.text = element_text(size=18,face="bold",color="black"),
    legend.title = element_text(size=18,face="bold"),
    legend.text = element_text(size=16,face="bold")
  )+
  guides(color="none")

# ggsave("EGFRlibWithoutGeneDrive.png",width=5,height=5)

# Condition 2
ggplot(df[df$condition==2,],aes(x=day,y=mean/1e6,color=population))+theme_bw()+
  geom_errorbar(aes(ymin=mean/1e6-sd/1e6,ymax=mean/1e6+sd/1e6),width=0.3,size=1.25)+
  geom_point(size=4)+
  geom_line(size=2)+
  xlab("Days")+ylab("Population Size (M)")+
  ggtitle("With Gene Drive")+
  scale_y_continuous(breaks=c(0,1,2),limits=c(0,2.5))+
  scale_x_continuous(breaks=c(0,6,12,18,21),labels=c(0,6,12,18,33))+ # hack, rather than using ggbreak, just going to plot day 33 time point as "day 21" and then fix in illustrator
  scale_color_manual("",values=c(lib.col,GD.col))+
  theme(
    plot.title = element_text(size=22,face="bold",hjust=0.5),
    axis.title = element_text(size=20,face="bold"),
    axis.text = element_text(size=18,face="bold",color="black"),
    legend.title = element_text(size=18,face="bold"),
    legend.text = element_text(size=16,face="bold")
  )+
  guides(color="none")

# ggsave("EGFRlibWithGeneDrive.png",width=5,height=5)
