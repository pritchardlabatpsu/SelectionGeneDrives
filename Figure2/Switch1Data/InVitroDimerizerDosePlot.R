
## Set up
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(grDevices)

## Import data
df = read.csv("BaF3DimerizerGrowthSummary_1101119.csv",header=T,stringsAsFactors=F) # from AnalyzeGrowthTracking1101119.R

df.sub = df[df$cell=="HME WT"&df$dose!=0.01&df$dose<=1000,]

## Plot data
ggplot(df.sub,aes(x=as.factor(dose),y=ave_rate))+theme_bw()+
  geom_errorbar(aes(ymin=ave_rate-sem,ymax=ave_rate+sem),
                width=0.25,size=1.5)+
  geom_bar(stat="identity",fill="#F99900",color="black")+
  xlab("Dimerizer Dose (nM)")+
  ylab("Growth Rate (/day)")+
  ggtitle("In Vitro Switch 1 vEGFRerl\nDimerizer-Dependent Growth")+
  scale_y_continuous(breaks=c(-0.5,0,.5))+
  theme(
    plot.title = element_text(size=22,face="bold",hjust=0.5),
    axis.title = element_text(size=20,face="bold"),
    axis.text = element_text(size=18,face="bold",color="black")
  )

# ggsave("InVitrooDoseResponseComplete.png",width=5,height=5)
