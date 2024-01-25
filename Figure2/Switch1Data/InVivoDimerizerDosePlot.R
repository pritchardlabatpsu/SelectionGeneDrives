
## Set up
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)

## Import data (already summarized in DimerizerDoseCheck_20220121.xlsx)
df = read.csv("DimerizerDoseCheck_20220121Summarized.csv",header=T,stringsAsFactors=F)

## Plot data
ggplot(df,aes(x=as.factor(dose.mg.kg),y=tumor.volume.mean))+theme_bw()+
  geom_errorbar(aes(ymin=tumor.volume.mean-tumor.volume.sem,ymax=tumor.volume.mean+tumor.volume.sem),
                width=0.25,size=1.5)+
  geom_bar(stat="identity",fill="#F99900",color="black")+
  xlab("Dose (mg/kg)")+
  ylab("Tumor Volume (mm3)")+
  ggtitle("Dimerizer-Dependent Growth")+
  scale_y_continuous(breaks=c(0,100,200))+
  theme(
    plot.title = element_text(size=22,face="bold",hjust=0.5),
    axis.title = element_text(size=20,face="bold"),
    axis.text = element_text(size=18,face="bold",color="black")
  )

# ggsave("InVivoDoseResponse.png",width=5,height=5)
