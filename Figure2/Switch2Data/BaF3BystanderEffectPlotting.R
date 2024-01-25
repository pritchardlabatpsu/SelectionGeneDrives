
## Set up
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)

df = read.csv("BaF3CDBystander.csv",header=T)

df$killing = 1-df$mean.viability

m = (df$killing[df$percent.CD==100]-df$killing[df$percent.CD==0])/(100-0)
b = df$killing[df$percent.CD==0]

df$no.bystndr = m*df$percent.CD+b

ggplot(df,aes(x=percent.CD,y=killing))+theme_bw()+
  geom_ribbon(aes(ymin=no.bystndr,ymax=killing),fill="#BA50FF",alpha=0.4)+
  geom_line(size=2,color="#8F00FF",)+
  geom_segment(x=0,y=b,xend=100,yend=b+100*m,size=2,color="gray35")+
  geom_errorbar(aes(ymin=killing-sd.viability,ymax=killing+sd.viability),width=4,size=1.5,color="black")+
  geom_point(size=5,shape=21,fill="#8F00FF",color="black")+
  geom_point(data=df[df$percent.CD%in%c(0,100),],size=5,shape=21,fill="gray35",color="black")+
  scale_size_manual(values = c(8,5,5,5,8))+
  scale_x_continuous(breaks=c(0,50,100))+
  scale_y_continuous(breaks=c(0,0.5,1))+
  xlab("Percent Expressing S2 vCyD")+ylab("Relative Drug Effect")+
  ggtitle("5-FC Bystander Effect")+
  theme(plot.title=element_text(hjust=0.5,size=22,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text = element_text(size=18,face="bold",color="black")
  )
# ggsave("BaF3CDBystander.png",width=6.25,height=5)
