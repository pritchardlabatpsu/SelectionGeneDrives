
## Set up
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)

## Import data 
df = read.csv("Switch2TumorRegressionSummary_030622.csv",header=T,stringsAsFactors=F)

## Plot data
ggplot(df,aes(x=day,y=mean,color=condition))+theme_bw()+
  geom_errorbar(aes(ymin=mean-sem,ymax=mean+sem),
                width=0.25,size=1.5)+
  geom_line(size=3)+
  geom_point(size=5)+
  # geom_bar(stat="identity",fill="#1AD31A",color="black")+
  xlab("Time (days)")+
  ylab("Fold-Change in Tumor Volume")+
  ggtitle("In Vivo Switch 2 Activity")+
  scale_y_continuous(breaks=c(0,5,10,15))+
  scale_x_continuous(limits=c(0,6.5))+
  scale_color_manual(values=c("#2ABAFC","#BA50FF","#8F00FF"))+
  theme(
    plot.title = element_text(size=22,face="bold",hjust=0.5),
    axis.title = element_text(size=20,face="bold"),
    axis.text = element_text(size=18,face="bold",color="black")
  )+
  guides(color="none")

# ggsave("InVivoSwitch2.png",width=8,height=6)

# Adjust plot formatting to match figure
# Plot standard deviation

ggplot(df,aes(x=day,y=mean,fill=condition,color=condition))+theme_bw()+
  geom_errorbar(aes(ymax=mean+sem*sqrt(6),ymin=mean),color="black",width=0.2,size=1.5)+
  geom_line(size=3)+
  geom_point(shape=21,size=5,color="black")+
  # geom_bar(stat="identity",fill="#1AD31A",color="black")+
  xlab("Days")+
  ylab("Fold-Change in Tumor Volume")+
  ggtitle("In Vivo S2 vCyD Activity")+
  scale_y_continuous(breaks=c(0,5,10,15))+
  scale_x_continuous(limits=c(0,6.5))+
  scale_color_manual(values=c("#2ABAFC","#BA50FF","#8F00FF"))+
  scale_fill_manual(values=c("#2ABAFC","#BA50FF","#8F00FF"))+
  theme(
    plot.title = element_text(size=22,face="bold",hjust=0.5),
    axis.title = element_text(size=20,face="bold"),
    axis.text = element_text(size=18,face="bold",color="black")
  )+
  guides(color="none",fill="none")

# ggsave("InVivoSwitch2_2.png",width=6.25,height=4.5)


# Calculate p-values for day 6

no.GD.plus.5FC = c(12.85891835,	17.46426501,	1.840838022,	25.58948097,	11.42758141, 8.163462655,	7.422257312,	5.743723712,	27.92492528,	3.225532362)
GD.no.5FC = c(6.752054501,	7.936517082,	10.08585203,	4.056271216,5.693571545,	3.096948576,	2.487911787,	11.78791659)
GD.plus.5FC = c(0.247689841,	0,	0,	0.084757245,	0.181845903,0.113698403,	0,	0,	0,	0)

t.test(no.GD.plus.5FC,GD.no.5FC) # 0.09
t.test(GD.plus.5FC,GD.no.5FC) # 0.00093

