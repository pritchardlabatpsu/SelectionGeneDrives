
## Set up
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(reshape2)

## Import data
df = read.csv("PC9BystanderFlow_012822Summarized.csv",header=T,stringsAsFactors=F)
df$condition = factor(df$condition,levels = unique(df$condition))
df$cell.type = factor(df$cell.type,levels=c("CD3+","CD19+","CD19-"))

## Plot data
ggplot(df[df$cell.type%in%c("CD19+","CD19-"),],
       aes(x=condition,y=mean,fill=cell.type))+
  theme_bw()+
  geom_errorbar(aes(ymin=mean-sem,ymax=mean+sem),position=position_dodge(0.9),width=0.5,size=1.25)+
  geom_bar(stat="identity",position=position_dodge(),color="black")+
  scale_fill_manual("",values=c("#8F00FF","#2ABAFC"),labels=c("Ag+","Ag-"))+
  ggtitle("Switch 2 vCD19 Bystander Effect")+
  xlab("Effector:Target Ratio")+
  ylab("Relative Viability")+
  scale_y_continuous(breaks=c(0,0.5,1))+
  theme(
    plot.title = element_text(size=22,face="bold",hjust=0.5),
    axis.title = element_text(size=18,face="bold"),
    axis.text = element_text(size=16,face="bold",color="black"),
    legend.text = element_text(size=16,face="bold",color="black"),
    legend.position = "bottom"
  )

# ggsave("ImmuneBystanderEffect.png",width=5.75,height=5)

ggplot(df[df$cell.type%in%c("CD19+","CD19-"),],
       aes(x=condition,y=mean))+
  theme_bw()+
  geom_errorbar(aes(ymin=mean-sem,ymax=mean+sem,alpha=cell.type),position=position_dodge(0.9),width=0.5,size=1.25)+
  geom_bar(aes(alpha=cell.type,fill=cell.type),stat="identity",position=position_dodge())+
  scale_fill_manual("",values=c("#1AD31A","#2ABAFC"),labels=c("Ag+","Ag-"))+
  scale_alpha_manual(values=c(1,0))+
  ggtitle("T Cell Bystander Effect")+
  xlab("Effector:Target Ratio")+
  ylab("Relative Viability")+
  scale_y_continuous(breaks=c(0,0.5,1))+
  theme(
    plot.title = element_text(size=22,face="bold",hjust=0.5),
    axis.title = element_text(size=18,face="bold"),
    axis.text = element_text(size=16,face="bold",color="black"),
    legend.text = element_text(size=14,face="bold",color="black"),
    legend.position = "bottom"
  )+
  guides(alpha="none")

# ggsave("ImmuneBystanderEffectCD19Only.png",width=5,height=5)

ggplot(df[df$cell.type%in%c("CD19+","CD19-"),],
       aes(x=condition,y=mean))+
  theme_bw()+
  geom_errorbar(aes(ymin=mean-sem,ymax=mean+sem,alpha=cell.type),position=position_dodge(0.9),width=0.5,size=1.25)+
  geom_bar(aes(fill=cell.type),stat="identity",position=position_dodge())+
  scale_fill_manual("",values=c("#b5f7b5","#2ABAFC"),labels=c("Ag+","Ag-"))+
  # scale_alpha_manual(values=c(1,0))+
  ggtitle("T Cell Bystander Effect")+
  xlab("Effector:Target Ratio")+
  ylab("Relative Viability")+
  scale_y_continuous(breaks=c(0,0.5,1))+
  theme(
    plot.title = element_text(size=22,face="bold",hjust=0.5),
    axis.title = element_text(size=18,face="bold"),
    axis.text = element_text(size=16,face="bold",color="black"),
    legend.text = element_text(size=14,face="bold",color="black"),
    legend.position = "bottom"
  )+
  guides(alpha="none")

# ggsave("ImmuneBystanderEffectCD19negEmphasized.png",width=5,height=5)

