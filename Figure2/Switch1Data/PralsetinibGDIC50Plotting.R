
## Set up
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(dr4pl)

## Import data
df = read.csv("Pralsetinib_GDIC50Summarized_042322.csv",header=T,stringsAsFactors=F)
df$dose = 10^df$log.dose

## Fit IC50 curves

cndtns = unique(df$condition)

max.dose = max(df$dose)
min.dose = min(df$dose)
curve.doses = 10^seq(log10(max.dose),log10(min.dose),length.out=50)

curves.df = as.data.frame(expand.grid(dose=curve.doses,condition=cndtns))
curves.df = curves.df[,c(2,1)]
curves.df$log.dose = log10(curves.df$dose)
curves.df$pred = NA

for (cndtn in cndtns) {
  
  # Fit Curve
  fit = dr4pl(mean~dose,
              data=df[df$condition==cndtn,],
              method.init="logistic",
              init.parm = dr4pl_theta(0.99,100,-2,0.01),
              upperl=c(1,Inf,Inf,0.02),lowerl=c(0.98,-Inf,-Inf,0))
  # plot(fit)
  
  # Predict response
  pred = MeanResponse(coef(fit),curve.doses)
  curves.df$pred[curves.df$condition==cndtn] = pred
  
}

df$condition = factor(df$condition,levels = c("WT","G810R","GD -dim","GD +dim"))
curves.df$condition = factor(curves.df$condition,levels = c("WT","G810R","GD -dim","GD +dim"))

ggplot()+theme_bw()+
  geom_line(data=curves.df,aes(x=log.dose,y=pred,color=condition),size=2.5,alpha=0.75)+
  geom_errorbar(data=df,aes(x=log.dose,ymin=mean-sem,ymax=mean+sem),color="black",width=0.15,size=1)+
  geom_point(data=df,aes(x=log.dose,y=mean,fill=condition),size=5,color="black",shape=21)+
  scale_color_manual(values=c("#2ABAFC","#EE1D23","#7C7C7C","#F99900"))+
  scale_fill_manual("",values=c("#2ABAFC","#EE1D23","#7C7C7C","#F99900"),
                    labels = c("WT","Resistant","S1 -Dim","S1 +Dim"))+
  xlab("Pralsetinib Dose (nM)")+ylab("Relative Viability")+
  ggtitle("TPC1 S1 vRETprals Dose Response")+
  theme(
    plot.title = element_text(size=22,face="bold",hjust=0.5),
    axis.title = element_text(size=20,face="bold"),
    axis.text = element_text(size=18,face="bold",color="black"),
    legend.title = element_text(size=16,face="bold",color="black"),
    legend.text = element_text(size=14,face="bold",color="black"))+
  guides(color="none")

# ggsave("TPC1_Switch1IC50.png",width=7,height=5)

