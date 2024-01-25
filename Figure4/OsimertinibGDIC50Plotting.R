
## Set up
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(dr4pl)

## Import and parse data
df = read.csv("HMEC_PC9C2IC50_100222Summarized.csv",header=T,stringsAsFactors=F)
df$dose = 10^df$log.dose.nM

df$condition = paste(df$genotype,"-dim")
df$condition[df$dimerizer] = paste(df$genotype[df$dimerizer],"+dim")

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

df$condition = factor(df$condition,levels = c("WT -dim","ELC -dim","HMEC -dim","HMEC +dim"))
clrs = c("#2ABAFC","#EE1D23","#7C7C7C","#F99900")

ggplot()+theme_bw()+
  # geom_line(data=curves.df,aes(x=log.dose,y=pred,color=condition),size=2.5,alpha=0.5)+
  geom_line(data=df,aes(x=log.dose.nM,y=mean,color=condition),stat="smooth",size=2.5,alpha=0.75,se=F)+
  geom_errorbar(data=df,aes(x=log.dose.nM,ymin=mean-sem,ymax=mean+sem),color="black",width=0.15,size=1)+
  geom_point(data=df,aes(x=log.dose.nM,y=mean,fill=condition),size=4,color="black",shape=21)+
  scale_fill_manual("",values=clrs,
                    labels = c("WT","EGFR L858R/C797S"," -Dim","+Dim"))+
  scale_color_manual(values=clrs)+
  xlab("Osimertinib Dose (nM)")+ylab("Relative Viability")+
  scale_y_continuous(breaks=c(0,0.5,1))+
  ggtitle("PC9 Complete Gene Drive\nOsimertinib Dose Response")+
  theme(
    plot.title = element_text(size=22,face="bold",hjust=0.5),
    axis.title = element_text(size=20,face="bold"),
    axis.text = element_text(size=18,face="bold",color="black"),
    legend.title = element_text(size=16,face="bold",color="black"),
    legend.text = element_text(size=14,face="bold",color="black"),
    legend.spacing.y = unit(0.2,'in'))+
  guides(color="none",
         fill = guide_legend(byrow = TRUE))

# ggsave("PC9OsimertinibGDIC50.png",width=8,height=5)

