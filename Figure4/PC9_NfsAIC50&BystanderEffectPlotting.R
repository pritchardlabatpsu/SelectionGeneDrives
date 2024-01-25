
## Set up
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(dr4pl)

### Plot IC50s ###

## Read in data
df = read.csv("MixedPC9_5FC_IC50_20220227_Elvin_SMLSummarized.csv",header=T)

## Plot IC50s

# Focus on 0% and 100% NfsA

sub.df = df[df$perc.GD%in%c(0,100),]
colnames(sub.df)[1] = "condition"
sub.df$dose = 10^sub.df$log.dose.uM

cndtns = c(0,100)

# Fit IC50 Curves

log.max.dose = 4
log.min.dose = 0
curve.doses = 10^seq(log.max.dose,log.min.dose,length.out=50)

curves.df = as.data.frame(expand.grid(dose=curve.doses,condition=cndtns))
curves.df = curves.df[,c(2,1)]
curves.df$log.dose = log10(curves.df$dose)
curves.df$pred = NA

for (cndtn in cndtns) {
  
  # Fit Curve
  fit = dr4pl(mean~dose,
              data=sub.df[sub.df$condition==cndtn,],
              method.init="logistic",
              init.parm = dr4pl_theta(0.99,100,-2,0.01),
              upperl=c(1,Inf,Inf,0.02),lowerl=c(0.98,-Inf,-Inf,0))
  plot(fit)
  
  # Predict response
  pred = MeanResponse(coef(fit),curve.doses)
  curves.df$pred[curves.df$condition==cndtn] = pred
  
}

curves.df$condition = as.factor(curves.df$condition)
sub.df$condition = as.factor(sub.df$condition)

GD.col = "#8F00FF"
WT.col = "#2ABAFC"

ggplot()+theme_bw()+
  geom_line(data=curves.df,aes(x=log.dose,y=pred,color=condition),size=2)+
  geom_errorbar(data=sub.df,aes(x=log.dose.uM,ymin=mean-sem,ymax=mean+sem),color="black",width=0.15,size=1)+
  geom_point(data=sub.df,aes(x=log.dose.uM,y=mean,fill=condition),size=4,color="black",shape=21)+
  scale_color_manual(values=c(WT.col,GD.col))+
  scale_fill_manual(values=c(WT.col,GD.col))+
  xlab("5-FC Dose (uM)")+ylab("Relative Viability")+
  scale_y_continuous(breaks=c(0,0.5,1),limits=c(0,NA))+
  ggtitle("PC9 Complete Gene Drive\n5-FC Dose Response")+
  theme(
    plot.title = element_text(size=22,face="bold",hjust=0.5),
    axis.title = element_text(size=20,face="bold"),
    axis.text = element_text(size=18,face="bold",color="black"))+
  guides(color="none",fill="none")

# ggsave("PC9_CompleteGD_5FCIC50.png",width=5,height=5)

### Plot bystander effect ###

rm(list=ls())
df = read.csv("MixedPC9_5FC_IC50_20220227_Elvin_SMLSummarized.csv",header=T)

# Focus on highest dose

sub.df = df[df$log.dose.uM==3.5&df$perc.GD!=12.5,]

sub.df$killing = 1-sub.df$mean

m = (sub.df$killing[sub.df$perc.GD==100]-sub.df$killing[sub.df$perc.GD==0])/(100-0)
b = sub.df$killing[sub.df$perc.GD==0]

sub.df$no.bystndr = m*sub.df$perc.GD+b

ggplot(sub.df,aes(x=perc.GD,y=killing))+theme_bw()+
  geom_ribbon(aes(ymin=no.bystndr,ymax=killing),fill="#BA50FF",alpha=0.4)+
  geom_line(size=2,color="#8F00FF",)+
  geom_segment(x=0,y=b,xend=100,yend=b+100*m,size=2,color="#8F00FF")+
  geom_errorbar(aes(ymin=killing-sem,ymax=killing+sem),width=4,size=1.5,color="black")+
  geom_point(size=5,shape=21,fill="#8F00FF",color="black")+
  scale_size_manual(values = c(8,5,5,5,8))+
  scale_x_continuous(breaks=c(0,50,100))+
  scale_y_continuous(breaks=c(0,0.3,0.6))+
  xlab("Percent Gene Drive")+ylab("Relative Drug Effect")+
  ggtitle("PC9 Complete Gene Drive\n5-FC Bystander Effect")+
  theme(plot.title=element_text(hjust=0.5,size=22,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text = element_text(size=18,face="bold",color="black")
  )
# ggsave("PC9_CompleteGD_Bystander.png",width=5,height=5)
