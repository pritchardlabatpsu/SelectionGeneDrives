
## Set up
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(reshape2)

## Read in and parse data

df = read.csv("CDCancerLinesIC50.csv",header=T,stringsAsFactors=F)

df.mlt = melt(df,id.vars = c("cell.line","condition","replicate"),
              variable.name = "log.dose",
              value.name = "lum.value")

cell.lines = unique(df.mlt$cell.line)
cndtns = unique(df.mlt$condition)

max.log.dose = 4.5
min.log.dose = 0.5
curve.doses = 10^seq(max.log.dose,min.log.dose,length.out=50)

for (cell.line in cell.lines) {
  
  cell.line.smmry = data.frame()
  
  # Initialize to fit IC50 curves
  
  curves.df = as.data.frame(expand.grid(dose=curve.doses,condition=cndtns))
  curves.df = curves.df[,c(2,1)]
  curves.df$log.dose = log10(curves.df$dose)
  curves.df$pred = NA
  
  for (cndtn in cndtns) {
    
    sub.df = df.mlt[df.mlt$cell.line==cell.line&df.mlt$condition==cndtn,]
    
    # Calculate relative viability
    vc.value = mean(sub.df$lum.value[sub.df$log.dose=="vc"])
    sub.df$rel.viab = sub.df$lum.value/vc.value
    
    # Aggregate data
    sub.agg = aggregate(sub.df$rel.viab,
                        by = list(
                          cell.line = sub.df$cell.line,
                          condition = sub.df$condition,
                          log.dose = sub.df$log.dose
                        ),
                        FUN = function(x) c(mean=mean(x),sem=sd(x)/sqrt(length(x))))
    sub.agg = do.call(data.frame,sub.agg)
    
    sub.agg = sub.agg[sub.agg$log.dose!="vc",]
    sub.agg$log.dose = as.numeric(gsub("X","",sub.agg$log.dose))
    
    cell.line.smmry = rbind(cell.line.smmry,sub.agg)
    
    # Fit IC50 curve
    
    sub.agg$dose = 10^sub.agg$log.dose
    
    # Fit Curve
    fit = dr4pl(x.mean~dose,
                data=sub.agg,
                method.init="logistic",
                init.parm = dr4pl_theta(0.99,100,-2,0.01),
                upperl=c(1,Inf,Inf,0.02),lowerl=c(0.98,-Inf,-Inf,0))
    # plot(fit)
    
    # Predict response
    pred = MeanResponse(coef(fit),curve.doses)
    curves.df$pred[curves.df$condition==cndtn] = pred
    
  }
  
  ## Plot results for each line
  
  cell.line.smmry$condition = factor(cell.line.smmry$condition,levels=c("PIG/5-FC","CD/5-FC","PIG/5-FU"))
  curves.df$condition = factor(curves.df$condition,levels=c("PIG/5-FC","CD/5-FC","PIG/5-FU"))
  
  clrs = c("#2ABAFC","#BA50FF","#8F00FF")
  
  g = ggplot()+theme_bw()+
    geom_line(data=curves.df,aes(x=log.dose,y=pred,color=condition),size=2.5,alpha=0.75)+
    geom_errorbar(data=cell.line.smmry,aes(x=log.dose,ymin=x.mean-x.sem,ymax=x.mean+x.sem),color="black",width=0.15,size=1)+
    geom_point(data=cell.line.smmry,aes(x=log.dose,y=x.mean,fill=condition),size=4,color="black",shape=21)+
    scale_color_manual(values=clrs)+
    scale_fill_manual("",values=clrs,
                      labels = c("WT +5-FC","S2 +5-FC","WT +5-FU"))+
    xlab("Dose (uM)")+ylab("Relative Viability")+
    scale_y_continuous(breaks=c(0,0.5,1),limits=c(0,NA))+
    ggtitle(paste0(cell.line,": S2 vCyD Dose Response"))+
    theme(
      plot.title = element_text(size=22,face="bold",hjust=0.5),
      axis.title = element_text(size=20,face="bold"),
      axis.text = element_text(size=18,face="bold",color="black"),
      legend.title = element_text(size=16,face="bold",color="black"),
      legend.text = element_text(size=14,face="bold",color="black"))+
    guides(color="none")
  
  print(g)
  
  ggsave(paste0(cell.line,"CDIC50.png"),width=7,height=5)  
}

### Run for BaF3 Data (separate experiment) ###

## Set up
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(reshape2)

## Read in and parse data

df = read.csv("BaF3CDIC50_091418.csv",header=T,stringsAsFactors=F)

df.mlt = melt(df,id.vars = c("cell.line","condition","replicate"),
              variable.name = "log.dose",
              value.name = "lum.value")

cell.lines = unique(df.mlt$cell.line)
cndtns = unique(df.mlt$condition)

max.log.dose = 2
min.log.dose = -2
curve.doses = 10^seq(max.log.dose,min.log.dose,length.out=50)

for (cell.line in cell.lines) {
  
  cell.line.smmry = data.frame()
  
  # Initialize to fit IC50 curves
  
  curves.df = as.data.frame(expand.grid(dose=curve.doses,condition=cndtns))
  curves.df = curves.df[,c(2,1)]
  curves.df$log.dose = log10(curves.df$dose)
  curves.df$pred = NA
  
  for (cndtn in cndtns) {
    
    sub.df = df.mlt[df.mlt$cell.line==cell.line&df.mlt$condition==cndtn,]
    
    # Calculate relative viability
    vc.value = mean(sub.df$lum.value[sub.df$log.dose=="vc"])
    sub.df$rel.viab = sub.df$lum.value/vc.value
    
    # Aggregate data
    sub.agg = aggregate(sub.df$rel.viab,
                        by = list(
                          cell.line = sub.df$cell.line,
                          condition = sub.df$condition,
                          log.dose = sub.df$log.dose
                        ),
                        FUN = function(x) c(mean=mean(x),sem=sd(x)/sqrt(length(x))))
    sub.agg = do.call(data.frame,sub.agg)
    
    sub.agg = sub.agg[sub.agg$log.dose!="vc",]
    sub.agg$log.dose = gsub("X\\.","-",sub.agg$log.dose)
    sub.agg$log.dose = gsub("X","",sub.agg$log.dose)
    sub.agg$log.dose = as.numeric(sub.agg$log.dose)
    
    cell.line.smmry = rbind(cell.line.smmry,sub.agg)
    
    # Fit IC50 curve
    
    sub.agg$dose = 10^sub.agg$log.dose
    
    # Fit Curve
    fit = dr4pl(x.mean~dose,
                data=sub.agg,
                method.init="logistic",
                init.parm = dr4pl_theta(0.99,100,-2,0.01),
                upperl=c(1,Inf,Inf,0.02),lowerl=c(0.98,-Inf,-Inf,0))
    # plot(fit)
    
    # Predict response
    pred = MeanResponse(coef(fit),curve.doses)
    curves.df$pred[curves.df$condition==cndtn] = pred
    
  }
  
  ## Plot results for each line
  
  cell.line.smmry$condition = factor(cell.line.smmry$condition,levels=c("PIG/5-FC","CD/5-FC","PIG/5-FU"))
  curves.df$condition = factor(curves.df$condition,levels=c("PIG/5-FC","CD/5-FC","PIG/5-FU"))
  
  clrs = c("#2ABAFC","#BA50FF","#8F00FF")
  
  g = ggplot()+theme_bw()+
    geom_line(data=curves.df,aes(x=log.dose,y=pred,color=condition),size=2.5,alpha=0.75)+
    geom_errorbar(data=cell.line.smmry,aes(x=log.dose,ymin=x.mean-x.sem,ymax=x.mean+x.sem),color="black",width=0.15,size=1)+
    geom_point(data=cell.line.smmry,aes(x=log.dose,y=x.mean,fill=condition),size=4,color="black",shape=21)+
    scale_color_manual(values=clrs)+
    scale_fill_manual(values=clrs)+
    xlab("Dose (uM)")+ylab("Relative Viability")+
    scale_y_continuous(breaks=c(0,0.5,1),limits=c(0,NA))+
    # ggtitle(paste(cell.line,"5-FC/5-FU\nDose Response"))+
    ggtitle("Switch 2 vCyD Dose Response")+
    theme(
      plot.title = element_text(size=22,face="bold",hjust=0.5),
      axis.title = element_text(size=20,face="bold"),
      axis.text = element_text(size=18,face="bold",color="black"))+
    guides(color="none",fill="none")
  
  print(g)
  
  # ggsave(paste0(cell.line,"CDIC50.png"),width=6,height=5)
}

