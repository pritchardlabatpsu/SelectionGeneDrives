
## Set up
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(dr4pl)
library(painter)

## Import data
df = read.csv("BypassIC50sIvan.csv",header=T,stringsAsFactors=F)

df$genotype[df$genotype=="pUltra"] = "WT"

## Normalize
# Unfortunately, Ivan didn't run an untreated control here
# So we'll normalize to the lowest dose (0.316 nM)

min.dose = min(df$concentration.uM)
gens = unique(df$genotype)
df$norm = NA

for (i in 1:length(gens)) {
  
  gen.i = gens[i]
  
  low.dose.value = mean(df$value[df$genotype==gen.i&df$concentration.uM==min.dose])
  df$norm[df$genotype==gen.i] = df$value[df$genotype==gen.i]/low.dose.value
  
}

## Aggregate data

agg.df = aggregate(df$norm,
                    by = list(condition=df$genotype,
                              dose=df$concentration.uM),
                    FUN = function(x) c(mean=mean(x),sem=sd(x)/sqrt(length(x))))
agg.df = do.call(data.frame,agg.df)

agg.df$log.dose = log10(agg.df$dose)

## Fit IC50 curves

cndtns = unique(agg.df$condition)

max.dose = max(agg.df$dose)
min.dose = min(agg.df$dose)
curve.doses = 10^seq(log10(max.dose),log10(min.dose),length.out=50)

curves.df = as.data.frame(expand.grid(dose=curve.doses,condition=cndtns))
curves.df = curves.df[,c(2,1)]
curves.df$log.dose = log10(curves.df$dose)
curves.df$pred = NA

IC50s = rep(NA,length(cndtns))
names(IC50s) = cndtns

for (cndtn in cndtns) {
  
  # Fit Curve
  fit = dr4pl(x.mean~dose,
              data=agg.df[agg.df$condition==cndtn,],
              method.init="logistic",
              init.parm = dr4pl_theta(0.99,100,-2,0.01),
              upperl=c(1,Inf,Inf,0.02),lowerl=c(0.98,-Inf,-Inf,0))
  # plot(fit)
  
  # Predict response
  pred = MeanResponse(coef(fit),curve.doses)
  curves.df$pred[curves.df$condition==cndtn] = pred
  
  IC50s[cndtn] = coef(fit)[2]
  
}

# Order by IC50
agg.df$condition = factor(agg.df$condition,levels = names(IC50s[order(IC50s)]))
curves.df$condition = factor(curves.df$condition,levels = names(IC50s[order(IC50s)]))

clrs = c("#2ABAFC",Palette("#f58e91","#EE1D23",7,mode="RGB",circular=F))

ggplot()+theme_bw()+
  geom_line(data=curves.df,aes(x=log.dose,y=pred,color=condition),size=2,alpha=0.4)+
  geom_errorbar(data=agg.df,aes(x=log.dose,ymin=x.mean-x.sem,ymax=x.mean+x.sem,color=condition),width=0.15,size=1)+
  geom_point(data=agg.df,aes(x=log.dose,y=x.mean,color=condition),size=4)+
  scale_color_manual("Genotype",values=clrs)+
  xlab("Osimertinib Dose (nM)")+ylab("Relative Viability")+
  ggtitle("Bypass Mutation Dose Response")+
  theme(
    plot.title = element_text(size=24,face="bold",hjust=0.5),
    axis.title = element_text(size=20,face="bold"),
    axis.text = element_text(size=18,face="bold",color="black"),
    legend.title = element_text(size=20,face="bold",color="black",hjust=0.5),
    legend.text = element_text(size=12,face="bold",color="black")
    )

# ggsave("BypassLibraryIC50s.png",width=8*0.85,height=6*0.85)

## Alternatively, plot viability just as bar plot (dose 100 nM)

agg.100 = agg.df[agg.df$dose==0.1,]
agg.100$condition = factor(agg.100$condition,levels = agg.100$condition[order(agg.100$x.mean)]) # order by viability

ggplot(agg.100,aes(x=condition,y=x.mean,fill=condition,color=condition))+theme_bw()+
  geom_errorbar(aes(ymin=x.mean-x.sem,ymax=x.mean+x.sem),width=0.25,size=1,color="black")+
  geom_bar(stat="identity",color="black")+
  scale_fill_manual("Genotype",values=clrs)+
  scale_color_manual("Genotype",values=clrs)+
  xlab("Genotype")+ylab("Relative Viability")+
  ggtitle("Bypass Oncogenes Drug Response")+
  theme(
    plot.title = element_text(size=20,face="bold",hjust=0.5),
    axis.title = element_text(size=18,face="bold"),
    axis.text = element_text(size=14,face="bold",color="black"),
    axis.text.x = element_text(angle=45,hjust=1)
  )+
  guides(
    color="none",
    fill="none"
  )
  
# ggsave("BypassLibraryViabilityAt100nM.png",width=6,height=5)
  
