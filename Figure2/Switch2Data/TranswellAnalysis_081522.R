
### Analyze results from anti-FasL experiment

## Setup
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(reshape2)

## Import and parse data
counts.df = read.csv("TranswellRawResults_081522.tsv",header=T,stringsAsFactors=F,sep="\t")

counts.df$Tcell = as.numeric(gsub(",","",counts.df$Tcell))
counts.df$CD19.pos = as.numeric(gsub(",","",counts.df$CD19.pos))
counts.df$CD19.neg = as.numeric(gsub(",","",counts.df$CD19.neg))

counts.df$condition = as.numeric(gsub("C","",counts.df$condition))
counts.df$replicate = as.numeric(gsub("R","",counts.df$replicate))

## Summarize data
cnts.agg = aggregate(counts.df[,3:5],
                    by = list(counts.df$condition),
                    FUN = function(x) c(mean=mean(x),sem=sd(x)/sqrt(length(x))))
cnts.agg = do.call(data.frame,cnts.agg)
colnames(cnts.agg)[1] = "condition"

## Visualize data
ggplot(cnts.agg,aes(x=as.factor(condition),y=Tcell.mean))+theme_bw()+
  geom_point()+
  geom_errorbar(aes(ymin=Tcell.mean-Tcell.sem,ymax=Tcell.mean+Tcell.sem))

ggplot(cnts.agg,aes(x=as.factor(condition),y=CD19.pos.mean))+theme_bw()+
  geom_point()+
  geom_errorbar(aes(ymin=CD19.pos.mean-CD19.pos.sem,ymax=CD19.pos.mean+CD19.pos.sem))

ggplot(cnts.agg,aes(x=as.factor(condition),y=CD19.neg.mean))+theme_bw()+
  geom_point()+
  geom_errorbar(aes(ymin=CD19.neg.mean-CD19.neg.sem,ymax=CD19.neg.mean+CD19.neg.sem))


## Quick stats tests
t.test(counts.df$CD19.neg[counts.df$condition==1],counts.df$CD19.neg[counts.df$condition==3])
t.test(counts.df$CD19.neg[counts.df$condition==2],counts.df$CD19.neg[counts.df$condition==4])
t.test(counts.df$CD19.neg[counts.df$condition==1],counts.df$CD19.neg[counts.df$condition==6])


## Calculate viability relative to -blin/T conditions
rel.counts = counts.df[counts.df$condition%in%c(1:4),c("condition","CD19.neg")]
rel.counts$rel.CD19.neg[rel.counts$condition%in%c(1,3)] = rel.counts$CD19.neg[rel.counts$condition%in%c(1,3)]/cnts.agg$CD19.neg.mean[cnts.agg$condition==1]
rel.counts$rel.CD19.neg[rel.counts$condition%in%c(2,4)] = rel.counts$CD19.neg[rel.counts$condition%in%c(2,4)]/cnts.agg$CD19.neg.mean[cnts.agg$condition==2]

rel.agg = aggregate(rel.counts[,3],
                    by = list(rel.counts$condition),
                    FUN = function(x) c(mean=mean(x),sem=sd(x)/sqrt(length(x))))
rel.agg = do.call(data.frame,rel.agg)
colnames(rel.agg)[1] = "condition"

ggplot(rel.agg,aes(x=as.factor(condition),y=x.mean))+theme_bw()+
  geom_errorbar(aes(ymin=x.mean-x.sem,ymax=x.mean+x.sem),width=0.25,size=1.5)+
  geom_bar(stat="identity",fill="#2ABAFC",color="black")+
  scale_y_continuous(breaks=c(0,0.5,1))+
  ylab("Relative Viability\nof Ag- Cells")+
  ggtitle("Immune Bystander\nTranswell Assay")+
  theme(
    plot.title = element_text(size=22,face="bold",hjust=0.5),
    axis.title = element_text(size=20,face="bold"),
    axis.text = element_text(size=18,face="bold",color="black"),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  )

# ggsave("TranswellExperiment.png",width=5,height=4)  
  
## T test

t.test(rel.counts$rel.CD19.neg[rel.counts$condition==3],rel.counts$rel.CD19.neg[rel.counts$condition==4])
