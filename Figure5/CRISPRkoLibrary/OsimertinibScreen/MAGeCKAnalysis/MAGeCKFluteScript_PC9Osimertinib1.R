
## Counts generated using PoolQScript_PC9Osimertinib1.R
## Using MAGeCK MLE to generate LFCs and p-values (following Wang et al Nat Prot 2019)

## Set up
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(MAGeCKFlute)
library(ggplot2)
library(ggrepel)

## Create and activate MAGeCK-VISPR environment - not sure if this will work in R, just ran commands on terminal
# system("conda create -n mageck-vispr mageck mageck-vispr python=3")
# system("source activate mageck-vispr")

## Generate input files for mageck mle command

samples = c("D0","DMSO_rep1","DMSO_rep2","osi_rep1","osi_rep2")

# Parse count table from PoolQ
counts.df = read.csv("../PoolQResults/counts.txt",sep="\t",stringsAsFactors=F)
colnames(counts.df) = c("sgRNA","Gene",samples)
counts.df$Gene[counts.df$Gene=="Non-Targeting Control"] = "NonTargeting" # MAGeCK does not play well with the dashes
write.table(counts.df,"GuideCounts.txt",sep="\t",row.names=F,quote=F)

# Create design matrix (see paper for explanation of format)
design.mat = data.frame(
  samples = samples,
  baseline = rep(1,5),
  DMSO = c(0,1,1,0,0),
  osi = c(0,0,0,1,1)
)
write.table(design.mat,"DesignMatrix.txt",sep="\t",row.names=F,quote=F)

# Create control guides file
# Note: previous version of guides file had some errors introduced when saving it in Excel (e.g. mistaking gene names as dates)
guides.df = read.csv("../BrunelloGuides.tsv",header=T,sep="\t",stringsAsFactors=F)
guides.df = guides.df[,c("Target.Gene.Symbol","sgRNA.Target.Sequence")]
NTC.guides = guides.df$sgRNA.Target.Sequence[guides.df$Target.Gene.Symbol=="Non-Targeting Control"]
write.table(NTC.guides,"NonTargetingControls.txt",row.names=F,col.names=F,quote=F)

# Ran MAGeCK MLE on the terminal (this present directory)
system("mageck mle --count-table GuideCounts.txt --design-matrix DesignMatrix.txt --norm-method control --control-sgrna NonTargetingControls.txt --output-prefix MLEResult") 
# Note: also trying a separate run where we set max target per gene to 1000 so that it can handle the nontargeting genes
# We'll see if that makes any difference in the analysis
system("mageck mle --count-table GuideCounts.txt --design-matrix DesignMatrix.txt --norm-method control --control-sgrna NonTargetingControls.txt --max-sgrnapergene-permutation 1000 --output-prefix MLEResultWithNTCs") 
# jk that didn't work - we'll stick with the original results

## MAGeCKFlute Analysis

FluteMLE(gene_summary="MLEResult.gene_summary.txt", ctrlname="DMSO", treatname="osi", organism="hsa")
# Okay, that wasn't terribly helpful

# Okay, so it looks like there's no way to get a p-value of differences in LFC values between conditions from MAGeCK MLE
# But mathematically, differences in LFC values (i.e. log(cond1/baseline)-log(cond2/baseline)) should be the same as the LFC of condition1 vs condition2 (log(condition1/condition2))
# Which is what MAGeCK RRA calculates, with appropriate p-values, so let's just use that

## MAGeCK RRA

# Remove D0 counts from counts file
cnts.ctrl.trt = counts.df[,-which(colnames(counts.df)=="D0")]
write.table(counts.df,"GuideCountsCtrlTrt.txt",sep="\t",row.names=F,quote=F)

system("mageck test -k GuideCountsCtrlTrt.txt -t osi_rep1,osi_rep2 -c DMSO_rep1,DMSO_rep2 -n RRAResult")

## Visualize results
RRA.results = ReadRRA("RRAResult.gene_summary.txt")
VolcanoView(RRA.results[RRA.results$Score>0,], x = "Score", y = "FDR", Label = "id",x_cutoff=1)+
  theme_bw()+
  scale_x_continuous(limits = c(-0.15,3.75))+
  scale_y_continuous(limits = c(-0.15,3.5))+
  theme(
    axis.text = element_text(color="black",face="bold")
  )

# Custom with ggplot2

LFC.cutoff = 1
pval.cutoff = 0.05

RRA.results$log10.adjp = -log10(RRA.results$FDR)
RRA.results$LFC.cutoff = RRA.results$Score>LFC.cutoff
RRA.results$pval.cutoff = RRA.results$log10.adjp > -log10(0.05)
RRA.results$hits = RRA.results$LFC.cutoff & RRA.results$pval.cutoff
RRA.results$label = RRA.results$hits & RRA.results$log10.adjp>2

max.hit = "KCTD5"
max.score = RRA.results$Score[RRA.results$id==max.hit]
max.pval = RRA.results$log10.adjp[RRA.results$id==max.hit]
RRA.results$dist2maxhit = sqrt((RRA.results$Score-max.score)^2+(RRA.results$log10.adjp-max.pval)^2)

ggplot(RRA.results[RRA.results$Score>0,],aes(x=Score,y=log10.adjp,color=dist2maxhit))+theme_bw()+
  geom_point(size=3.5,alpha=0.8)+
  ylim(0,3.5)+
  geom_label_repel(data=RRA.results[RRA.results$label,],aes(label=id),fontface="bold")+
  scale_color_gradient(high="#2ABAFC",low="#EE1D23")+
  xlab("Log2 Fold Change")+
  ylab("-Log10 Adjusted Pvalue")+
  ggtitle("Genome-wide Osimertinib Screen")+
  theme(
    plot.title = element_text(size=18,face="bold",hjust=0.5),
    axis.title = element_text(size=18,face="bold"),
    axis.text = element_text(size=16,face="bold",color="black")
  )+
  guides(color="none")

# ggsave("VolcanoPlot.png",width=5,height=5)
# ggsave("VolcanoPlot_021523.png",width=8,height=5)

