
###  Plot QC for Twist EGFR library

## Set up

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(pals)

## Import and parse QC Report

qc.df = read.csv("TwistQC.csv",header=T,stringsAsFactors=F)
qc.df$variant_proportion = as.numeric(gsub("%","",qc.df$variant_proportion))/100

## Analyze abundances

qc.df$present = qc.df$variant_proportion>0

variant.number = sum(qc.df$present) # 2717 variants
fraction.present = mean(qc.df$present) # 94% representation

max.proportion = max(qc.df$variant_proportion[qc.df$present])
min.proportion = min(qc.df$variant_proportion[qc.df$present])
fold.diff.prop = max.proportion/min.proportion # 459-fold

## Plot variant proportions

# Fix values so that they sum to one (current proportions don't quite sum to 1 because of rounding errors)
positions = unique(qc.df$customer_aa_position)
qc.df$variant_proportion_scaled = 0 
for (position in positions) {
  
  sum.prop = sum(qc.df$variant_proportion[qc.df$customer_aa_position==position])
  if (sum.prop>0) {
    qc.df$variant_proportion_scaled[qc.df$customer_aa_position==position] = qc.df$variant_proportion[qc.df$customer_aa_position==position]/sum.prop
  }
  
}

ggplot(qc.df,aes(x=customer_aa_position,y=variant_proportion_scaled,fill=variant_aa))+theme_bw()+
  geom_bar(stat="identity")+
  # scale_fill_manual("Variant AA",values=rainbow(20))+
  scale_fill_manual("Variant AA",values=jet(20))+
  ggtitle("EGFR Variant Library")+
  xlab("Amino Acid Position")+
  ylab("Variant Proportion")+
  scale_y_continuous(breaks=c(0,0.5,1))+
  theme(
    plot.title = element_text(hjust=0.5,size=20,face="bold"),
    axis.title = element_text(size=18,face="bold"),
    axis.text = element_text(size=16,face="bold",color="black"),
    legend.title = element_text(size=14,face="bold"),
    legend.text = element_text(size=12,face="bold",color="black"),
    legend.position = "bottom"
  )
# ggsave("LibraryQC.png",width=8,height=6)

## Plot histogram

ggplot(qc.df,aes(x=variant_proportion*100))+theme_bw()+
  geom_histogram(binwidth=1,fill="#7EA086")+
  xlab("Variant Proportion (%)")+
  ylab("Count")+
  ggtitle("EGFR Library Variant Distribution")+
  theme(
    plot.title = element_text(hjust=0.5,size=20,face="bold"),
    axis.title = element_text(size=18,face="bold"),
    axis.text = element_text(size=16,face="bold",color="black")
  )
# ggsave("EGFRVariantDistribution.png",width=7,height=5)

mean(qc.df$variant_proportion<0.1)
