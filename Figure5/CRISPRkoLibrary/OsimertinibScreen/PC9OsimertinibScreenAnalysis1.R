
## Set up
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(reshape2)
library(viridis)

## Import and parse data

# Read in count data
counts.df = read.csv("PoolQResults/counts.txt",sep="\t",header=T,stringsAsFactors=F)

# Check for coverage
plot(colSums(counts.df[3:ncol(counts.df)]))

## Calculate LFC values

cndtns = colnames(counts.df[3:ncol(counts.df)])

# Normalize reads

nrmlz = function(x) {
  nrm = log2(x/sum(x)*1e6+1)
  return(nrm)
}

norm.df = cbind(counts.df[,1:2],
                apply(counts.df[,cndtns],2,nrmlz))


# Calculate LFC relative to D0 (for appropriate replicate)

LFC.df = norm.df[,1:2]
trtmnts = cndtns[!grepl("D0",cndtns)]

for (trtmnt in trtmnts) {
  
  LFC.df[,trtmnt] = norm.df[,trtmnt] - norm.df$D0
  
}


# QC for DMSO controls

DMSO.df = LFC.df[,grepl("DMSO|Row",colnames(LFC.df))]

# Label guides essential, nonessential, conditionally essential, or non-targeting
DMSO.df$essentiality = "conditionally essential"
DMSO.df$essentiality[DMSO.df$Row.Barcode.IDs=="Non-Targeting Control"] = "NTC"

nonessntl.df = read.csv("nonessentials.csv",header=T,stringsAsFactors=F)
nonessentials = unlist(lapply(strsplit(nonessntl.df$gene," "),'[[',1))
DMSO.df$essentiality[DMSO.df$Row.Barcode.IDs%in%nonessentials] = "nonessential"

essntl.df = read.csv("common_essentials.csv",header=T,stringsAsFactors=F)
essentials = unlist(lapply(strsplit(essntl.df$gene," "),'[[',1))
DMSO.df$essentiality[DMSO.df$Row.Barcode.IDs%in%essentials] = "essential"

DMSO.mlt = melt(DMSO.df,id.vars = c("Row.Barcode","Row.Barcode.IDs","essentiality"),
                variable.name="condition",value.name="LFC")

DMSO.mlt$essentiality = factor(DMSO.mlt$essentiality,levels=c("essential","conditionally essential","nonessential","NTC"))

ggplot(DMSO.mlt[grepl("rep1",DMSO.mlt$condition)&DMSO.mlt$essentiality!="conditionally essential",],aes(x=LFC,fill=essentiality))+theme_bw()+
  geom_histogram(position="identity",alpha=0.8)+
  scale_fill_manual("Essentiality",values=c("#ED6A5A","#9BC1BC","#5D576B"),labels=c("Common Essential","Common Nonessential","Non-targeting Control"))+
  ggtitle("DMSO Replicate 1")+
  ylab("Count")+
  theme(
    plot.title = element_text(size=20,hjust=0.5,face="bold"),
    axis.text = element_text(size=16,face="bold",color="black"),
    axis.title = element_text(size=18,face="bold"),
    legend.title = element_text(size=18,face="bold"),
    legend.text = element_text(size=12,face="bold")
  )
# ggsave("PC9BrunelloScreen_NNMDrep1.png",width=7,height=5)
ggplot(DMSO.mlt[grepl("rep2",DMSO.mlt$condition)&DMSO.mlt$essentiality!="conditionally essential",],aes(x=LFC,fill=essentiality))+theme_bw()+
  geom_histogram(position="identity",alpha=0.8)+
  scale_fill_manual("Essentiality",values=c("#ED6A5A","#9BC1BC","#5D576B"),labels=c("Common Essential","Common Nonessential","Non-targeting Control"))+
  ggtitle("DMSO Replicate 2")+
  ylab("Count")+
  theme(
    plot.title = element_text(size=20,hjust=0.5,face="bold"),
    axis.text = element_text(size=16,face="bold",color="black"),
    axis.title = element_text(size=18,face="bold"),
    legend.title = element_text(size=18,face="bold"),
    legend.text = element_text(size=12,face="bold")
  )
# ggsave("PC9BrunelloScreen_NNMDrep2.png",width=7,height=5)


cndtns = c("rep1","rep2")
mean.essntl = rep(NA,2)
names(mean.essntl) = cndtns
mean.NTC = mean.guide
sd.NTC = mean.guide
for (i in 1:length(cndtns)) {
  cndtn = cndtns[i]
  mean.essntl[i] = mean(DMSO.df[DMSO.df$essentiality=="essential",grep(cndtn,colnames(DMSO.df),value=T)])
  mean.NTC[i] = mean(DMSO.df[DMSO.df$essentiality=="NTC",grep(cndtn,colnames(DMSO.df),value=T)])
  sd.NTC[i] = sd(DMSO.df[DMSO.df$essentiality=="NTC",grep(cndtn,colnames(DMSO.df),value=T)])
}

NNMD = (mean.essntl-mean.NTC)/sd.NTC
# NNMD = -2.4 for both replicates, indicating good separation of controls


# Calculate average LFC for each gene

genes = unique(LFC.df$Row.Barcode.IDs)
genes = genes[genes!="Non-Targeting Control"]

LFC.avg = data.frame(gene = genes)

for (gene in genes) {
  for (trtmnt in trtmnts) {
    
    avg = mean(LFC.df[LFC.df$Row.Barcode.IDs==gene,trtmnt])
    LFC.avg[LFC.avg$gene==gene,trtmnt] = avg
    
  }
}

# Check correlations

cor(LFC.avg$DMSO_rep1,LFC.avg$DMSO_rep2) # 0.63
cor(LFC.avg$osi_rep1,LFC.avg$osi_rep2) # 0.60

# Plot correlations

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

LFC.avg$DMSO.density = get_density(LFC.avg$DMSO_rep1,LFC.avg$DMSO_rep2,n=100)

ggplot(LFC.avg,aes(x=DMSO_rep1,y=DMSO_rep2,color=DMSO.density))+theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",color="gray")+
  geom_point()+
  scale_color_viridis()+
  xlab("Replicate 1 LFC")+ylab("Replicate 2 LFC")+
  ggtitle("DMSO Conditions")+
  xlim(-2.75,1.25)+ylim(-2.75,1.25)+
  theme(
    plot.title = element_text(size=20,hjust=0.5,face="bold"),
    axis.text = element_text(size=16,face="bold",color="black"),
    axis.title = element_text(size=18,face="bold"),
    legend.title = element_text(size=18,face="bold"),
    legend.text = element_text(size=12,face="bold")
  )+
  guides(color="none")
# ggsave("DMSOCorrelationPlot.png",width=5,height=5)

LFC.avg$osi.density = get_density(LFC.avg$osi_rep1,LFC.avg$osi_rep2,n=100)

ggplot(LFC.avg,aes(x=osi_rep1,y=osi_rep2,color=osi.density))+theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",color="gray")+
  geom_point()+
  scale_color_viridis()+
  xlab("Replicate 1 LFC")+ylab("Replicate 2 LFC")+
  ggtitle("Osimertinib Conditions")+
  xlim(-3.5,4)+ ylim(-3.5,4)+
  theme(
    plot.title = element_text(size=20,hjust=0.5,face="bold"),
    axis.text = element_text(size=16,face="bold",color="black"),
    axis.title = element_text(size=18,face="bold"),
    legend.title = element_text(size=18,face="bold"),
    legend.text = element_text(size=12,face="bold")
  )+
  guides(color="none")
# ggsave("OsimertinibCorrelationPlot.png",width=5,height=5)


# For now, collapse LFCs by averaging

LFC.avg$DMSO = rowMeans(LFC.avg[,c("DMSO_rep1","DMSO_rep2")])
LFC.avg$osi = rowMeans(LFC.avg[,c("osi_rep1","osi_rep2")])

# Difference in LFCs

LFC.avg$diff = LFC.avg$osi-LFC.avg$DMSO

plot(LFC.avg$diff)

top.hits = LFC.avg$gene[order(LFC.avg$diff)]
