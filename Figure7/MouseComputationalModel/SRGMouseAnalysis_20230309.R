### Analyze SRG Simulation Results
### 9/23/21

## Set up

rm(list=ls())
library(ggplot2)
library(reshape2)
library(viridis)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Import and summarize files

csv.files = list.files(pattern = "\\.csv$")

smmry.df = as.data.frame(matrix(nrow=length(csv.files),ncol=5))
colnames(smmry.df) = c("gd.frac","swtch.dly","p.erad","med.PFS","med.t2erad")

for (i in 1:length(csv.files)) {
  
  file.i = csv.files[i]
  fname.splt = strsplit(file.i,"_")[[1]]
  
  smmry.df$gd.frac[i] = 10^(-as.numeric(gsub("negloggdfrac","",fname.splt[2])))
  smmry.df$swtch.dly[i] = as.numeric(gsub("swtchdly","",fname.splt[3]))
  
  
  df.i = read.csv(file.i,header=F)
  colnames(df.i) = c("erad","t2end")
  df.i = df.i[!is.na(df.i$erad),]
  
  smmry.df$p.erad[i] = mean(df.i$erad)
  
  PFS = df.i$t2end[!df.i$erad]
  smmry.df$med.PFS[i] = median(PFS)
  
  t2erad = df.i$t2end[df.i$erad]
  smmry.df$med.t2erad[i] = median(t2erad)
  
  print(i)
  
}

gd.fracs = unique(smmry.df$gd.frac) # gotta call the gene drive fractions because they were log-transformed and the unlog-transformed in the files, so values have some rounding error

## Visualize summarized data

# sub.df = smmry.df[smmry.df$poptreat==1e9&smmry.df$mut==1e-8&smmry.df$swtchdly==365*1&smmry.df$a_B==0.04&smmry.df$trnovr==14&smmry.df$ngr==0.01,]

# |smmry.df$gd.frac==0

ggplot(smmry.df[smmry.df$gd.frac>0.1,],aes(x=as.factor(gd.frac),y=swtch.dly,fill=med.PFS))+theme_bw()+
  geom_tile()+
  ggtitle("In Vivo Treatment\nSchedule Optimization")+
  xlab("Initial Gene Drive Proportion")+
  ylab("Switch Overlap (days)")+
  scale_fill_gradientn("Time to\nFailure",colors=viridis(10),breaks=c(75,100,125))+
  scale_x_discrete(breaks=gd.fracs[c(1,3,5,7,11)],labels = c("50%","40%","30%","20%","0%"))+
  # scale_y_continuous(breaks=c(-0.14,0),labels=c("",""))+
  theme(
    plot.title = element_text(size=22,hjust=0.5,face="bold"),
    axis.text = element_text(size=16,face="bold",color="black"),
    axis.title = element_text(size=18,face="bold"),
    legend.title = element_text(size=16,face="bold",hjust=0.5),
    legend.text = element_text(size=12,face="bold")
  )

# ggsave("TreatmentScheduleOptimization.png",width=6,height=5)
