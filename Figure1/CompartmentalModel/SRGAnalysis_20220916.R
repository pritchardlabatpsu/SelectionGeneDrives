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

smmry.df = as.data.frame(matrix(nrow=length(csv.files),ncol=11))
colnames(smmry.df) = c("poptreat","mut","infeff","swtchdly","dimres","a_B","trnovr","ngr","p.erad","med.PFS","med.t2erad")

for (i in 1:length(csv.files)) {
  
  file.i = csv.files[i]
  fname.splt = strsplit(file.i,"_")[[1]]
  
  smmry.df$poptreat[i] = 10^as.numeric(gsub("logpoptreat","",fname.splt[2]))
  smmry.df$mut[i] = 10^(-as.numeric(gsub("neglogmut","",fname.splt[3])))
  smmry.df$infeff[i] = 10^(-as.numeric(gsub("negloginfeff","",fname.splt[4])))
  smmry.df$swtchdly[i] = as.numeric(gsub("swtchdly","",fname.splt[5]))
  smmry.df$dimres[i] = as.numeric(gsub("dimres","",fname.splt[6]))
  smmry.df$a_B[i] = as.numeric(gsub("aB","",fname.splt[7]))
  smmry.df$trnovr[i] = as.numeric(gsub("turnovr","",fname.splt[8]))
  smmry.df$ngr[i] = as.numeric(gsub("ngr","",fname.splt[9]))
  
  
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

## Translate dimres parameter into fitness (net growth rate) relative to native resistant population
# The variable dimres scales the targeted therapy kill rate (a_T) in the original simulation.
# The overall drug kill rate for the gene drive population is defined as a_T*dim_res, so
# dim_res = 0 is completely resistant and
# dim_res = 1 is completely sensitive
# for the parameters used, a neutral fitness (no net growth/death) was dim_res = 0.25
# So gene drive cells would retain a positive net growth for dim_res between 0 and 0.25
# We'd like to present it in terms of fitness relative to resistant population.
# Resistant cells in drug have a net growth of (b-d-a_T) = (0.14-0.13-0) = 0.01 (default)
# Gene drive cells in drug have a net growth of (b-d-a_T*dim_res) = (0.14-0.13-0.04*dim_res)
# So their relative net growth rate is (0.01-0.04*dim_res)/0.01
# or 1 - (0.04*dim_res)/0.01 = 1 - 4*dim_res

# More generally, gene drive relative netgrowth will be: (ngr-a_T*dim_res)/ngr = 1-a_T/ngr*dim_res

# smmry.df$gd.relfitness = 1 - 4*smmry.df$dimres
smmry.df$gd.relfitness = 1 - 0.04/smmry.df$ngr*smmry.df$dimres

## Visualize summarized data

sub.df = smmry.df[smmry.df$poptreat==1e9&smmry.df$mut==1e-8&smmry.df$swtchdly==365*1&smmry.df$a_B==0.04&smmry.df$trnovr==14&smmry.df$ngr==0.01,]

ggplot(sub.df,aes(x=as.factor(infeff),y=gd.relfitness,fill=p.erad))+theme_bw()+
  geom_tile()+
  ggtitle("Compartmental Model Results")+
  xlab("Gene Delivery Efficiency")+
  ylab("Gene Drive Fitness\n(Relative to Native Resistance)")+
  scale_fill_gradientn("",colors=viridis(10),breaks=c(0,1),limits=c(0,1))+
  scale_x_discrete(breaks=c(0,.001,0.01,0.1),labels=parse(text=c("0","10^-3","10^-2","10^-1")))+
  # scale_y_continuous(breaks=c(-0.14,0),labels=c("",""))+
  theme(
    plot.title = element_text(size=22,hjust=0.5,face="bold"),
    axis.text = element_text(size=16,face="bold",color="black"),
    axis.title = element_text(size=18,face="bold")
  )

# ggsave("CompartmentalModelResults.png",width=5.5,height=5)

## Analyze population structure at start of treatment

pop.struct.files = list.files(path="../SRG_PopulationStructure_20220915/",pattern = "\\.csv$")

pop.smmry.df = as.data.frame(matrix(nrow=length(pop.struct.files),ncol=5))
colnames(pop.smmry.df) = c("infeff","dimres","mean.gd","mean.res","gd2res.ratio")

for (i in 1:length(pop.struct.files)) {
  
  file = pop.struct.files[i]
  
  fname.splt = strsplit(file,"_")[[1]]
  pop.smmry.df$infeff[i] = 10^(-as.numeric(gsub("negloginfeff","",fname.splt[4])))
  pop.smmry.df$dimres[i] = as.numeric(gsub("dimres","",fname.splt[6]))
  
  df.i = read.csv(paste0("../SRG_PopulationStructure_20220915/",file),header=F)
  colnames(df.i) = c("time","S","R","M","D","G","Q","L","C")
  
  df.i$res.pops = rowSums(df.i[,c("R","M","D","Q","L","C")])
  
  pop.smmry.df$mean.gd[i] = mean(df.i$G)
  pop.smmry.df$mean.res[i] = mean(df.i$res.pops)
  pop.smmry.df$gd2res.ratio[i] = mean(df.i$G/df.i$res.pops)
  
}

mean(pop.smmry.df$gd2res.ratio[pop.smmry.df$infeff==max(pop.smmry.df$infeff)]) # 21882.4
mean(pop.smmry.df$gd2res.ratio[pop.smmry.df$infeff==0.001]) # 96.3

## Sensitivity Analysis

base.pars = list(
  poptreat = 1e9,
  mut = 1e-8,
  swtchdly = 365,
  a_B = 0.04,
  trnovr = 14,
  ngr = 0.01
)

min.PFS = min(smmry.df$med.PFS,na.rm=T)
max.PFS = max(smmry.df$med.PFS,na.rm=T)

# Vary pop treat

pop.treat.vec = unique(smmry.df$poptreat)

for (i in 1:length(pop.treat.vec)) {
  
  poptreat.i = pop.treat.vec[i]
  
  sub.df = smmry.df[smmry.df$poptreat==poptreat.i&
                      smmry.df$mut==base.pars$mut&
                      smmry.df$swtchdly==base.pars$swtchdly&
                      smmry.df$a_B==base.pars$a_B&
                      smmry.df$trnovr==base.pars$trnovr&
                      smmry.df$ngr==base.pars$ngr,]
  
  ggplot(sub.df,aes(x=as.factor(infeff),y=gd.relfitness,fill=p.erad))+theme_bw()+
    geom_tile()+
    ggtitle(paste("Detection Size:",poptreat.i,"Cells"))+
    xlab("Gene Delivery Efficiency")+
    ylab("Gene Drive Fitness\n(Relative to Native Resistance)")+
    scale_fill_gradientn("Erad\nProb",colors=viridis(10),breaks=c(0,1),limits=c(0,1))+
    scale_x_discrete(breaks=c(0,.001,0.01,0.1),labels=parse(text=c("0","10^-3","10^-2","10^-1")))+
    # scale_y_continuous(breaks=c(-0.14,0),labels=c("",""))+
    theme(
      plot.title = element_text(size=22,hjust=0.5,face="bold"),
      axis.text = element_text(size=16,face="bold",color="black"),
      axis.title = element_text(size=18,face="bold"),
      legend.title = element_text(size=14,face="bold",hjust = 0.5),
      legend.text = element_text(size=12,face="bold")
    )
    # ggsave(paste0("SensitivityAnalysisFigs/EradicationPoptreat",poptreat.i,".png"),width=5.5,height=5)
    
    ggplot(sub.df,aes(x=as.factor(infeff),y=gd.relfitness,fill=med.PFS))+theme_bw()+
      geom_tile()+
      ggtitle(paste("Detection Size:",poptreat.i,"Cells"))+
      xlab("Gene Delivery Efficiency")+
      ylab("Gene Drive Fitness\n(Relative to Native Resistance)")+
      scale_fill_gradientn("Median\nPFS",colors=viridis(10),limits=c(min.PFS,max.PFS))+
      scale_x_discrete(breaks=c(0,.001,0.01,0.1),labels=parse(text=c("0","10^-3","10^-2","10^-1")))+
      # scale_y_continuous(breaks=c(-0.14,0),labels=c("",""))+
      theme(
        plot.title = element_text(size=22,hjust=0.5,face="bold"),
        axis.text = element_text(size=16,face="bold",color="black"),
        axis.title = element_text(size=18,face="bold"),
        legend.title = element_text(size=14,face="bold",hjust = 0.5),
        legend.text = element_text(size=12,face="bold")
      )
    # ggsave(paste0("SensitivityAnalysisFigs/PFSPoptreat",poptreat.i,".png"),width=5.5,height=5)
}


# Vary mut

mut.vec = unique(smmry.df$mut)

for (i in 1:length(mut.vec)) {
  
  mut.i = mut.vec[i]
  
  sub.df = smmry.df[smmry.df$poptreat==base.pars$poptreat&
                      smmry.df$mut==mut.i&
                      smmry.df$swtchdly==base.pars$swtchdly&
                      smmry.df$a_B==base.pars$a_B&
                      smmry.df$trnovr==base.pars$trnovr&
                      smmry.df$ngr==base.pars$ngr,]
  
  ggplot(sub.df,aes(x=as.factor(infeff),y=gd.relfitness,fill=p.erad))+theme_bw()+
    geom_tile()+
    ggtitle(paste("Mutation Rate:",mut.i,"/div"))+
    xlab("Gene Delivery Efficiency")+
    ylab("Gene Drive Fitness\n(Relative to Native Resistance)")+
    scale_fill_gradientn("Erad\nProb",colors=viridis(10),breaks=c(0,1),limits=c(0,1))+
    scale_x_discrete(breaks=c(0,.001,0.01,0.1),labels=parse(text=c("0","10^-3","10^-2","10^-1")))+
    # scale_y_continuous(breaks=c(-0.14,0),labels=c("",""))+
    theme(
      plot.title = element_text(size=22,hjust=0.5,face="bold"),
      axis.text = element_text(size=16,face="bold",color="black"),
      axis.title = element_text(size=18,face="bold"),
      legend.title = element_text(size=14,face="bold",hjust = 0.5),
      legend.text = element_text(size=12,face="bold")
    )
  # ggsave(paste0("SensitivityAnalysisFigs/EradicationMut",mut.i,".png"),width=5.5,height=5)

  ggplot(sub.df,aes(x=as.factor(infeff),y=gd.relfitness,fill=med.PFS))+theme_bw()+
    geom_tile()+
    ggtitle(paste("Mutation Rate:",mut.i,"/div"))+
    xlab("Gene Delivery Efficiency")+
    ylab("Gene Drive Fitness\n(Relative to Native Resistance)")+
    scale_fill_gradientn("Median\nPFS",colors=viridis(10),limits=c(min.PFS,max.PFS))+
    scale_x_discrete(breaks=c(0,.001,0.01,0.1),labels=parse(text=c("0","10^-3","10^-2","10^-1")))+
    # scale_y_continuous(breaks=c(-0.14,0),labels=c("",""))+
    theme(
      plot.title = element_text(size=22,hjust=0.5,face="bold"),
      axis.text = element_text(size=16,face="bold",color="black"),
      axis.title = element_text(size=18,face="bold"),
      legend.title = element_text(size=14,face="bold",hjust = 0.5),
      legend.text = element_text(size=12,face="bold")
    )
  # ggsave(paste0("SensitivityAnalysisFigs/PFSMut",mut.i,".png"),width=5.5,height=5)
}

# Vary turnover

trnovr.vec = unique(smmry.df$trnovr)

for (i in 1:length(trnovr.vec)) {
  
  trnovr.i = trnovr.vec[i]
  
  sub.df = smmry.df[smmry.df$poptreat==base.pars$poptreat&
                      smmry.df$mut==base.pars$mut&
                      smmry.df$swtchdly==base.pars$swtchdly&
                      smmry.df$a_B==base.pars$a_B&
                      smmry.df$trnovr==trnovr.i&
                      smmry.df$ngr==base.pars$ngr,]
  
  ggplot(sub.df,aes(x=as.factor(infeff),y=gd.relfitness,fill=p.erad))+theme_bw()+
    geom_tile()+
    ggtitle(paste("Turnover Rate:",trnovr.i))+
    xlab("Gene Delivery Efficiency")+
    ylab("Gene Drive Fitness\n(Relative to Native Resistance)")+
    scale_fill_gradientn("Erad\nProb",colors=viridis(10),breaks=c(0,1),limits=c(0,1))+
    scale_x_discrete(breaks=c(0,.001,0.01,0.1),labels=parse(text=c("0","10^-3","10^-2","10^-1")))+
    # scale_y_continuous(breaks=c(-0.14,0),labels=c("",""))+
    theme(
      plot.title = element_text(size=22,hjust=0.5,face="bold"),
      axis.text = element_text(size=16,face="bold",color="black"),
      axis.title = element_text(size=18,face="bold"),
      legend.title = element_text(size=14,face="bold",hjust = 0.5),
      legend.text = element_text(size=12,face="bold")
    )
  # ggsave(paste0("SensitivityAnalysisFigs/EradicationTurnover",trnovr.i,".png"),width=5.5,height=5)
  
  ggplot(sub.df,aes(x=as.factor(infeff),y=gd.relfitness,fill=med.PFS))+theme_bw()+
    geom_tile()+
    ggtitle(paste("Turnover Rate:",trnovr.i))+
    xlab("Gene Delivery Efficiency")+
    ylab("Gene Drive Fitness\n(Relative to Native Resistance)")+
    scale_fill_gradientn("Median\nPFS",colors=viridis(10),limits=c(min.PFS,max.PFS))+
    scale_x_discrete(breaks=c(0,.001,0.01,0.1),labels=parse(text=c("0","10^-3","10^-2","10^-1")))+
    # scale_y_continuous(breaks=c(-0.14,0),labels=c("",""))+
    theme(
      plot.title = element_text(size=22,hjust=0.5,face="bold"),
      axis.text = element_text(size=16,face="bold",color="black"),
      axis.title = element_text(size=18,face="bold"),
      legend.title = element_text(size=14,face="bold",hjust = 0.5),
      legend.text = element_text(size=12,face="bold")
    )
  # ggsave(paste0("SensitivityAnalysisFigs/PFSTurnover",trnovr.i,".png"),width=5.5,height=5)
}

# Vary net growth rate

ngr.vec = unique(smmry.df$ngr)

for (i in 1:length(ngr.vec)) {
  
  ngr.i = ngr.vec[i]
  
  sub.df = smmry.df[smmry.df$poptreat==base.pars$poptreat&
                      smmry.df$mut==base.pars$mut&
                      smmry.df$swtchdly==base.pars$swtchdly&
                      smmry.df$a_B==base.pars$a_B&
                      smmry.df$trnovr==base.pars$trnovr&
                      smmry.df$ngr==ngr.i,]
  
  ggplot(sub.df,aes(x=as.factor(infeff),y=gd.relfitness,fill=p.erad))+theme_bw()+
    geom_tile()+
    ggtitle(paste("Net Growth Rate:",ngr.i,"/day"))+
    xlab("Gene Delivery Efficiency")+
    ylab("Gene Drive Fitness\n(Relative to Native Resistance)")+
    scale_fill_gradientn("Erad\nProb",colors=viridis(10),breaks=c(0,1),limits=c(0,1))+
    scale_x_discrete(breaks=c(0,.001,0.01,0.1),labels=parse(text=c("0","10^-3","10^-2","10^-1")))+
    # scale_y_continuous(breaks=c(-0.14,0),labels=c("",""))+
    theme(
      plot.title = element_text(size=22,hjust=0.5,face="bold"),
      axis.text = element_text(size=16,face="bold",color="black"),
      axis.title = element_text(size=18,face="bold"),
      legend.title = element_text(size=14,face="bold",hjust = 0.5),
      legend.text = element_text(size=12,face="bold")
    )
  # ggsave(paste0("SensitivityAnalysisFigs/EradicationNetgrowth",ngr.i,".png"),width=5.5,height=5)
  
  ggplot(sub.df,aes(x=as.factor(infeff),y=gd.relfitness,fill=med.PFS))+theme_bw()+
    geom_tile()+
    ggtitle(paste("Turnover Rate:",ngr.i,"/day"))+
    xlab("Gene Delivery Efficiency")+
    ylab("Gene Drive Fitness\n(Relative to Native Resistance)")+
    scale_fill_gradientn("Median\nPFS",colors=viridis(10),limits=c(min.PFS,max.PFS))+
    scale_x_discrete(breaks=c(0,.001,0.01,0.1),labels=parse(text=c("0","10^-3","10^-2","10^-1")))+
    # scale_y_continuous(breaks=c(-0.14,0),labels=c("",""))+
    theme(
      plot.title = element_text(size=22,hjust=0.5,face="bold"),
      axis.text = element_text(size=16,face="bold",color="black"),
      axis.title = element_text(size=18,face="bold"),
      legend.title = element_text(size=14,face="bold",hjust = 0.5),
      legend.text = element_text(size=12,face="bold")
    )
  # ggsave(paste0("SensitivityAnalysisFigs/PFSNetgrowth",ngr.i,".png"),width=5.5,height=5)
}

# Vary switch delay

swtchdly.vec = unique(smmry.df$swtchdly)

sub.df = smmry.df[smmry.df$poptreat==base.pars$poptreat&
                    smmry.df$mut==base.pars$mut&
                    smmry.df$infeff==0.01&
                    smmry.df$dimres==0&
                    smmry.df$a_B==base.pars$a_B&
                    smmry.df$trnovr==base.pars$trnovr&
                    smmry.df$ngr==base.pars$ngr,]

sub.df$swtchdly.mnths = round(sub.df$swtchdly/30.4167)

sub.df = sub.df[sub.df$swtchdly.mnths<=40,]

# Draw lines between points
line.plot = approx(sub.df$swtchdly.mnths,sub.df$p.erad,n=1000)
line.plot.df = data.frame(x=line.plot$x,y=line.plot$y)

ggplot(sub.df,aes(x=swtchdly.mnths,y=p.erad*100))+theme_bw()+
  geom_point(data=line.plot.df,aes(x=x,y=y*100,color=y),size=1.25)+
  geom_point(aes(fill=p.erad),shape=21,color="black",size=4)+
  xlab("Switch Delay (months)")+
  ylab("Eradication Probability (%)")+
  ggtitle("Optimizing Switch Scheduling")+
  scale_color_gradientn("Erad\nProb",colors=viridis(10),breaks=c(0,1),limits=c(0,1))+
  scale_fill_gradientn("Erad\nProb",colors=viridis(10),breaks=c(0,1),limits=c(0,1))+
  theme(
    plot.title = element_text(size=22,hjust=0.5,face="bold"),
    axis.text = element_text(size=16,face="bold",color="black"),
    axis.title = element_text(size=18,face="bold"),
    legend.title = element_text(size=14,face="bold",hjust = 0.5),
    legend.text = element_text(size=12,face="bold")
  )
# ggsave("SwitchScheduling.png",width=5.5,height=5)
