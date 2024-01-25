
## Set up

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(reshape2)
library(ggbreak)

## Import and parse data

# Import
counts.df = read.csv("BrunelloGeneDriveCounts.csv",header=T,stringsAsFactors=F)

counts.df$day = as.numeric(gsub("D","",counts.df$day))
counts.df$gene.drive.status = gsub("'","",counts.df$gene.drive.status)
counts.df$replicate = gsub("rep","",counts.df$replicate)

counts.df$gd = as.numeric(gsub(",","",counts.df$gd))
counts.df$Brnlo = as.numeric(gsub(",","",counts.df$Brnlo))

counts.df$gd.conc = counts.df$gd/counts.df$volume
counts.df$Brnlo.conc = counts.df$Brnlo/counts.df$volume

# Reshape

counts.rs = melt(counts.df[,c(1:4,8,9)],id.vars=c("osi.conc","gene.drive.status","replicate","day"))

## Loop through conditions and plot

osi.concs = c("30nM","50nM")
GD.status = c("-GD","+GD")

cndtns = expand.grid(osi.conc = osi.concs,GD.status = GD.status)

# Scaling factor
# .conc values are k/mL and samples are aliquots from 12 mL
# so scaling to M cells in 12 mL
sclng.fctr = 12/1000

for (i in 1:nrow(cndtns)) {
  
  osi.conc.i = cndtns$osi.conc[i]
  GD.stat.i = cndtns$GD.status[i]
  
  sub.counts = counts.rs[counts.rs$osi.conc==osi.conc.i&counts.df$gene.drive.status==GD.stat.i,]
  
  # Aggregate counts
  sub.agg = aggregate(sub.counts$value,
                      by = list(osi.conc=sub.counts$osi.conc,
                                GD.stat=sub.counts$gene.drive.status,
                                day=sub.counts$day,
                                population = sub.counts$variable),
                      FUN = function(x) c(mean=mean(x),sem=sd(x)/sqrt(length(x))))
  sub.agg = do.call(data.frame,sub.agg)
  
  # Remove GFP counts for -GD
  if (GD.stat.i=="-GD") {
    sub.agg = sub.agg[sub.agg$population!="gd.conc",]
    clrs = c("#EE1D23")
    # Plot
    # g = ggplot(sub.agg,aes(x=day,y=x.mean*sclng.fctr,color=population))+theme_bw()+
    #   geom_errorbar(aes(ymin=(x.mean-x.sem)*sclng.fctr,ymax=(x.mean+x.sem)*sclng.fctr),width=0.5,size=1.25,color="black")+
    #   geom_line(size=2)+
    #   geom_point(size=5)+
    #   geom_point(data=sub.agg[sub.agg$population=="gd.conc",],size=4)+ # hack to get the green points over the red ones
    #   scale_color_manual(values=clrs)+
    #   scale_y_continuous(breaks=c(0,2,4,6),limits=c(0,6))+
    #   scale_x_continuous(breaks=c(0,7,14,21),limits=c(-0.5,28))+
    #   xlab("Days")+ylab("Population Size (M)")+
    #   ggtitle(paste(osi.conc.i,GD.stat.i))+
    #  theme(
    #   plot.title = element_text(size=22,face="bold",hjust=0.5),
    #   axis.title = element_text(size=20,face="bold"),
    #   axis.text = element_text(size=18,face="bold",color="black"),
    #   legend.title = element_text(size=18,face="bold"),
    #   legend.text = element_text(size=16,face="bold")
    # )+
    # guides(color="none")
    
  } else {
    clrs = c("#EE1D23","#1AD31A")
    sub.agg$population = factor(sub.agg$population,levels=c("Brnlo.conc","gd.conc")) # So gene drive is on top in plot
    
    # Plot
    # g = ggplot(sub.agg,aes(x=day,y=x.mean*sclng.fctr,color=population))+theme_bw()+
    #   geom_errorbar(aes(ymin=(x.mean-x.sem)*sclng.fctr,ymax=(x.mean+x.sem)*sclng.fctr),width=0.5,size=1.25,color="black")+
    #   geom_line(size=2)+
    #   geom_point(size=4)+
    #   geom_point(data=sub.agg[sub.agg$population=="gd.conc",],size=4)+ # hack to get the green points over the red ones
    #   scale_color_manual(values=clrs)+
    #   scale_y_continuous(breaks=c(0,2,4,6),limits=c(0,6))+
    #   scale_x_break(c(18,29))+
    #   scale_x_continuous(breaks=c(0,7,14,21,30),limits=c(-0.5,31))+
    #   xlab("Days")+ylab("Population Size (M)")+
    #   ggtitle(paste(osi.conc.i,GD.stat.i))+
    #  theme(
    #   plot.title = element_text(size=22,face="bold",hjust=0.5),
    #   axis.title = element_text(size=20,face="bold"),
    #   axis.text = element_text(size=18,face="bold",color="black"),
    #   axis.text.y.right = element_blank(),
    #   axis.text.x.top = element_blank(),
    #   axis.ticks.x.top = element_blank(),
    #   legend.title = element_text(size=18,face="bold"),
    #   legend.text = element_text(size=16,face="bold")
    # )+
    # guides(color="none")
  }
  
  g = ggplot(sub.agg,aes(x=day,y=x.mean*sclng.fctr,color=population))+theme_bw()+
    geom_errorbar(aes(ymin=(x.mean-x.sem)*sclng.fctr,ymax=(x.mean+x.sem)*sclng.fctr),width=0.5,size=1.25)+
    geom_line(size=2)+
    geom_point(size=4)+
    geom_point(data=sub.agg[sub.agg$population=="gd.conc",],size=4)+ # hack to get the green points over the red ones
    scale_color_manual(values=clrs)+
    scale_y_continuous(breaks=c(0,2,4,6),limits=c(0,6))+
    scale_x_continuous(breaks=c(0,9,18,27),limits=c(-1,28))+ # Note here: rather than plotting the day 40 time point here with a break in the x-axis, we're just going to relabel the D28 in illustrator as D40 and introduce the axis break there, since the D28 and D40 points are both zeros
    xlab("Days")+ylab("Population Size (M)")+
    ggtitle(paste(osi.conc.i,GD.stat.i))+
    theme(
      plot.title = element_text(size=22,face="bold",hjust=0.5),
      axis.title = element_text(size=20,face="bold"),
      axis.text = element_text(size=18,face="bold",color="black"),
      legend.title = element_text(size=18,face="bold"),
      legend.text = element_text(size=16,face="bold")
    )+
    guides(color="none")
  
  print(g)
  
  # ggsave(paste0("PC9BrunelloGD",osi.conc.i,GD.stat.i,".png"),width=5,height=5)
  
  
}

## Alt: plot individually because ggbreak doesn't handle for-loops well

# Aggregate

agg.df = aggregate(counts.rs$value,
                    by = list(osi.conc=counts.rs$osi.conc,
                              GD.stat=counts.rs$gene.drive.status,
                              day=counts.rs$day,
                              population = counts.rs$variable),
                    FUN = function(x) c(mean=mean(x),sem=sd(x)/sqrt(length(x))))
agg.df = do.call(data.frame,agg.df)

# Without gene drive

ggplot(agg.df[agg.df$osi.conc=="30nM"&agg.df$GD.stat=="-GD"&agg.df$population=="Brnlo.conc",], # don't plot 0 counts for gene drive
       aes(x=day,y=x.mean*sclng.fctr,color=population))+theme_bw()+
  geom_errorbar(aes(ymin=(x.mean-x.sem)*sclng.fctr,ymax=(x.mean+x.sem)*sclng.fctr),width=0.5,size=1.25)+
  geom_line(size=2)+
  geom_point(size=4)+
  scale_color_manual(values=clrs)+
  scale_y_continuous(breaks=c(0,2,4,6),limits=c(0,6))+
  scale_x_continuous(breaks=c(0,12,24),limits=c(-1,28))+
  xlab("Days")+ylab("Population Size (M)")+
  ggtitle("Without Gene Drive")+
  theme(
    plot.title = element_text(size=22,face="bold",hjust=0.5),
    axis.title = element_text(size=20,face="bold"),
    axis.text = element_text(size=18,face="bold",color="black"),
    legend.title = element_text(size=18,face="bold"),
    legend.text = element_text(size=16,face="bold")
  )+
  guides(color="none")

# ggsave("BrunelloGrowthTracking_WithoutGeneDrive.png",width=5,height=5)

ggplot(agg.df[agg.df$osi.conc=="50nM"&agg.df$GD.stat=="+GD",], # don't plot 0 counts for gene drive
       aes(x=day,y=x.mean*sclng.fctr,color=population))+theme_bw()+
  geom_errorbar(aes(ymin=(x.mean-x.sem)*sclng.fctr,ymax=(x.mean+x.sem)*sclng.fctr),width=0.5,size=1.25)+
  geom_line(size=2)+
  geom_point(size=4)+
  geom_point(data=agg.df[agg.df$osi.conc=="50nM"&agg.df$GD.stat=="+GD"&agg.df$population=="gd.conc",],size=4)+ # hack to get the green points over the red ones
  geom_line(data=agg.df[agg.df$osi.conc=="50nM"&agg.df$GD.stat=="+GD"&agg.df$population=="gd.conc",],size=2)+ # hack to get the green points over the red ones
  scale_color_manual(values=c("#1AD31A","#EE1D23"))+
  scale_y_continuous(breaks=c(0,2,4,6),limits=c(0,6))+
  scale_x_continuous(breaks=c(0,12,24,27),limits=c(-1,28),labels=c(0,12,24,36))+ # Note here: rather than plotting the day 36 time point here with a break in the x-axis, we're just going to relabel the D28 in illustrator as D40 and introduce the axis break there, since the D28 and D40 points are both zeros
  xlab("Days")+ylab("Population Size (M)")+
  ggtitle("With Gene Drive")+
  theme(
    plot.title = element_text(size=22,face="bold",hjust=0.5),
    axis.title = element_text(size=20,face="bold"),
    axis.text = element_text(size=18,face="bold",color="black"),
    legend.title = element_text(size=18,face="bold"),
    legend.text = element_text(size=16,face="bold")
  )+
  guides(color="none")

# ggsave("BrunelloGrowthTracking_WithGeneDrive.png",width=5,height=5)
