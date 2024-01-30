
### Analyze results from gene drive spatial agent-based model ###

## Set up

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(viridis)

## Loop through results

csv.files = list.files(pattern="*.csv$")

sims.df = data.frame()

for (i in 1:length(csv.files)) {
  
  file.i = csv.files[i]
  
  # Get simulation attributes
  file.splt = strsplit(file.i,"_")[[1]]
  size0 = as.numeric(gsub("size0","",file.splt[2])) # initial population size
  tight = as.numeric(gsub("tight","",file.splt[3])) # measure of dispersion of gene drive cells (0 indicates even spatial distribution, large numbers indicate tight spatial gathering into foci)
  nfoci = as.numeric(gsub("nfoci","",file.splt[4])) # number of foci of gene drive cells
  killrad = as.numeric(gsub("killrad","",file.splt[5])) # radius of bystander effect (in cell lengths)
  swtchdly = as.numeric(gsub("swtchdly","",file.splt[6])) # delay between turning on switch 2 and turning off switch 1
  sim = as.numeric(gsub("iter","",file.splt[7])) # simulation number
  
  # Read in file
  df.i = read.csv(file.i,header=F)
  colnames(df.i) = c("iter","time","nS","nR","nG","nAS","nAR")
  # nS: sensitive
  # nR: resistant
  # nG: gene drive
  # nAS: sensitive, within kill radius of gene drive
  # nAR: resistant, within kill radius of gene drive
  
  # Check if eradication
  final.pop = df.i[nrow(df.i),3:7]
  erad = sum(final.pop)==0 # note: some simulations reach time limit before all gene drive cells are cleared (but no other cells remain) - note as eradicated
  
  sims.i = data.frame(
    size0 = size0,
    tight = tight,
    nfoci = nfoci,
    killrad = killrad,
    swtchdly = swtchdly,
    sim = sim,
    erad = erad
  )
  
  sims.df = rbind(sims.df,sims.i)
  
}

## Summarize results
sims.agg = aggregate(sims.df$erad,
                    by = list(size0 = sims.df$size0,
                              tight = sims.df$tight,
                              nfoci = sims.df$nfoci,
                              killrad = sims.df$killrad,
                              swtchdly = sims.df$swtchdly),
                    FUN = mean)
colnames(sims.agg)[6] = "p.erad"

## Replace "tight" metric with measure of average distance between gene drive cells
# Values obtained from GeneDriveDispersionMeasurements
disp.df = read.csv("../GeneDriveDispersionMeasurements/GeneDriveDispersionMetric_091422.csv",header=F)
colnames(disp.df) = c("size0","tight","sim","mean.dist")

disp.agg = aggregate(disp.df$mean.dist,
                     by=list(
                       size0 = disp.df$size0,
                       tight = disp.df$tight
                     ),
                     FUN = mean)
colnames(disp.agg) = c("size0","tight","mean.mean.dist")

# Calculate mean.mean.dist relative to minimum and maximum possible values
disp.agg$rel.mean.dist = NA
sims.agg$mean.dist = NA
sims.df$mean.dist = NA
sims.agg$rel.mean.dist = NA
sims.df$rel.mean.dist = NA

for (i in 1:nrow(disp.agg)) {
  
  size0.i = disp.agg$size0[i]
  tight.i = disp.agg$tight[i]
  mean.mean.dist.i = disp.agg$mean.mean.dist[i]
  
  sims.agg$mean.dist[sims.agg$size0==size0.i&sims.agg$tight==tight.i] = mean.mean.dist.i
  sims.df$mean.dist[sims.df$size0==size0.i&sims.df$tight==tight.i] = mean.mean.dist.i
  
  # Calculate mean.mean.dist relative to minimum and maximum possible values
  max.mean.dist.i = max(disp.agg$mean.mean.dist[disp.agg$size0==size0.i])
  min.mean.dist.i = min(disp.agg$mean.mean.dist[disp.agg$size0==size0.i])
  rel.mean.dist.i = ((mean.mean.dist.i-min.mean.dist.i)/(max.mean.dist.i-min.mean.dist.i))
  sims.agg$rel.mean.dist[sims.agg$size0==size0.i&sims.agg$tight==tight.i] = rel.mean.dist.i
  disp.agg$rel.mean.dist[disp.agg$size0==size0.i&disp.agg$tight==tight.i] = rel.mean.dist.i
  sims.df$rel.mean.dist[sims.df$size0==size0.i&sims.df$tight==tight.i] = rel.mean.dist.i
  
}

## Visualize results

ggplot(sims.agg[sims.agg$size0==10000&sims.agg$killrad!=1.5,],aes(x=as.factor(killrad),y=as.factor(round(rel.mean.dist,2)),fill=p.erad))+theme_bw()+
  geom_raster()+
  scale_fill_gradientn("",colors=viridis(10),breaks=c(0,1),limits=c(0,1))+
  xlab("Bystander Activity Distance\n(Cell Lengths)")+
  # ylab("Initial Gene Drive Cell\nDispersion (Cell Lengths)")+
  ylab("Relative Dispersion")+
  # scale_y_discrete(labels=c("3.6","4.9","6.1","7.5","8.0"))+
  ggtitle("Spatial ABM Results")+
  theme(
    plot.title = element_text(size=22,hjust=0.5,face="bold"),
    axis.text = element_text(size=16,face="bold",color="black"),
    axis.title = element_text(size=18,face="bold"),
    legend.text = element_blank()
  )

# ggsave("GeneDriveABMResults.png",width=5.5,height=5)

# Check relationship between tightness and probability of eradication

killrad.vec = unique(sims.df$killrad)

p.vals = rep(NA,length(rel.mean.dist.vec))

for (i in 1:length(killrad.vec)) {
  
  killrad.i = killrad.vec[i]
  
  x = sims.df$tight[sims.df$killrad==killrad.i]
  y = sims.df$erad[sims.df$killrad==killrad.i]
  fit = lm(y~x)
  smmry = summary(fit)
  p.vals[i] = smmry$coefficients[2,4]
  
}

plot.pval.df = data.frame(killrad = killrad.vec,pval = p.vals)
plot.pval.df$significant = plot.pval.df$pval<=0.05

ggplot(plot.pval.df,aes(x=killrad,y=pval,color=significant))+theme_bw()+
  geom_point(size=3)+
  geom_hline(yintercept=0.05,linetype='dashed',color='gray',size=1)+
  scale_color_manual(values=c("gray15","red"))+
  xlab("Bystander Activity Distance")+
  ylab("P-value")+
  ggtitle("Linear Regression Summary:\nEradication Probability vs.\nDispersion Parameter")+
  theme(
    plot.title = element_text(size=22,hjust=0.5,face="bold"),
    axis.text = element_text(size=16,face="bold",color="black"),
    axis.title = element_text(size=18,face="bold")
  )+
  guides(color="none")

ggsave("RegressionSummary.png",width=5,height=5)

