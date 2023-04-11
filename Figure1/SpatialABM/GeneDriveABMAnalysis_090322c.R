
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
disp.df = read.csv("../GeneDriveDispersionMeasurements/GeneDriveDispersionMetric.csv",header=F)
colnames(disp.df) = c("size0","tight","sim","mean.dist")

disp.agg = aggregate(disp.df$mean.dist,
                     by=list(
                       size0 = disp.df$size0,
                       tight = disp.df$tight
                     ),
                     FUN = mean)
colnames(disp.agg) = c("size0","tight","mean.mean.dist")

sims.agg$mean.dist = NA
for (i in 1:nrow(disp.agg)) {
  
  size0.i = disp.agg$size0[i]
  tight.i = disp.agg$tight[i]
  mean.mean.dist.i = disp.agg$mean.mean.dist[i]
  
  sims.agg$mean.dist[sims.agg$size0==size0.i&sims.agg$tight==tight.i] = round(mean.mean.dist.i,1)
  
}

## Visualize results

ggplot(sims.agg[sims.agg$size0==5000&sims.agg$nfoci==1&sims.agg$swtchdly==180&!sims.agg$killrad%in%c(1,1.75),],aes(x=as.factor(killrad),y=rev(as.factor(mean.dist)),fill=p.erad))+theme_bw()+
  geom_raster()+
  scale_fill_gradientn("",colors=viridis(10),breaks=c(0,0.5,1),limits=c(0,1))+
  xlab("Activated Prodrug Diffusion\nDistance (Cell Lengths)")+
  ylab("Initial Gene Drive Cell\nDispersion (Cell Lengths)")+
  scale_y_discrete(labels=c("3.6","4.9","6.1","7.5","8.0"))+
  ggtitle("Spatial ABM Results")+
  theme(
    plot.title = element_text(size=22,hjust=0.5,face="bold"),
    axis.text = element_text(size=16,face="bold",color="black"),
    axis.title = element_text(size=18,face="bold")
  )

# ggsave("GeneDriveABMResults.png",width=5.5,height=5)

