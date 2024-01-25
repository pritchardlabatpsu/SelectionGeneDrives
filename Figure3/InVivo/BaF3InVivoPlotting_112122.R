### Plot tumor volumes for PC9C2 in vivo gene drive experiment

## Set up
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(dplyr)
library(ggbreak)

## Import and parse data
msrmnts.df = read.csv("InVivoBaF3MeasurementsRaw_112022.csv",header=T,stringsAsFactors=F)
msrmnts.df$mouse = as.factor(msrmnts.df$mouse)
msrmnts.df$mouse.flank = paste(msrmnts.df$mouse,msrmnts.df$flank)

## Intermediate calculations

msrmnts.df$volume1 = pi*msrmnts.df$length1*msrmnts.df$width1*msrmnts.df$height1/6
msrmnts.df$volume2 = pi*msrmnts.df$length2*msrmnts.df$width2*msrmnts.df$height2/6
msrmnts.df$volume3 = pi*msrmnts.df$length3*msrmnts.df$width3*msrmnts.df$height3/6
msrmnts.df$volume = rowSums(msrmnts.df[,c("volume1","volume2","volume3")],na.rm=T)

# For missing counts, replace zeros from rowSums with NAs
msrmnts.df$volume[is.na(msrmnts.df$length1)] = NA

# Find relative change in volume
msrmnts.df$rel.volume = NA
mouse.flanks = unique(msrmnts.df$mouse.flank)
for (mouse.flank in mouse.flanks) {
  
  init.vol = msrmnts.df$volume[msrmnts.df$day==0&msrmnts.df$mouse.flank==mouse.flank]
  
  msrmnts.df$rel.volume[msrmnts.df$mouse.flank==mouse.flank] = msrmnts.df$volume[msrmnts.df$mouse.flank==mouse.flank]/init.vol
  
}

## Plot data

ggplot(msrmnts.df[msrmnts.df$condition==1,],aes(x=day,y=volume,color=mouse.flank))+theme_bw()+
  geom_point()+
  geom_line()

## Summarize data

# msrmnts.df = msrmnts.df[msrmnts.df$mouse.flank!="3 right",]

msrmnts.agg = msrmnts.df %>%
  group_by(condition,day) %>%
  summarize(vol.mean = mean(volume,na.rm=T),
            vol.sem = sd(volume,na.rm=T)/sqrt(length(volume)),
            rel.vol.mean = mean(rel.volume,na.rm=T),
            rel.vol.sem = sd(rel.volume,na.rm=T)/sqrt(length(rel.volume)))

# Omit condition once we began sacking mice (day 17 for condition 1; day 23 for condition 2)

msrmnts.agg = msrmnts.agg[(msrmnts.agg$condition==1&msrmnts.agg$day<=17)|
                            (msrmnts.agg$condition==2)|
                            (msrmnts.agg$condition==3&msrmnts.agg$day<=21),]
msrmnts.agg$condition = as.factor(msrmnts.agg$condition)

# Absolute volume
ggplot(msrmnts.agg,aes(x=day,y=vol.mean,color=condition))+theme_bw()+
  geom_point(size=3)+
  geom_line(size=1.5)+
  geom_errorbar(aes(ymin=vol.mean-vol.sem,ymax=vol.mean+vol.sem),width=0.5,size=1.5)

# Change in volume
ggplot(msrmnts.agg,aes(x=day,y=rel.vol.mean,color=condition))+theme_bw()+
  geom_point(size=3)+
  geom_line(size=1.5)+
  geom_errorbar(aes(ymin=rel.vol.mean-rel.vol.sem,ymax=rel.vol.mean+rel.vol.sem),width=0.5,size=1.5)

## Plot change in volume by condition
# Not doing in a for loop because scale_x_break has a hard time handling ambiguity

# Condition 1
cndtn.i = 1
ggplot(msrmnts.agg[msrmnts.agg$condition==cndtn.i,],
       aes(x=day,y=rel.vol.mean))+theme_bw()+
  geom_line(size=2)+
  geom_point(size=4)+
  geom_errorbar(aes(ymin=rel.vol.mean-rel.vol.sem,ymax=rel.vol.mean+rel.vol.sem),width=0.5,size=1.25)+
  ggtitle(paste("Condition",cndtn.i))+
  xlab("Days")+
  ylab("Change in Tumor Volume")+
  ylim(0,1.75)+
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size=18,face="bold"),
    axis.text = element_text(size=16,color="black",face="bold"),
    axis.text.x.top = element_blank(),
    axis.ticks.x.top = element_blank()
  )
# ggsave(paste0("InVivo_Condition",cndtn.i,"_112122.png"),width=5,height=5)

# Condition 2
cndtn.i = 2
ggplot(msrmnts.agg[msrmnts.agg$condition==cndtn.i,],
       aes(x=day,y=rel.vol.mean))+theme_bw()+
  geom_line(size=2)+
  geom_point(size=4)+
  geom_errorbar(aes(ymin=rel.vol.mean-rel.vol.sem,ymax=rel.vol.mean+rel.vol.sem),width=0.5,size=1.25)+
  ggtitle(paste("Condition",cndtn.i))+
  xlab("Days")+
  ylab("Change in Tumor Volume")+
  ylim(0,1.75)+
  scale_x_break(c(17,39))+
  scale_x_continuous(breaks=c(0,5,10,15,40),limits = c(0,41))+
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size=18,face="bold"),
    axis.text = element_text(size=16,color="black",face="bold"),
    axis.text.x.top = element_blank(),
    axis.ticks.x.top = element_blank()
  )
# ggsave(paste0("InVivo_Condition",cndtn.i,"_112122.png"),width=5,height=5)


# Condition 3
cndtn.i = 3
ggplot(msrmnts.agg[msrmnts.agg$condition==cndtn.i,],
       aes(x=day,y=rel.vol.mean))+theme_bw()+
  geom_line(size=2)+
  geom_point(size=4)+
  geom_errorbar(aes(ymin=rel.vol.mean-rel.vol.sem,ymax=rel.vol.mean+rel.vol.sem),width=0.5,size=1.25)+
  ggtitle(paste("Condition",cndtn.i))+
  xlab("Days")+
  ylab("Change in Tumor Volume")+
  ylim(0,1.75)+
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size=18,face="bold"),
    axis.text = element_text(size=16,color="black",face="bold"),
    axis.text.x.top = element_blank(),
    axis.ticks.x.top = element_blank()
  )
# ggsave(paste0("InVivo_Condition",cndtn.i,"_112122.png"),width=5,height=5)


# Change in volume
ggplot(msrmnts.agg[msrmnts.agg$condition%in%c(1,2),],aes(x=day,y=rel.vol.mean,color=condition))+theme_bw()+
  geom_point(size=3)+
  geom_line(size=1.5)+
  geom_errorbar(aes(ymin=rel.vol.mean-rel.vol.sem,ymax=rel.vol.mean+rel.vol.sem),width=0.5,size=1.5)+
  xlab("Time (days)")+
  ylab("Change in Tumor Volume")+
  scale_x_break(c(17,39))+
  scale_x_continuous(breaks=c(0,5,10,15,40),limits = c(0,41))+
  scale_color_manual(values=c("#EE1D23","#1AD31A"))+
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size=18,face="bold"),
    axis.text = element_text(size=16,color="black",face="bold"),
    axis.text.x.top = element_blank(),
    axis.ticks.x.top = element_blank()
  )

## Determine when to sack mice

# Sack when one tumor exceeds 150 or total exceeds 250
single.flank.max = 150
total.max = 250

last.day = max(msrmnts.df$day)
msrmnts.last = msrmnts.df[msrmnts.df$day==last.day,]

total.burden = aggregate(msrmnts.df$volume,
                         by = list(mouse = msrmnts.df$mouse,
                                   day = msrmnts.df$day),
                         FUN = sum)
colnames(total.burden)[3] = "total.volume"
total.burden$total.volume.exceeded = total.burden$total.volume>total.max

total.burden$single.flank.exceeded = F
for (i in 1:nrow(total.burden)) {
  mouse.i = total.burden$mouse[i]
  day.i = total.burden$day[i]
  flank.volumes = msrmnts.df$volume[msrmnts.df$mouse==mouse.i&msrmnts.df$day==day.i]
  total.burden$single.flank.exceeded[i] = sum(flank.volumes>single.flank.max)>0
}

total.burden$sack = total.burden$total.volume.exceeded|total.burden$single.flank.exceeded
total.burden$sack[is.na(total.burden$sack)] = F
total.burden$mouse[total.burden$sack]

## Plot Kaplan-Meier curve

treat.failure.df = data.frame(
  mouse = 1:15,
  condition = c(rep(1,5),rep(2,5),rep(3,5)),
  day.fail = rep(NA,15)
)

for (mouse in 1:15) {
  
  sub.burden = total.burden[total.burden$mouse==mouse,]
  day.fail = sub.burden$day[which(sub.burden$sack)[1]]
  treat.failure.df$day.fail[treat.failure.df$mouse==mouse] = day.fail
  
}

treat.failure.df$day.fail[treat.failure.df$mouse==13] = 26 # tumors didn't quite reach threshold for sacking, but mouse was clearly sick on D26; after sacking, quick dissection showed very enlarged spleen, possibly from metastases

KM.df = data.frame()

for (condition in 1:3) {
  
  sub.treat.fail = treat.failure.df[treat.failure.df$condition==condition,]
  day.fails = sub.treat.fail$day.fail
  
  # Order failure days
  day.fails = day.fails[order(day.fails)]
  
  sub.KM = data.frame(condition=condition,day=0,frac.surv=1)
  
  n.mice = length(day.fails)
  frac.surv = 1
  for (i in 1:n.mice) {
    sub.KM = rbind(sub.KM,
                   data.frame(
                     condition=rep(condition,2),
                     day=rep(day.fails[i],2),
                     frac.surv=c(frac.surv,frac.surv-1/n.mice))
                   )
    frac.surv = frac.surv-1/n.mice
  }
  sub.KM$frac.surv = round(sub.KM$frac.surv,digits=1)
  
  KM.df = rbind(KM.df,sub.KM)
}

KM.df$condition = as.factor(KM.df$condition)

# Include point for condition 2
KM.df = rbind(KM.df,
              data.frame(condition = 2,
                         day = 45,
                         frac.surv = 1)
              )

# Plot

ggplot(KM.df,aes(x=day,y=frac.surv,color=condition))+theme_bw()+
  geom_line()
  