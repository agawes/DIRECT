module load gcc/6.2.0
module load intel/redist/2017.2.174
module load intel/perflibs/64/2017_update2
module load R/3.4.0-ICC-MKL

R

library(ggplot2)
library(data.table)

options(width=150)
setwd("~/Archetypes_Age/")
archetypes=read.table("archetypes.coefficients.wp2.2.txt",h=T,sep="\t")

### read in psGRSs - generated by Anubha
psGRS=read.table("/home/Teams/teamVIP/Analysis/Share_Results/Scores.for.Agata.30thMay2019AM.txt",h=T)
names(psGRS)[1]="studyid"
psGRS=merge(archetypes, psGRS, by="studyid")

GRS_by_archetype=data.frame()
for (a in 1:4){
  df=data.frame(archetype=rep(colnames(psGRS)[a+1],7), GRS=colnames(psGRS[6:12]), 
                p=sapply(1:7, function(x) summary(lm(psGRS[,1+a] ~ psGRS[,5+x]))$coefficients[8]),
                beta=sapply(1:7, function(x) summary(lm(psGRS[,1+a] ~ psGRS[,5+x]))$coefficients[2]) )
  GRS_by_archetype=rbind(GRS_by_archetype, df)
}

GRS_by_archetype$logP=-log10(GRS_by_archetype$p) * GRS_by_archetype$beta/abs(GRS_by_archetype$beta)
GRS_by_archetype$star=cut(GRS_by_archetype$p, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", "")) 

pdf("Fig3A.GRS_by_archetype.regression.pdf",width=6,height=4)
p=ggplot(aes(x=GRS, y=archetype, fill=logP), data=GRS_by_archetype)
p + geom_tile() + scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C") + 
  geom_text(aes(label=star), color="black", size=3) + scale_y_discrete(limits = rev(levels(GRS_by_archetype$archetype))) +
  labs(y="Archetype", x=NULL, fill="-log10P") + theme_bw()
dev.off()
