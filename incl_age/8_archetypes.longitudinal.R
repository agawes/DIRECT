module load gcc/6.2.0
module load intel/redist/2017.2.174
module load intel/perflibs/64/2017_update2
module load R/3.4.0-ICC-MKL

R

library(archetypes)
library(RColorBrewer)
library(ggfortify) # for the PCA plot
library(fossil)

setwd("~/Archetypes_Age/longitudinal/")
phe=read.table("~/Archetypes_Age/WP2_2.m0_18_36.rnt_resid_flt.011019.txt",h=T,sep="\t")
phe$m=gsub(".+_","",rownames(phe))
phe$m=factor(phe$m)
table(phe$m)
# M0 M18 M36 
# 726 591 423 
table(table(phe$studyid))
#   1   2   3 
# 155 263 353 

#### run the archetypes analysis on the full longitudinal dataset with 4 archetypes
m=phe[,c(4:35)]

seed=55555

nrep=100
aa=list()
set.seed(seed)
aa<-stepArchetypes(m, k=4, nrep=nrep, method=robustArchetypes)

### get the best model 
best=bestModel(aa)

### coefficients (soft clustering)
coef=coef(best)

### hard clusters:
clusters=max.col(coef)

### explore the correlations of phenotypes with archetypes:
data.frame(sort(cor(cbind(coef, m))[1,])) ## MOD
data.frame(sort(cor(cbind(coef, m))[2,]))  ## MARD
data.frame(sort(cor(cbind(coef, m))[3,]))  ## SIDD
data.frame(sort(cor(cbind(coef, m))[4,]))	## SIRD

## rename to compare with baseline clustering:

### reorder the coef matrix:
coef=coef[,c(2,1,4,3)]
clusters_copy=rep(0,length(clusters))
clusters_copy[which(clusters==1)]<-2
clusters_copy[which(clusters==2)]<-1
clusters_copy[which(clusters==3)]<-4
clusters_copy[which(clusters==4)]<-3
clusters=clusters_copy

### test at cut-off of 0.5 how many subjects are in each archetype
##### test how many subjects are still defined by the same major archetype
cut=0.6
longit_df=data.frame(studyid=as.character(unique(phe$studyid)) )

select=which(phe$m=="M0")
arch_M0=data.frame(coef[select,])
arch_M0$studyid=phe[select,]$studyid
arch_M0$arch_M0=sapply(1:nrow(arch_M0), function(x) 
	ifelse(arch_M0[x,(clusters[select][x])]>cut, clusters[select][x], "mix") )

colnames(arch_M0)[1:4]=c("A_M0","B_M0","C_M0","D_M0")
longit_df=merge(longit_df,arch_M0, by="studyid", all=T)

select=which(phe$m=="M18")
arch_M18=data.frame(coef[select,])
arch_M18$studyid=phe[select,]$studyid
arch_M18$arch_M18=sapply(1:nrow(arch_M18), function(x) 
	ifelse(arch_M18[x,(clusters[select][x])]>cut, clusters[select][x], "mix") )
colnames(arch_M18)[1:4]=c("A_M18","B_M18","C_M18","D_M18")
longit_df=merge(longit_df,arch_M18, by="studyid", all=T)

select=which(phe$m=="M36")
arch_M36=data.frame(coef[select,])
arch_M36$studyid=phe[select,]$studyid
arch_M36$arch_M36=sapply(1:nrow(arch_M36), function(x) 
	ifelse(arch_M36[x,(clusters[select][x])]>cut, clusters[select][x], "mix") )
colnames(arch_M36)[1:4]=c("A_M36","B_M36","C_M36","D_M36")
longit_df=merge(longit_df,arch_M36, by="studyid", all=T)

longit_df[ longit_df == 1] <- "A"
longit_df[ longit_df == 2] <- "B"
longit_df[ longit_df == 3] <- "C"
longit_df[ longit_df == 4] <- "D"

write.table(longit_df, file="WP2_2.archetypes_longitudinal.txt",sep="\t",quote=F)

### 3 little pie's 
df_m0=data.frame(table(longit_df$arch_M0))
colnames(df_m0)=c("Archetype","Count")
df_m0$label=paste0(df_m0$Archetype," ",sapply(df_m0$Count, function(x) round(x/sum(df_m0$Count)*100)),"%")

df_m18=data.frame(table(longit_df$arch_M18))
colnames(df_m18)=c("Archetype","Count")
df_m18$label=paste0(df_m18$Archetype," ",sapply(df_m18$Count, function(x) round(x/sum(df_m18$Count)*100)),"%")

df_m36=data.frame(table(longit_df$arch_M36))
colnames(df_m36)=c("Archetype","Count")
df_m36$label=paste0(df_m36$Archetype," ",sapply(df_m36$Count, function(x) round(x/sum(df_m36$Count)*100)),"%")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols=c(gg_color_hue(4)[c(1,3,2,4)],"grey")

pdf("Fig4A.pie_charts.pdf", height=4, width=12)
par(mfrow=(c(1,3)))
pie(df_m0$Count, labels=df_m0$label, col=cols,clockwise=T, main="M0")
pie(df_m18$Count, labels=df_m18$label, col=cols,clockwise=T, main="M18")
pie(df_m36$Count, labels=df_m36$label, col=cols,clockwise=T, main="M36")
dev.off()

### what percentage stays in the same main cluster:
m0_m18=table(longit_df[,c("arch_M0","arch_M18")])/rowSums(table(longit_df[,c("arch_M0","arch_M18")]))
write.table(m0_m18, file="longitudinal.m0_m18.txt",sep="\t",quote=F)

m0_m36=table(longit_df[,c("arch_M0","arch_M36")])/rowSums(table(longit_df[,c("arch_M0","arch_M36")]))
write.table(m0_m36, file="longitudinal.m0_m36.txt",sep="\t",quote=F)

