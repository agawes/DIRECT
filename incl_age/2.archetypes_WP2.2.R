module load gcc/6.2.0
module load intel/redist/2017.2.174
module load intel/perflibs/64/2017_update2
module load R/3.4.0-ICC-MKL

R

library("archetypes")
library(RColorBrewer)
library(ggfortify) # for the PCA plot
library("grDevices")

setwd("~/Archetypes_Age/")

phe=read.table("WP2_2.archetype.pheno_matrix.200919.txt",sep="\t",h=T)

seed=55555

nrep=100
aa=list()
set.seed(seed)
aa<-stepArchetypes(phe, k=4, nrep=nrep, method=robustArchetypes)

### get the best model 
best=bestModel(aa)

### coefficients (soft clustering)
coef=coef(best)

### hard clusters:
### use archetype membership >0.5 as cut-off:
# clusters_0.5=sapply(1:nrow(coef), function(x) 
#	ifelse(length(which(coef[x,]>0.5))==1,which(coef[x,]>0.5),0))
# table(clusters_0.5)
# clusters_0.5
#   0   1   2   3   4 
# 268 151  41  88 178 

# clusters_0.55=sapply(1:nrow(coef), function(x) 
#	ifelse(length(which(coef[x,]>0.55))==1,which(coef[x,]>0.55),0))
#table(clusters_0.55)
# clusters_0.55
#   0   1   2   3   4 
# 367 117  34  63 145 

clusters_0.6=sapply(1:nrow(coef), function(x) 
	ifelse(length(which(coef[x,]>0.6))==1,which(coef[x,]>0.6),0))
table(clusters_0.6)
# clusters_0.6
#   0   1   2   3   4 
# 472  84  22  45 103 

### plot hard-clusters on PCA
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
gg_cols=c("grey",gg_color_hue(4))

pdf("Fig1A.archetypes_hard.PCA.pdf")
#col_0.5=sapply(clusters_0.5, function(x) gg_cols[x+1])
#autoplot(prcomp(phe), col=col_0.5, size=2) + theme_bw()
#col_0.55=sapply(clusters_0.55, function(x) gg_cols[x+1])
#autoplot(prcomp(phe), col=col_0.55, size=2) + theme_bw()
col_0.6=sapply(clusters_0.6, function(x) gg_cols[x+1])
autoplot(prcomp(phe), col=col_0.6, size=2) + theme_bw()
dev.off()

#archetypes=data.frame(studyid=rownames(phe), arch_0.5=clusters_0.5,
#	arch_0.55=clusters_0.55, arch_0.6=clusters_0.6)
archetypes=data.frame(studyid=rownames(phe),  arch_0.6=clusters_0.6)

### explore the correlations of phenotypes with archetypes:
data.frame(sort(cor(cbind(coef, phe))[1,])) ## IR  / SIRD
data.frame(sort(cor(cbind(coef, phe))[2,]))  ## BMI-driven / MOD
data.frame(sort(cor(cbind(coef, phe))[3,]))  ## BC-driven (glu.sens) /SIDD
data.frame(sort(cor(cbind(coef, phe))[4,]))	## IS-group, low insulin, PFR1, older /MARD

tmp=sapply(archetypes[,2], function(x) replace(x, x==4,"A"))
tmp=sapply(tmp, function(x) replace(x, x==3, "D"))
tmp=sapply(tmp, function(x) replace(x, x==2, "B"))
tmp=sapply(tmp, function(x) replace(x, x==1, "C"))
tmp=sapply(tmp, function(x) replace(x, x==0, "mix"))
archetypes$arch_0.6=tmp

colnames(coef)=c("C","B","D","A")
archetypes=cbind(archetypes, coef[,c(4,2,1,3)])
write.table(archetypes, "archetypes.wp2.2.txt",sep="\t",quote=F, row.names=F)


coef=data.frame(coef)
coef$studyid=rownames(phe)
coef=coef[,c(5,4,2,1,3)]
write.table(coef, "archetypes.coefficients.wp2.2.txt",sep="\t",quote=F, row.names=F)

save.image("WP2.2_archetypes.Rdata")

### PCA plot - shaded by archetype membership
### plot - soft clustering PCA - gradient of color represents contribution of archetype
coef=coef[,c("A","B","C","D")]
clusters=max.col(coef)

p1=c(rep("lightgrey",25),colorRampPalette(brewer.pal(9,"Reds"))(76))
p2=c(rep("lightgrey",25),colorRampPalette(brewer.pal(9,"Blues"))(76))
p3=c(rep("lightgrey",25),colorRampPalette(brewer.pal(9,"Greens"))(76))
p4=c(rep("lightgrey",25),colorRampPalette(brewer.pal(9,"Purples"))(76))

col_df=data.frame(a1=sapply(round(coef[,1], digits=2), function(x) p1[x*100+1] ), 
	a2=sapply(round(coef[,2], digits=2), function(x) p2[x*100+1] ),
	a3=sapply(round(coef[,3], digits=2), function(x) p3[x*100+1] ),
	a4=sapply(round(coef[,4], digits=2), function(x) p4[x*100+1] ))
col_gradient=as.character(sapply(1:length(clusters), function(x) col_df[x,clusters[x]]))

col_outline=rep("lightgrey",length(clusters))
col_outline[which(clusters_0.6==4)]="darkred"
col_outline[which(clusters_0.6==2)]="blue"
col_outline[which(clusters_0.6==1)]="darkgreen"
col_outline[which(clusters_0.6==3)]="purple"


pdf("Fig1A.archetypes.PCA.pdf")
autoplot(prcomp(phe), col=col_outline, size=3, shape=21, fill=col_gradient, stroke=3) + theme_bw() +
theme(axis.text=element_text(size=14), axis.title=element_text(size=14))
dev.off()

pdf("Fig1A.legend.pdf",height=3)
## plot archetype colour legends separately:
layout(matrix(1:4,ncol=4), width = c(1,1,1,1),height = c(1,1))

legend_a1 <- as.raster(matrix(rev(p1), ncol=1))
legend_a2 <- as.raster(matrix(rev(p2), ncol=1))
legend_a3 <- as.raster(matrix(rev(p3), ncol=1))
legend_a4 <- as.raster(matrix(rev(p4), ncol=1))

plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Archetype A')
text(x=1.5, y = seq(0,1,l=5), labels = seq(0,1,l=5), cex=1.5)
rasterImage(legend_a1, 0, 0, 1,1)

plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Archetype B')
text(x=1.5, y = seq(0,1,l=5), labels = seq(0,1,l=5),cex=1.5)
rasterImage(legend_a2, 0, 0, 1,1)

plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Archetype C')
text(x=1.5, y = seq(0,1,l=5), labels = seq(0,1,l=5), cex=1.5)
rasterImage(legend_a3, 0, 0, 1,1)

plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Archetype D')
text(x=1.5, y = seq(0,1,l=5), labels = seq(0,1,l=5), cex=1.5)
rasterImage(legend_a4, 0, 0, 1,1)
dev.off()

### small pie chart with numbers of individualss in subgroups
df=data.frame(table(archetypes$arch_0.6))
colnames(df)=c("Archetype","Count")
df$label=paste0(df$Archetype," ",sapply(df$Count, function(x) round(x/sum(df$Count)*100)),"%")



pdf("Fig1B.pie_chart.pdf", height=4, width=4)
pie(df$Count, labels=df$label, col=c("darkred", "blue", "darkgreen","purple","grey"),
	clockwise=T)

dev.off()

### compare to Ahlqvist-like clusters:

#km_res=read.table("~/LGroop_clustering/K-means.DIRECT.WP2.2.LG_clusters.txt",h=T)
km_res=read.table("/home/Teams/teamVIP/Analysis/Clustering/K-means.DIRECT.WP2.2.LG_clusters.txt",h=T)

km_res$cluster_all=factor(km_res$cluster_all, levels=c("SIDD","SIRD","MOD","MARD"))

km_res = merge(km_res,archetypes,by="studyid")

## add GADA+ individuals:
gada=read.table("/home/Teams/teamVIP/Analysis/Clustering/WP2.2_GADA_pos.txt")
archetypes$GADA=sapply(as.character(archetypes$studyid), function(x) ifelse(x %in% gada$x,1,0))

enrichment=matrix(,nrow=5,ncol=5)
for (i in 1:length(levels(km_res$cluster_all))){
  for(j in 1:length(levels(km_res$arch_0.6))){
    A=table(km_res[,c("cluster_all","arch_0.6")])[i,j]
    B=as.numeric(rowSums(table(km_res[,c("cluster_all", "arch_0.6")]))[i])
    C=nrow(km_res)
    D=as.numeric(colSums(table(km_res[,c("cluster_all","arch_0.6")]))[j])
    enrichment[i,j]=phyper(A-1,B,C-B,D,lower.tail=F)
  }
}

archetypes$GADA=factor(archetypes$GADA)
for(j in 1:length(levels(archetypes$arch_0.6))){
    A=table(archetypes[,c("GADA","arch_0.6")])[2,j]
    B=as.numeric(rowSums(table(archetypes[,c("GADA", "arch_0.6")])))[2]
    C=nrow(archetypes)
    D=as.numeric(colSums(table(archetypes[,c("GADA","arch_0.6")]))[j])
    enrichment[5,j]=phyper(A-1,B,C-B,D,lower.tail=F)
  }
colnames(enrichment)=c("A","B","C","D","mix")
rownames(enrichment)=c("SIDD","SIRD","MOD","MARD","GADA")


write.table(enrichment, file="SupTable.LG_enrichment.txt",sep="\t",quote=F)
