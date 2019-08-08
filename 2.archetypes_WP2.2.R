module load gcc/6.2.0
module load intel/redist/2017.2.174
module load intel/perflibs/64/2017_update2
module load R/3.4.0-ICC-MKL

R

library("archetypes")
library(RColorBrewer)
library(ggfortify) # for the PCA plot

setwd("~/archetypes/")

phe=read.table("../FINAL_baseline_analysis_2019/M0_FINAL.pheno_matrix.wp2_2.flt_rnt_resid.txt",sep="\t",h=T)
m=phe[,c(5:35)]

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
data.frame(sort(cor(cbind(coef, m))[1,])) ## IS-group, low insulin, PFR1?
data.frame(sort(cor(cbind(coef, m))[2,]))  ## BMI-driven
data.frame(sort(cor(cbind(coef, m))[3,]))  ## IR-driven
data.frame(sort(cor(cbind(coef, m))[4,]))	## BC-driven (glu.sens)

### plot - soft clustering PCA - gradient of color represents contribution of archetype

p1=c(rep("lightgrey",25),colorRampPalette(brewer.pal(9,"Reds"))(76))
p2=c(rep("lightgrey",25),colorRampPalette(brewer.pal(9,"Blues"))(76))
p3=c(rep("lightgrey",25),colorRampPalette(brewer.pal(9,"Greens"))(76))
p4=c(rep("lightgrey",25),colorRampPalette(brewer.pal(9,"Purples"))(76))

col_df=data.frame(a1=sapply(round(coef[,1], digits=2), function(x) p1[x*100+1] ), 
	a2=sapply(round(coef[,2], digits=2), function(x) p2[x*100+1] ),
	a3=sapply(round(coef[,3], digits=2), function(x) p3[x*100+1] ),
	a4=sapply(round(coef[,4], digits=2), function(x) p4[x*100+1] ))
col_gradient=as.character(sapply(1:length(clusters), function(x) col_df[x,clusters[x]]))

pdf("Fig1A.archetypes.PCA.pdf")
autoplot(prcomp(phe[,c(5:35)]), col=col_gradient, size=2) + theme_bw()

## plot archetype colour legends separately:
layout(matrix(1:4,ncol=4), width = c(1,1,1,1),height = c(1,1))

legend_a1 <- as.raster(matrix(rev(p1), ncol=1))
legend_a2 <- as.raster(matrix(rev(p2), ncol=1))
legend_a3 <- as.raster(matrix(rev(p3), ncol=1))
legend_a4 <- as.raster(matrix(rev(p4), ncol=1))

plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Archetype 1')
text(x=1.5, y = seq(0,1,l=5), labels = seq(0,1,l=5), cex=1.5)
rasterImage(legend_a1, 0, 0, 1,1)

plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Archetype 2')
text(x=1.5, y = seq(0,1,l=5), labels = seq(0,1,l=5),cex=1.5)
rasterImage(legend_a2, 0, 0, 1,1)

plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Archetype 3')
text(x=1.5, y = seq(0,1,l=5), labels = seq(0,1,l=5), cex=1.5)
rasterImage(legend_a3, 0, 0, 1,1)

plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Archetype 4')
text(x=1.5, y = seq(0,1,l=5), labels = seq(0,1,l=5), cex=1.5)
rasterImage(legend_a4, 0, 0, 1,1)
dev.off()


archetypes=data.frame(coef)
colnames(archetypes)=c("A1_IS","A2_BMI", "A3_IR", "A4_BC")
rownames(archetypes)=rownames(m)

save.image("WP2.2_archetypes.Rdata")
