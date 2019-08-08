module load gcc/6.2.0
module load intel/redist/2017.2.174
module load intel/perflibs/64/2017_update2
module load R/3.4.0-ICC-MKL

R

library(archetypes)
library(RColorBrewer)
library(ggfortify) # for the PCA plot
library(fossil)

setwd("~/archetypes/longitudinal/")
phe=read.table("~/pheno_tables/FINAL_221118/M0_18_36.rnt.resid.221118.txt",h=T,sep="\t")
phe$m=gsub(".+_","",rownames(phe))
phe$m=factor(phe$m)

### remove any subjects with NA values:
na=numeric()
for (i in 1:ncol(phe)){
	na = c(na, which(is.na(phe[,i])))
}
na=unique(na)
phe=phe[-na,]

## exclude 28-2-394_M18 - inconsistent BMI information at 18M:
phe=phe[which(rownames(phe) != "28-2-394_M18"),]
table(phe$m)
# M0 M18 M36 
# 726 591 423 
table(table(phe$studyid))

write.table(phe, file="FINAL.pheno_matrix.wp2_2.m0_m18_m36.txt",sep="\t",quote=F)

#### run the archetypes analysis on the full longitudinal dataset

#### what's the optimal archetype number:
m=phe[,c(5:35)]

### let's run this for 10 different seeds, and see if we always get similar answers
seeds=sample(1:1000000, 10)

nrep=10
max_k=10
aa=list()
for (s in seeds){
  print(paste0("Seed: ",s))
  set.seed(s)
  ### which method to use: archetypes, wighted Archetypes, or robustArchetypes
  aa[[s]]<-stepArchetypes(m, k=1:max_k, nrep=nrep, method=robustArchetypes)
}

pdf("pick_archetypes_number.screeplot.1.pdf")
sapply(seeds, function(s) print(screeplot(aa[[s]], main=s)))
dev.off()

  ### evaluate stability of the clusters:
 ind=list()
 adj.ind=list()

for (s in seeds){
	ind[[s]]=list()
	adj.ind[[s]]=list()
  for (i in 1:max_k){ ## for each k
    ind[[s]][[i]]=rep(NA,nrep)
    adj.ind[[s]][[i]]=rep(NA,nrep)
    
    best=as.numeric(which(rss(aa[[s]])[i,]==min(rss(aa[[s]])[i,]))[1]) ## find best solution for each seed, and each k to use as reference
    if (!is.na(best)){
    m_best=max.col(coef(aa[[s]][[i]][[best]]))
    for (j in 1:nrep){ ## for all other models
      print(paste(s,i,j))
      if (best != j){
        m2=max.col(coef(aa[[s]][[i]][[j]]))
        adj.ind[[s]][[i]][j]=adj.rand.index(m_best,m2)   ### find better way of evaluating stability here!
        ind[[s]][[i]][j]=rand.index(m_best,m2)   
        
      }
    }
    }
  }
}
 
pdf("rand_index.pdf")
for (s in seeds){
  plot(1:max_k, sapply(ind[[s]],function(x) mean(x, na.rm=T)), type="l", ylim=c(0,1), xlab="Archetypes",ylab="Rand index", main=s)
  lines(1:max_k, sapply(adj.ind[[s]],function(x) mean(x, na.rm=T)), col="red")
  legend("topright", c("Standard","Adjusted"), col=c("black","red"),lty=1, bty="n")
}
 dev.off() 


 ### it's not as convincing as with the intitial dataset - different seeds point to 3-5 groups
 	### might be better to rerun this with different archetype membership thresholds
### but let's continue with 4 to be consistent

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
data.frame(sort(cor(cbind(coef, m))[1,])) ## BC-driven (glu.sens)
data.frame(sort(cor(cbind(coef, m))[2,]))  ## IR-driven
data.frame(sort(cor(cbind(coef, m))[3,]))  ## IS-group, low insulin, PFR1?
data.frame(sort(cor(cbind(coef, m))[4,]))	## BMI-driven

## to compare with baseline clustering:
# A1 ## IS-group, low insulin, PFR1? -> red
# A2 ## BMI-driven -> blue
# A3 ## IR-driven -> green
# A4 ## BC-driven (glu.sens) --> purple

### reorder the coef matrix:
coef=coef[,c(3,4,2,1)]
clusters_copy=rep(0,length(clusters))
clusters_copy[which(clusters==1)]<-4
clusters_copy[which(clusters==2)]<-3
clusters_copy[which(clusters==3)]<-1
clusters_copy[which(clusters==4)]<-2
clusters=clusters_copy

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

pdf("Fig4A.longitudinal_archetypes.PCA.pdf")
autoplot(prcomp(m), col=col_gradient, size=2) + theme_bw()

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

##### test how many subjects are still defined by the same major archetype
longit_df=data.frame(studyid=as.character(unique(phe$studyid)) )

arch_M0=data.frame(coef[which(phe$m=="M0"),])
arch_M0$studyid=phe[which(phe$m=="M0"),]$studyid
arch_M0$dom_arch_M0=clusters[which(phe$m=="M0")]
colnames(arch_M0)[1:4]=c("A1_M0","A2_M0","A3_M0","A4_M0")
longit_df=merge(longit_df,arch_M0, by="studyid", all=T)

arch_M18=data.frame(coef[which(phe$m=="M18"),])
arch_M18$studyid=phe[which(phe$m=="M18"),]$studyid
arch_M18$dom_arch_M18=clusters[which(phe$m=="M18")]
colnames(arch_M18)[1:4]=c("A1_M18","A2_M18","A3_M18","A4_M18")
longit_df=merge(longit_df,arch_M18, by="studyid", all=T)

arch_M36=data.frame(coef[which(phe$m=="M36"),])
arch_M36$studyid=phe[which(phe$m=="M36"),]$studyid
arch_M36$dom_arch_M36=clusters[which(phe$m=="M36")]
colnames(arch_M36)[1:4]=c("A1_M36","A2_M36","A3_M36","A4_M36")
longit_df=merge(longit_df,arch_M36, by="studyid", all=T)

table(longit_df$dom_arch_M0)/length(which(!(is.na(longit_df$dom_arch_M0))))
#         1         2         3         4 
# 0.3071625 0.2079890 0.3677686 0.1170799 
table(longit_df$dom_arch_M18)/length(which(!(is.na(longit_df$dom_arch_M18))))
#         1         2         3         4 
# 0.2707276 0.1116751 0.4145516 0.2030457 
table(longit_df$dom_arch_M36)/length(which(!(is.na(longit_df$dom_arch_M36))))
#         1         2         3         4 
# 0.2434988 0.1323877 0.3900709 0.2340426 

save.image("longitudinal.Rdata")

### what percentage stays in the same main cluster:
table(longit_df[,c("dom_arch_M0","dom_arch_M18")])/rowSums(table(longit_df[,c("dom_arch_M0","dom_arch_M18")]))
#            dom_arch_M18
# dom_arch_M0          1          2          3          4
#           1 0.66091954 0.08045977 0.14367816 0.11494253
#           2 0.12264151 0.26415094 0.40566038 0.20754717
#           3 0.04950495 0.09405941 0.71782178 0.13861386
#           4 0.15942029 0.05797101 0.11594203 0.66666667

table(longit_df[,c("dom_arch_M0","dom_arch_M36")])/rowSums(table(longit_df[,c("dom_arch_M0","dom_arch_M36")]))
#           dom_arch_M36
# dom_arch_M0          1          2          3          4
#           1 0.59375000 0.10937500 0.12500000 0.17187500
#           2 0.16216216 0.22972973 0.39189189 0.21621622
#           3 0.06040268 0.12751678 0.65100671 0.16107383
#           4 0.06976744 0.04651163 0.13953488 0.74418605

library(GGally)
pdf("M0_M18.archetypes_cor.pdf")
ggpairs(longit_df[,c(2:5,7:10)])
dev.off()

pdf("M0_M36.archetypes_cor.pdf")
ggpairs(longit_df[,c(2:5,12:15)])
dev.off()


### calculate Euclidean distance between the same individual in 4-dim space:
rownames(coef)=rownames(m)
d=as.matrix(dist(coef))


# pdf("delta_archetypes_over_time.hist.pdf")
# par(mfrow=c(2,2))
# plot(density(longit_df$A1_M18 - longit_df$A1_M0, na.rm=T), col="red")
# lines(density(longit_df$A1_M36 - longit_df$A1_M0, na.rm=T), col="darkred")

# plot(density(longit_df$A2_M18 - longit_df$A2_M0, na.rm=T), col="turquoise")
# lines(density(longit_df$A2_M36 - longit_df$A2_M0, na.rm=T), col="blue")
# plot(density(longit_df$A3_M18 - longit_df$A3_M0, na.rm=T), col="green")
# lines(density(longit_df$A3_M36 - longit_df$A3_M0, na.rm=T), col="darkgreen")
# plot(density(longit_df$A4_M18 - longit_df$A4_M0, na.rm=T), col="purple")
# lines(density(longit_df$A4_M36 - longit_df$A4_M0, na.rm=T), col="darkorchid3")
# par(mfrow=c(1,2))

# plot(density(longit_df$A1_M18 - longit_df$A1_M0, na.rm=T), col="red", ylim=c(0,6), xlim=c(-1,1))
# lines(density(longit_df$A2_M18 - longit_df$A2_M0, na.rm=T), col="turquoise")
# lines(density(longit_df$A3_M18 - longit_df$A3_M0, na.rm=T), col="green")
# lines(density(longit_df$A4_M18 - longit_df$A4_M0, na.rm=T), col="purple")

# plot(density(longit_df$A1_M36 - longit_df$A1_M0, na.rm=T), col="red", ylim=c(0,6), xlim=c(-1,1))
# lines(density(longit_df$A2_M36 - longit_df$A2_M0, na.rm=T), col="turquoise")
# lines(density(longit_df$A3_M36 - longit_df$A3_M0, na.rm=T), col="green")
# lines(density(longit_df$A4_M36 - longit_df$A4_M0, na.rm=T), col="purple")

## boxplot of the deltas - does not look very convincing despite Kruskal p-val
pdf("delta_archetypes_over_time.hist.pdf")
delta_archetype=data.frame(group=c(rep("A1_M0-M18",nrow(longit_df)),rep("A2_M0-M18",nrow(longit_df)),
	rep("A3_M0-M18",nrow(longit_df)),rep("A4_M0-M18",nrow(longit_df)),rep("A1_M0-M36",nrow(longit_df)),
	rep("A2_M0-M36",nrow(longit_df)),rep("A3_M0-M36",nrow(longit_df)),rep("A4_M0-M36",nrow(longit_df))), 
	delta=c(longit_df$A1_M18 - longit_df$A1_M0, longit_df$A2_M18 - longit_df$A2_M0,
		longit_df$A3_M18 - longit_df$A3_M0, longit_df$A4_M18 - longit_df$A4_M0,
		longit_df$A1_M36 - longit_df$A1_M0, longit_df$A2_M36 - longit_df$A2_M0,
		longit_df$A3_M36 - longit_df$A3_M0, longit_df$A4_M36 - longit_df$A4_M0))
boxplot(delta_archetype$delta~delta_archetype$group)
ggplot(delta_archetype, aes(x=group, y=delta))+geom_boxplot() +theme_bw()
dev.off()

kruskal.test(delta_archetype$delta~delta_archetype$group)

# 	Kruskal-Wallis rank sum test
# data:  delta_archetype$delta by delta_archetype$group
# Kruskal-Wallis chi-squared = 235.17, df = 7, p-value < 2.2e-16
library(dunn.test)
dunn.test(delta_archetype$delta,delta_archetype$group)

