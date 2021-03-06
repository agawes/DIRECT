module load gcc/6.2.0
module load intel/redist/2017.2.174
module load intel/perflibs/64/2017_update2
module load R/3.4.0-ICC-MKL


R

library(ggplot2)

options(width=150)
setwd("~/Archetypes_Age/")
archetypes=read.table("archetypes.wp2.2.txt",h=T,sep="\t")
archetypes=cbind(archetypes, model.matrix(~arch_0.5+0, data=archetypes))
archetypes=cbind(archetypes, model.matrix(~arch_0.55+0, data=archetypes))
archetypes=cbind(archetypes, model.matrix(~arch_0.6+0, data=archetypes))

### read in psGRSs - generated by Anubha
psGRS=read.table("/home/Teams/teamVIP/Analysis/Share_Results/Scores.for.Agata.30thMay2019AM.txt",h=T)
names(psGRS)[1]="studyid"
psGRS=merge(archetypes, psGRS, by="studyid")

### first test overall Kruskal-Wallis differences in GRSs between groups
GRS_res=data.frame(GRS=colnames(psGRS)[24:30], 
	KW_p_0.5=sapply(24:30, function(x) kruskal.test(psGRS[,x], psGRS$arch_0.5)$p.value ),
	KW_p_0.55=sapply(24:30, function(x) kruskal.test(psGRS[,x], psGRS$arch_0.55)$p.value ),
	KW_p_0.6=sapply(24:30, function(x) kruskal.test(psGRS[,x], psGRS$arch_0.6)$p.value )
	)

### compare each subgroup to remaining individuals
grs=colnames(psGRS)[c(24:30)]

grs_wilcox = data.frame()
for (a in colnames(psGRS)[9:23]){
	for (g in grs){
		df=data.frame(archetype=a, phenotype=g, p=wilcox.test(as.numeric(psGRS[,g])~psGRS[,a])$p.value,
			est=-wilcox.test(as.numeric(psGRS[,g])~psGRS[,a], conf.int=T)$estimate)
# instead of Wilcox pseudomean estimate - use difference in means
#			est=mean(ps[groups[,gr]==1,g])-mean(groups[groups[,gr]==0,g]))
		grs_wilcox=rbind(grs_wilcox,df)
	}
}

grs_wilcox$star=cut(grs_wilcox$p, breaks=c(-Inf, 0.001, 0.01, 0.05,Inf), label=c("***", "**", "*","")) 
write.table(grs_wilcox, file="WP2.2.psGRS_AM.wilcox.txt",sep="\t",quote=F)

grs_wilcox$logP=-log10(grs_wilcox$p) * grs_wilcox$est/abs(grs_wilcox$est)

pdf("Fig3A.GRS_by_archetype.pdf",width=7,height=4)
p=ggplot(aes(x=phenotype, y=archetype, fill=logP), data=grs_wilcox)
p + geom_tile() + scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C") + 
  geom_text(aes(label=star), color="black", size=3) + scale_y_discrete(limits = rev(levels(grs_wilcox$archetype))) +
  labs(y="Archetype", x=NULL, fill="-log10P") + theme_bw()
dev.off()

## psGRS - regression

GRS_df=data.frame(GRS=grs,
	p_MARD=sapply(grs, function(x) summary(lm(psGRS[,x] ~ psGRS$MARD))$coefficients[8]),
	beta_MARD=sapply(grs, function(x) summary(lm(psGRS[,x] ~ psGRS$MARD))$coefficients[2]),
	p_MOD=sapply(grs, function(x) summary(lm(psGRS[,x] ~ psGRS$MOD))$coefficients[8]),
	beta_MOD=sapply(grs, function(x) summary(lm(psGRS[,x] ~ psGRS$MOD))$coefficients[2]),
	p_SIRD=sapply(grs, function(x) summary(lm(psGRS[,x] ~ psGRS$SIRD))$coefficients[8]),
	beta_SIRD=sapply(grs, function(x) summary(lm(psGRS[,x] ~ psGRS$SIRD))$coefficients[2]),
	p_SIDD=sapply(grs, function(x) summary(lm(psGRS[,x] ~ psGRS$SIDD))$coefficients[8]),
	beta_SIDD=sapply(grs, function(x) summary(lm(psGRS[,x] ~ psGRS$SIDD))$coefficients[2])
)
write.table(GRS_df, file="S_Table.psGRS.regression.txt",sep="\t", quote=F)

### odds-ratio style plot:
arch=c("MARD","MOD","SIRD","SIDD")

labels=apply(expand.grid(factor(arch), grs), 1, function(x) paste0(x[2],":",x[1]))
lm_coefs=as.numeric(sapply(grs, function(g) sapply(arch, function(x) lm(psGRS[,g] ~ psGRS[,x])$coef[2])))
CIlow=as.numeric(sapply(grs, function(g) sapply(arch, function(x) confint(lm(psGRS[,g] ~ psGRS[,x]))[2,1])))
CIhigh=as.numeric(sapply(grs, function(g) sapply(arch, function(x) confint(lm(psGRS[,g] ~ psGRS[,x]))[2,2])))
lm_signif=as.numeric(sapply(grs, function(g) sapply(arch, function(x) summary(lm(psGRS[,g] ~ psGRS[,x]))$coef[2,4])))
psGRS_lm_res=data.frame(label=labels, coef=lm_coefs, CIlow=CIlow, CIhigh=CIhigh, lm_signif=lm_signif)

psGRS_lm_res$col=c(rep("darkgrey",4), rep("grey",4), rep("darkgrey",4),rep("grey",4),
    	rep("darkgrey",4), rep("grey",4),rep("darkgrey",4))
psGRS_lm_res$col[which(psGRS_lm_res$lm_signif<0.05)]="darkred"

ylab=c(rev(levels(psGRS_lm_res$label))[1:4],"",rev(levels(psGRS_lm_res$label))[5:8],"",
	rev(levels(psGRS_lm_res$label))[9:12],"",rev(levels(psGRS_lm_res$label))[13:16],"",
	rev(levels(psGRS_lm_res$label))[17:20],"",rev(levels(psGRS_lm_res$label))[21:24],"",
	rev(levels(psGRS_lm_res$label))[25:28])
pdf("psGRS.regression_by_archetype.pdf", width=5)
p <- ggplot(psGRS_lm_res, aes(x = lm_coefs, y = label)) + 
    geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed") + 
    geom_errorbarh(aes(xmax = CIhigh, xmin = CIlow), size = .5, height = 0.2, color = "gray50") +
    geom_point(size = 3.5, color = psGRS_lm_res$col) +
    theme_bw()+ 
	scale_y_discrete(limits = ylab) +
    #scale_y_discrete(limits = rev(levels(psGRS_lm_res$label))) +

    theme(panel.grid.minor = element_blank()) +
    ylab("") +
    xlab("psGRS - estimate (95% CI)") 
p
dev.off()