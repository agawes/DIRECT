module load gcc/6.2.0
module load intel/redist/2017.2.174
module load intel/perflibs/64/2017_update2
module load R/3.4.0-ICC-MKL

R

library(ggplot2)

options(width=150)
setwd("~/Archetypes_Age/")
archetypes=read.table("archetypes.coefficients.wp2.2.txt",h=T,sep="\t")

slopes=read.table("/home/Teams/teamVIP/Analysis/Roberto/20181008_WP22_HbA1c_progression_parameters_Bizzotto_FINAL/20181008_final_slopes_Bizzotto.txt",h=T, sep="\t")
names(slopes)[1]="studyid"

slopes=merge(slopes,archetypes,by="studyid")

phe=read.table("WP2_2.archetype.pheno_matrix.200919.txt",sep="\t",h=T)
phe$studyid=rownames(phe)
merged=merge(slopes, phe, by="studyid",all=T)

test=colnames(merged)[4:39]
## rescale all variables between -1 and 1:
library("scales")
merged[,test]=apply(merged[,test], 2, function(x) rescale(x, to=c(-1,1)))
slopes_df=data.frame(test=test,
	p=sapply(test, function(x) summary(lm(merged[,x] ~ merged$slope))$coefficients[8]),
	beta=sapply(test, function(x) summary(lm(merged[,x] ~ merged$slope))$coefficients[2]),
	CIlow=sapply(test, function(x) confint(lm(merged[,x] ~ merged$slope))[2,1]),
	CIhigh=sapply(test, function(x) confint(lm(merged[,x] ~ merged$slope))[2,2])
)
slopes_df$col=rep("grey",nrow(slopes_df))
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
slopes_df$col[1:4]=gg_color_hue(4)[c(1,3,2,4)]
slopes_df$test=factor(slopes_df$test, 
	levels=rev(as.character(slopes_df[order(slopes_df$beta),]$test)))

pdf("progression_estimates_by_archetype_and_phenotype.pdf",width=4)
p <- ggplot(slopes_df, aes(x = beta, y = test)) + 
    geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed") + 
    geom_errorbarh(aes(xmax = CIhigh, xmin = CIlow), size = .5, height = 0.2, color = "gray50") +
    geom_point(size = 3.5, color = slopes_df$col, 
    	alpha=ifelse(slopes_df$p<0.05,1,0.5)) +
    theme_bw()+
    theme(panel.grid.minor = element_blank()) +
    theme(axis.text.y = element_text(hjust=0)) +
    ylab("") +
    xlab("Slope - estimate (95% CI)") #+
p
dev.off()
