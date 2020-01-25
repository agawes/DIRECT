module load gcc/6.2.0
module load intel/redist/2017.2.174
module load intel/perflibs/64/2017_update2
module load R/3.4.0-ICC-MKL

R

options(width=150)
setwd("~/Archetypes_Age/")
archetypes=read.table("archetypes.wp2.2.txt",h=T,sep="\t")
archetypes=cbind(archetypes, model.matrix(~arch_0.6+0, data=archetypes))

slopes=read.table("/home/Teams/teamVIP/Analysis/Roberto/20181008_WP22_HbA1c_progression_parameters_Bizzotto_FINAL/20181008_final_slopes_Bizzotto.txt",h=T, sep="\t")
names(slopes)[1]="studyid"

slopes_df=merge(slopes,archetypes,by="studyid")

subset_changed=slopes_df[slopes_df$group=="changed",]
subset_unchanged=slopes_df[slopes_df$group=="unchanged",]
subset_untreated=slopes_df[slopes_df$group=="untreated",]
subset_treated=slopes_df[slopes_df$group!="untreated",]

### create a table with differences in progression
prog_df = data.frame(group=c("all","changed","unchanged","untreated","treated"),
	KW_p_0.6=c(kruskal.test(slopes_df$slope, slopes_df$arch_0.6)$p.value, # p=0.004161871
		kruskal.test(subset_changed$slope, subset_changed$arch_0.6)$p.value,	# p=0.686204576
		kruskal.test(subset_unchanged$slope, subset_unchanged$arch_0.6)$p.value, # p=0.163958933
		kruskal.test(subset_untreated$slope, subset_untreated$arch_0.6)$p.value, # p=0.029831210
		kruskal.test(subset_treated$slope, subset_treated$arch_0.6)$p.value)) # p=0.326914443
#
#write.table(prog_df, file="Sup.Table.HbA1c_slopes.KW.txt",quote=F,sep="\t")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols=c(gg_color_hue(4)[c(1,3,2,4)],"grey")

pdf("Fig2A.WP2_2.HbA1c_slopes.boxplot.pdf",width=10)
par(mar=c(10,5,1,1),mfrow=c(1,3), bty="l")

boxplot(slopes_df$slope ~ slopes_df$arch_0.6, col=cols, cex.axis=2, lty=1, staplewex=0, boxwex=0.8, boxlwd=1, medlwd=1,
	ylab="T2D progression in all patients", cex.lab=2, bty="l", ylim=c(min(slopes$slope),max(slopes$slope)))
abline(h=0, lty=2, col="grey")
boxplot(subset_untreated$slope ~ subset_untreated$arch_0.6, col=cols, cex.axis=2, lty=1, staplewex=0, boxwex=0.8, boxlwd=1, medlwd=1,
	ylab="T2D progression in untreated patients", cex.lab=2, bty="l", ylim=c(min(slopes$slope),max(slopes$slope)))
abline(h=0, lty=2, col="grey")
boxplot(subset_treated$slope ~ subset_treated$arch_0.6, col=cols, cex.axis=2,lty=1, staplewex=0, boxwex=0.8, boxlwd=1, medlwd=1,
	ylab="T2D progression in treated patients", cex.lab=2, bty="l", ylim=c(min(slopes$slope),max(slopes$slope)) )
abline(h=0, lty=2, col="grey")
dev.off()

slopes_df$treated=rep("no",nrow(slopes_df))
slopes_df[slopes_df$group!="untreated",]$treated="yes"
slopes_df$treated=factor(slopes_df$treated)

groups=c("A","B","C","D","mix")
colnames(slopes_df)=c("studyid","group","slope","arch_0.6","coef_A","coef_B","coef_C","coef_D","A","B","C","D","mix","treated")

prog_df_2  = data.frame(group=c(paste0(groups,"_all"),paste0(groups,"_treated"),
	paste0(groups,"_untreated"),paste0(groups,"_tr_vs_untr"), 
	paste0(groups,"_tr_vs_untr+greater")),

	p=c(sapply(groups, function(x) 
	wilcox.test(slopes_df$slope~slopes_df[,x])$p.value), ## all
 		sapply(groups, function(x) 
 			wilcox.test(slopes_df[slopes_df$group!="untreated",]$slope~slopes_df[slopes_df$group!="untreated",x])$p.value),
 		sapply(groups, function(x) 
 			wilcox.test(slopes_df[slopes_df$group=="untreated",]$slope~slopes_df[slopes_df$group=="untreated",x])$p.value),
 		sapply(groups, function(x) 
 			wilcox.test(slopes_df[slopes_df$arch_0.6==x,]$slope
 				~slopes_df[slopes_df$arch_0.6==x,]$treated)$p.value),
 		sapply(groups, function(x) 
 			wilcox.test(slopes_df[slopes_df$arch_0.6==x,]$slope
 				~slopes_df[slopes_df$arch_0.6==x,]$treated, alternative="greater")$p.value)
 	),
 	est=-(c(sapply(groups, function(x) 
	wilcox.test(slopes_df$slope~slopes_df[,x], conf.int=T)$estimate), ## all
 		sapply(groups, function(x) 
 			wilcox.test(slopes_df[slopes_df$group!="untreated",]$slope~slopes_df[slopes_df$group!="untreated",x], conf.int=T)$estimate),
 		sapply(groups, function(x) 
 			wilcox.test(slopes_df[slopes_df$group=="untreated",]$slope~slopes_df[slopes_df$group=="untreated",x], conf.int=T)$estimate),
 		sapply(groups, function(x) 
 			wilcox.test(slopes_df[slopes_df$arch_0.6==x,]$slope
 				~slopes_df[slopes_df$arch_0.6==x,]$treated, conf.int=T)$estimate),
 		sapply(groups, function(x) 
 			wilcox.test(slopes_df[slopes_df$arch_0.6==x,]$slope
 				~slopes_df[slopes_df$arch_0.6==x,]$treated, alternative="greater", conf.int=T)$estimate)
 	))
 )
write.table(prog_df, file="Sup.Table2.progression_by_archetype.kruskal.txt",sep="\t",col.names=NA, quote=F)

write.table(prog_df_2, file="Sup.Table2.progression_by_archetype.wilcox.txt",sep="\t",col.names=NA, quote=F)

### make one at each cut-off:

## make heatmap bars with differences in progression, like in Figure 1
#library(ggplot2)
#prog_df$star=cut(prog_df$p, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", "")) 
#prog_df=prog_df[1:15,]
#prog_df$treatment=factor(gsub(".+_","",prog_df$group),levels=c("treated","untreated","all"))
#prog_df$subgroup=gsub("_.+","",prog_df$group)

#pdf("Fig2E.progression_heatmap.by_group.pdf",width=4,height=2)
#p <- ggplot(aes(x=subgroup, y=treatment, fill=est), data=prog_df)
#p + geom_tile() + scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C") + 
#  geom_text(aes(label=star), color="black", size=3) + 
#  labs(y=NULL, x=NULL, fill="estimate") + theme_bw() 
#dev.off()

#### differences in medication between groups
meds_by_group=read.table("~/pheno_tables/FINAL_221118/medication.wp2.2.21122018.txt",h=T,sep="\t")
meds_by_group=meds_by_group[,-2]
meds_by_group = merge(archetypes,meds_by_group, by="studyid")

### test differences between groups with Fisher's exact test:
meds_to_test=which(apply(meds_by_group[,c(12:39)],2,function(x) mean(x, na.rm=T)) != 0) + 11

fish=sapply(meds_to_test, function(x) fisher.test(table(meds_by_group[,c(2,x)]),workspace=2e8)$p.value)
names(fish)=colnames(meds_by_group[meds_to_test])

## run wilcox for each group against rest:
meds_df=data.frame(meds=names(meds_to_test), fisher=fish,
	A_p=sapply(meds_to_test, function(x) 
		fisher.test(table(meds_by_group[,c(7,x)]),workspace=2e8)$p.value),
	A_OR=sapply(meds_to_test, function(x) 
		ifelse(is.null(fisher.test(table(meds_by_group[,c(7,x)]),workspace=2e8)$estimate),
			NA,fisher.test(table(meds_by_group[,c(7,x)]),workspace=2e8)$estimate)),
	B_p=sapply(meds_to_test, function(x) 
		fisher.test(table(meds_by_group[,c(8,x)]),workspace=2e8)$p.value),
	B_OR=sapply(meds_to_test, function(x) 
		ifelse(is.null(fisher.test(table(meds_by_group[,c(8,x)]),workspace=2e8)$estimate),
			NA,fisher.test(table(meds_by_group[,c(8,x)]),workspace=2e8)$estimate)),
	C_p=sapply(meds_to_test, function(x) 
		fisher.test(table(meds_by_group[,c(9,x)]),workspace=2e8)$p.value),
	C_OR=sapply(meds_to_test, function(x) 
		ifelse(is.null(fisher.test(table(meds_by_group[,c(9,x)]),workspace=2e8)$estimate),
			NA,fisher.test(table(meds_by_group[,c(9,x)]),workspace=2e8)$estimate)),
	D_p=sapply(meds_to_test, function(x) 
		fisher.test(table(meds_by_group[,c(10,x)]),workspace=2e8)$p.value),
	D_OR=sapply(meds_to_test, function(x) 
		ifelse(is.null(fisher.test(table(meds_by_group[,c(10,x)]),workspace=2e8)$estimate),
			NA,fisher.test(table(meds_by_group[,c(10,x)]),workspace=2e8)$estimate)),
	mix_p=sapply(meds_to_test, function(x) 
		fisher.test(table(meds_by_group[,c(11,x)]),workspace=2e8)$p.value),
	mix_OR=sapply(meds_to_test, function(x) 
		ifelse(is.null(fisher.test(table(meds_by_group[,c(11,x)]),workspace=2e8)$estimate),
			NA,fisher.test(table(meds_by_group[,c(11,x)]),workspace=2e8)$estimate))
)
write.table(meds_df, file="Sup.Table.meds_by_group.fisher.txt",sep="\t",quote=F)

# plot
pdf("Fig2B.meds_by_group.incl_gray.pdf",width=11)
par(mar=c(5,5,1,1),mfrow=c(1,2), oma=c(0,0,0,10), xpd=NA)
## arch_0.6
plot(1,type="n",xlab="",ylab="% patients on T2D medication",xlim=c(0.5,3.5),ylim=c(0,90),xaxt="n",cex.lab=2,
	cex.axis=2, cex=2,bty="l")
axis(1,at=c(1,2,3),labels=c("0m","18m","36m"), cex.axis=2)
for (i in 35:37){
  points(rep(i-34,5),100*as.numeric(table(meds_by_group[,c(2,i)])[,2]/rowSums(table(meds_by_group[,c(2,i)]))),
  	col="black",bg=cols, pch=21, cex=2)
}

plot(1,type="n",xlab="",ylab="% patients with T2D medication change",xlim=c(1.5,3.5),ylim=c(0,90),xaxt="n",cex.lab=2,
	cex.axis=2, cex=2, bty="l")
axis(1,at=c(2,3),labels=c("18m","36m"), cex.axis=2)
for (i in 38:39){
  points(rep(i-36,5),100*as.numeric(table(meds_by_group[,c(2,i)])[,2]/rowSums(table(meds_by_group[,c(2,i)]))),col="black",
  	bg=cols, pch=21, cex=2)
}

legend(3.5,85, c("A","B","C","D","mix"), title="Archetype",col="black",pt.bg=cols,pch=21,
 pt.cex=2, cex=2, bty="n")

dev.off()

##### statins
statins=read.table("~/pheno_tables/WP2.2_statins-2018-05-18.txt",h=T,sep="\t")
names(statins)[1]="studyid"

archetypes$statins=rep(0,nrow(archetypes))
for (i in 1:nrow(statins)){
  if ((statins[i,]$crf == "scr" | statins[i,]$crf == "0m") & as.character(statins[i,]$studyid) %in% archetypes$studyid){
  	archetypes[which(archetypes$studyid==as.character(statins[i,]$studyid)),]$statins=1 }
}
fisher.test(table(archetypes$statins, archetypes$arch_0.5), workspace=2e8)  
table(archetypes$statins, archetypes$arch_0.5)
## at cut-0.6 not significant

#### differences in diabetes complications
comp=read.table("~/pheno_tables/diabetes_complications.CRFs.140518.txt",h=T,sep="\t")
comp_m=merge(comp, archetypes, by="studyid")
comp_m$arch_0.5=factor(comp_m$arch_0.5)

### test differences between groups with Fisher's exact test:
fish_comp=sapply(c(3:6), function(x) fisher.test(table(comp_m[,c(10,x)]))$p.value)
names(fish_comp)=colnames(comp_m)[3:6]
table(comp_m[,c(10,3)])

#         retinopathy_0m
# arch_0.5   0   1
#       0 265   3
#       1 143   8
#       2  41   0
#       3  80   8
#       4 176   2

