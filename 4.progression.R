module load gcc/6.2.0
module load intel/redist/2017.2.174
module load intel/perflibs/64/2017_update2
module load R/3.4.0-ICC-MKL

R

options(width=150)
setwd("~/archetypes/")
load("WP2.2_archetypes.Rdata")

slopes=read.table("/home/Teams/teamVIP/Analysis/Roberto/20181008_WP22_HbA1c_progression_parameters_Bizzotto_FINAL/20181008_final_slopes_Bizzotto.txt",h=T, sep="\t")
names(slopes)[1]="studyid"

archetypes$studyid=rownames(archetypes)
slopes=merge(slopes,archetypes,by="studyid")

slopes_df=data.frame(archetype=colnames(archetypes)[1:4],
	p_all=sapply(1:4, function(x) summary(lm(slopes[,3+x] ~ slopes$slope))$coefficients[8]),
	beta_all=sapply(1:4, function(x) summary(lm(slopes[,3+x] ~ slopes$slope))$coefficients[2]),
	p_untreated=sapply(1:4, function(x) 
		summary(lm(slopes[which(slopes$group=="untreated"),3+x] 
			~ slopes[which(slopes$group=="untreated"),]$slope))$coefficients[8]),
	beta_untreated=sapply(1:4, function(x) 
		summary(lm(slopes[which(slopes$group=="untreated"),3+x] 
			~ slopes[which(slopes$group=="untreated"),]$slope))$coefficients[2]),
	p_changed=sapply(1:4, function(x) 
		summary(lm(slopes[which(slopes$group=="changed"),3+x] 
			~ slopes[which(slopes$group=="changed"),]$slope))$coefficients[8]),
	beta_changed=sapply(1:4, function(x) 
		summary(lm(slopes[which(slopes$group=="changed"),3+x] 
			~ slopes[which(slopes$group=="changed"),]$slope))$coefficients[2]),
	p_unchanged=sapply(1:4, function(x) 
		summary(lm(slopes[which(slopes$group=="unchanged"),3+x] 
			~ slopes[which(slopes$group=="unchanged"),]$slope))$coefficients[8]),
	beta_unchanged=sapply(1:4, function(x) 
		summary(lm(slopes[which(slopes$group=="unchanged"),3+x] 
			~ slopes[which(slopes$group=="unchanged"),]$slope))$coefficients[2]),
	p_treated=sapply(1:4, function(x) 
		summary(lm(slopes[which(slopes$group!="untreated"),3+x] 
			~ slopes[which(slopes$group!="untreated"),]$slope))$coefficients[8]),
	beta_treated=sapply(1:4, function(x) 
		summary(lm(slopes[which(slopes$group!="untreated"),3+x] 
			~ slopes[which(slopes$group!="untreated"),]$slope))$coefficients[2]))



#### differences in medication:
meds=read.table("~/pheno_tables/FINAL_221118/medication.wp2.2.21122018.txt",h=T,sep="\t")
meds=meds[,-2]

meds_by_archetype = merge(meds,archetypes, by="studyid")

### test differences between treatment groups with Wilcox test:
meds_df=data.frame(archetype=colnames(archetypes)[1:4])

med_groups=colnames(meds[c(2:4,6:7,9,12,13, 16,19,21,22,25:29)])

for (med in med_groups){
	p=sapply(colnames(archetypes)[1:4], function(a)
		wilcox.test(meds_by_archetype[,a]~meds_by_archetype[,med])$p.value)
	est=-sapply(colnames(archetypes)[1:4], function(a)
		wilcox.test(meds_by_archetype[,a]~meds_by_archetype[,med], conf.int=T)$estimate)
	meds_df=cbind(meds_df, data.frame(p=p,est=est))
	names(meds_df)[ncol(meds_df)-1]=paste0(med ,"_p")
	names(meds_df)[ncol(meds_df)]=paste0(med,"_est")
}



