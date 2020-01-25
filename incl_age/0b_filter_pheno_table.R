module load gcc/6.2.0
module load intel/redist/2017.2.174
module load intel/perflibs/64/2017_update2
module load R/3.4.0-ICC-MKL

R

library("archetypes")
library(RColorBrewer)
library(ggfortify) # for the PCA plot

setwd("~/archetypes_cutoff/")

## read in the main filtered pheno matrix; use the same individuals
phe=read.table("../FINAL_baseline_analysis_2019/M0_FINAL.pheno_matrix.wp2_2.flt_rnt_resid.txt",sep="\t",h=T)
include_subjects=as.character(phe$studyid)
include_phenos=colnames(phe)[5:35]
## test the archetype solution, when including age as phenotype
## pheno table: /home/agatawa/pheno_tables/m0.rnt.center_sex_resid.txt

phe=read.table("~/pheno_tables/M0.rnt.resid.age_as_covariate.190919.txt",sep="\t",h=T)
rownames(phe)=phe$studyid
m=phe[include_subjects,c("age", include_phenos)]
write.table(m, file="WP2_2.archetype.pheno_matrix.200919.txt",sep="\t",quote=F)
