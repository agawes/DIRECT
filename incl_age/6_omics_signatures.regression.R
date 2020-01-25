module load gcc/6.2.0
module load intel/redist/2017.2.174
module load intel/perflibs/64/2017_update2
module load R/3.4.0-ICC-MKL

R

options(width=150)
options(stringsAsFactors = FALSE)

library(ggplot2)
library(data.table)
library(annotables)
library(reshape2)
library(plyr)
library(stats)
library(GenABEL) 

options(width=150)
setwd("~/Archetypes_Age/")
out.path = "~/Archetypes_Age/omics_regression/"

archetypes=read.table("archetypes.coefficients.wp2.2.txt",h=T,sep="\t")

load("~/Archetypes_Age/omics/omics.Rdata")

	### merging omics and subgroup info
	omics <- lapply(omics, function(x) merge(archetypes,x, by.y=0, by.x="studyid"))
	lapply(omics, dim)

	## run lm for each group
	arch=colnames(archetypes)[2:5]
	for (a in arch){
		print(a)
		for (i in c(1,2,4,5)){ ## do transcriptomics & proteomics separately to annotate genes
			print(i)
			p=apply(omics[[i]][,6:ncol(omics[[i]])],2, function(x) summary(lm(omics[[i]][,a]~x))$coefficients[2,4] )
			q=p.adjust(p, method="fdr")
			beta=apply(omics[[i]][,6:ncol(omics[[i]])],2, function(x) summary(lm(omics[[i]][,a]~x))$coefficients[2,1] )
			df=data.frame(var=colnames(omics[[i]])[6:ncol(omics[[i]])], p=p, q=q, beta=beta)
			write.table(df, file=paste0(out.path,a,".",names(omics)[i],".txt"),sep="\t", quote=F)
		}
	}
	## proteomics separately
	for (a in arch){
		print(a)
		i=3
		print(i)
		p=apply(omics[[i]][,6:ncol(omics[[i]])],2, function(x) summary(lm(omics[[i]][,a]~x))$coefficients[2,4] )
		q=p.adjust(p, method="fdr")
		beta=apply(omics[[i]][,6:ncol(omics[[i]])],2, function(x) summary(lm(omics[[i]][,a]~x))$coefficients[2,1] )
		df=data.frame(var=colnames(omics[[i]])[6:ncol(omics[[i]])], p=p, q=q, beta=beta, 
			protein=sapply(colnames(omics[[i]])[6:ncol(omics[[i]])], function(x) prot.info[x,]$gene))
		write.table(df, file=paste0(out.path,a,".",names(omics)[i],".txt"),sep="\t", quote=F)
		}
	}

## transcriptomics separately
	for (a in arch){
		print(a)
		i=6
		print(i)
		p=apply(omics[[i]][,6:ncol(omics[[i]])],2, function(x) summary(lm(omics[[i]][,a]~x))$coefficients[2,4] )
		q=p.adjust(p, method="fdr")
		beta=apply(omics[[i]][,6:ncol(omics[[i]])],2, function(x) summary(lm(omics[[i]][,a]~x))$coefficients[2,1] )
		df=data.frame(var=colnames(omics[[i]])[6:ncol(omics[[i]])], p=p, q=q, beta=beta, 
			gene=sapply(colnames(omics[[i]])[6:ncol(omics[[i]])], function(x) gene.ann[x,]$gene))
		write.table(df, file=paste0(out.path,a,".",names(omics)[i],".txt"),sep="\t", quote=F)
		}
	}

	save.image(paste0(out.path,"omics.Rdata"))
