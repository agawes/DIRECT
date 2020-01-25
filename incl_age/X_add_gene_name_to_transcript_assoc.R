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
setwd("~/Archetypes_Age/omics")

## add gene annotation to transcript asscoaition results
files <- list.files(path=".", pattern = "*.transcript.txt")

for (f in files){
	print(f)
	df=read.table(f,h=T,sep="\t",row.names=1)
	df$gene=sapply(strsplit(df$var, split="\\."), function(x) 
		grch37[which(grch37$ensgene == x[1]),]$symbol[1])
	print(head(df))
	write.table(df, file=f,sep="\t", quote=F)
}


	
