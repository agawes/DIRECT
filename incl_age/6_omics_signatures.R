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
out.path = "~/Archetypes_Age/omics/"

archetypes=read.table("archetypes.wp2.2.txt",h=T,sep="\t")
archetypes=cbind(archetypes, model.matrix(~arch_0.5+0, data=archetypes))
archetypes=cbind(archetypes, model.matrix(~arch_0.55+0, data=archetypes))
archetypes=cbind(archetypes, model.matrix(~arch_0.6+0, data=archetypes))

### path to omics files: using residuals adjusted for technical + biological covariates and rank normalized
filepath.targeted.metab = "/home/Teams/teamVIP/Analysis/MultiOmics/Data/Targeted_Metabolites/Residuals_Technical_AgeSex_DIRECT_WP2.2_TargetedMetabolites_12032019.rank.tab"
filepath.targeted.metab.old = "/home/Teams/teamVIP/Analysis/MultiOmics/Data/Targeted_Metabolites/Residuals_Technical_AgeSex_DIRECT_WP2.2_TargetedMetabolites_09022018.rank.tab"

filepath.untargeted.metab = "/home/Teams/teamVIP/Analysis/MultiOmics/Data/Untargeted_Metabolites/Metabolomic_WP2.2_DIRECT_BiolTechCenter_Release2_3Aug17_11122017_inverseNormalise.txt"
filepath.proteomics = "/home/Teams/teamVIP/Analysis/MultiOmics/Data/Proteomics/Proteomic_WP22_DIRECT_BiolTechCenter_Release2_3Aug17_02032018_ResidualsInverseNormalised.txt"
filepath.myriad = "~/monocle_0m_18m/omics_local/WP2.2-Myriadx_hsCRP.log.lm_age_sex_center.resid.tsv"
filepath.transcriptomics = "/home/Teams/teamVIP/Analysis/MultiOmics/Data/Transcriptomics/RPKM_WP2.2_DIRECT_TechAgeGenderResiduals_24082017.rank.bedv6.gz"
filepath.transcriptomics.key = "/home/Teams/teamVIP/Analysis/MultiOmics/Metadata/DIRECT_Releasev2_3August2017.tab"
filepath.olink = "/home/Data/Repository/Proteomics/WP2/Olink/"

## metadata:
metadata=read.table("~/pheno_tables/FINAL_221118/M0.raw.221118.txt",h=T,sep="\t")

##### Caroline's code - processing Olink:
### OLINK Gene identifiers
# make list of all OLINK metadata files and read data into a list of dataframes using plyr package
my_files = list.files(path=filepath.olink, pattern = "*_binder.csv", include.dirs = TRUE, full.names=TRUE)
my_files=my_files[-1]
metadata.olink = llply(my_files, function(x) read.delim(x, header=TRUE, check.names=FALSE, sep=","))

# merge all 
metadata.olink <- do.call("rbind", metadata.olink)
# make new identifier
metadata.olink$assayID <- paste(metadata.olink$id, metadata.olink$Assay, sep="_")
# rename columns with empty spaces 
names(metadata.olink) <- c("id","Panel","Version","Assay","Uniprot_ID","LOD","MissingFreq","assayID")
# collect a list of binders with too many missing data
remove.olink <- subset(metadata.olink, MissingFreq > 0.99)$id	# 16binders

### create separate list for Olink data
my_files <- list.files(path=filepath.olink, pattern = "*_NPX.csv", include.dirs = TRUE, full.names=TRUE)
my_files=my_files[-1]
olink <- llply(my_files, function(x) read.delim(x, header=TRUE, check.names=FALSE, sep=","))
names(olink) <- c("olink.cam","olink.cvd2","olink.cvd3","olink.dev","olink.met")	# 3100, 3098, 3100, 3100 and 96 vars

# remove samples which did not pass QC
olink <- lapply(olink, function(x) subset(x, `QC warning` == "Pass"))	
lapply(olink, function(x) dim(x))	# 3063,3032,3073,3076,3052 and 96 vars
# change StudyID to studyid and set as rownames
for (i in 1:5){
names(olink[[i]])[names(olink[[i]]) == 'StudyID'] <- 'studyid'
}
olink <- lapply(olink, function(x){ row.names(x) <- x$studyid; x})

olink_panel=data.frame(ab=character(), panel=character())
for (i in 1:5){
	df=data.frame(ab=names(olink[[i]])[-c(1:4)],panel=rep(names(olink)[i],ncol(olink[[i]])-4))
	olink_panel=rbind(olink_panel, df)
}
olink_panel$panel=gsub("olink.","",olink_panel$panel)
olink_panel$Assay=as.character(sapply(olink_panel$ab, function(x) subset(metadata.olink, id==x, Assay)))
#write.table(olink_panel,file="olink.panel_info.txt",sep="\t",quote=F,row.names=F)

# collect all panels in one df (remove technical columns)
olink.all <- Reduce(function(x,y) merge(x,y, by = "studyid", all = TRUE), 
list(olink[[1]][-c(2:4)],olink[[2]][-c(2:4)],olink[[3]][-c(2:4)],olink[[4]][-c(2:4)]))	# 3100  369	

# remove binders with too many missing data

olink.all =olink.all[,which(!(names(olink.all) %in% remove.olink))]

	# rename assay id with gene names in metadata.olink
	names(olink.all) <- metadata.olink$Assay[match(names(olink.all), metadata.olink$id)]
	names(olink.all) <- c('studyid', names(olink.all[-1]))	# put back studyid name
	row.names(olink.all) <- olink.all$studyid
	olink.all=olink.all[,-1]

### rntransform & residualize
	olink.all=olink.all[which(rownames(olink.all) %in% metadata$studyid),]
	olink.rnt=apply(data.frame(olink.all),2,rntransform)
	rownames(olink.rnt)=rownames(olink.all)
	covs=data.frame(studyid=rownames(olink.rnt), center=metadata$center[match(rownames(olink.rnt), metadata$studyid)], 
		age=metadata$age[match(rownames(olink.rnt), metadata$studyid)],
		sex=metadata$CRF_2_100_eligibility.sex[match(rownames(olink.rnt), metadata$studyid)] )
	olink.res=apply(olink.rnt,2, function(x) resid(lm(x~covs$center+covs$age+covs$sex, na.action=na.exclude)))

omics <- list()

	## metabolites. 
	omics[[1]] <- read.table(filepath.targeted.metab, header=TRUE, check.names=FALSE)
	# set metabolite names as row.names
	row.names(omics[[1]]) <- omics[[1]]$Metabolite
	# save LOD % for each metab for later
	lods <- omics[[1]][,c("Metabolite","PercentageLOD")]		
	# transpose metabolite object to get individuals as rows and metabs as columns
	omics[[1]] <- t(omics[[1]][-c(1:2)])	# 797 ind, 119 metabs
	
	## untargeted metabolites
	omics[[2]] <- read.table(filepath.untargeted.metab, header=TRUE, check.names=FALSE)
	# set BIOCHEMICAL as row.names (biochemical is long!)
	row.names(omics[[2]]) <- omics[[2]]$BIOCHEMICAL
	# save IDs and missing info for later
	metab.un.info <- omics[[2]][,c("COMP_ID","BIOCHEMICAL","N_samples_noNAs")]
	# transpose 
	omics[[2]] = t(omics[[2]][-c(1:3)])	# 772 238
	
	## proteomics
	# antibody panel
	omics[[3]] <- read.delim(filepath.proteomics, header=TRUE, check.names=FALSE)
	# set AB as row.names
	row.names(omics[[3]]) <- omics[[3]]$AB
	# save  gene IDs and antibody info for later
	prot.info <- omics[[3]][,c("ENSEMBL_ID","AB")]
	# transpose 
	omics[[3]] <- t(omics[[3]][-c(1:2)])	# 789 377
	# add HGNC name:
	prot.info$gene=sapply(strsplit(prot.info$ENSEMBL_ID, split=";"), 
		function(x) paste(sapply(x, function(y) grch37[which(grch37$ensgene == y),]$symbol[1]),collapse=","))
	# change name to ensembl id
	## these are not unique!
	#colnames(omics[[3]]) <- prot.info$gene[match(colnames(omics[[3]]), prot.info$AB)]

	## myriad data
	omics[[4]] <- read.delim(filepath.myriad, header=TRUE, check.names=FALSE, sep="\t") # 3101, 17 (studyid in V2)
	row.names(omics[[4]]) <- omics[[4]]$studyid
	# remove sampleid variable
	omics[[4]] <- omics[[4]][,-c(1:2)]

	## OLINK
	omics[[5]] <- olink.res 	# remove studyid

	## transcriptomics
	omics[[6]] = read.table(gzfile(filepath.transcriptomics), header=FALSE, check.names=FALSE)	# 16208   801
	id.key = read.table(filepath.transcriptomics.key, header=TRUE, check.names=FALSE)
	# keep only WP2.2 individuals that pass inclusion/exclusion criteria
	id.key = subset(id.key, WP == "WP2_2")	# 795  24 (fits with 795 + 6 info variables = 801 in transcriptomics data
	# no studyid included in file, so need to retrieve from id.key	
	names(omics[[6]]) = c("CHR","START","STOP","ENG1","ENG2","STRAND",id.key$NewGenotypeID_Match2)
	# set ENSgeneID as row.names	
	row.names(omics[[6]]) = omics[[6]]$ENG1	
	# transpose
	omics[[6]] <- t(omics[[6]][-c(1:6)])	# 795 16209
	# keep individuals that pass inclusion/exclusion criteria
	id.ok = subset(id.key, Pass_InclusionExclusion == "Yes")	# 789  24
	omics[[6]] = subset(omics[[6]], rownames(omics[[6]]) %in% id.ok$NewGenotypeID_Match2) # 789 16209
	gene.ann=data.frame(ens=colnames(omics[[6]]), gene=sapply(strsplit(colnames(omics[[6]]), split="\\."),
	 function(x) grch37[which(grch37$ensgene == x[1]),]$symbol[1]))
	rownames(gene.ann)=gene.ann$ens

	names(omics) <- c("metab.targ","metab.untarg","prot","myriad", "OLINK","transcript")	
	lapply(omics, function(x) dim(x))

save(omics, file=paste0(out.path,"omics.Rdata"))

	### merging omics and subgroup info
	omics <- lapply(omics, function(x) merge(archetypes,x, by.y=0, by.x="studyid"))
	lapply(omics, dim)

	## run wilcox for each group - do proteomics and transcriptomics separately; need to add ID's
	arch=colnames(archetypes)[9:23]
	for (a in arch){
		print(a)
		for (i in c(1,2,4,5)){
			print(i)
			p=sapply(24:ncol(omics[[i]]), function(x) 
				wilcox.test(as.numeric(omics[[i]][,x])~omics[[i]][,a])$p.value)
			q=p.adjust(p, method="fdr")
			est=sapply(24:ncol(omics[[i]]), function(x) 
				-wilcox.test(as.numeric(omics[[i]][,x])~omics[[i]][,a], conf.int=T)$estimate)
			df=data.frame(var=colnames(omics[[i]])[24:ncol(omics[[i]])], p=p, q=q, est=est)
			write.table(df, file=paste0(out.path,a,".",names(omics)[i],".txt"),sep="\t", quote=F)
		}
	}

	## proteomics separately
	for (a in arch){
		print(a)
		i=3
		print(i)
		p=sapply(24:ncol(omics[[i]]), function(x) 
				wilcox.test(as.numeric(omics[[i]][,x])~omics[[i]][,a])$p.value)
		q=p.adjust(p, method="fdr")
		est=sapply(24:ncol(omics[[i]]), function(x) 
				-wilcox.test(as.numeric(omics[[i]][,x])~omics[[i]][,a], conf.int=T)$estimate)
		df=data.frame(var=colnames(omics[[i]])[24:ncol(omics[[i]])], p=p, q=q, est=est,
			protein=sapply(colnames(omics[[i]])[24:ncol(omics[[i]])], function(x) prot.info[x,]$gene))

		write.table(df, file=paste0(out.path,a,".",names(omics)[i],".txt"),sep="\t", quote=F)
	}

## transcriptomics separately
	for (a in arch){
		print(a)
		i=6
		print(i)
		p=sapply(14:ncol(omics[[i]]), function(x) 
				wilcox.test(as.numeric(omics[[i]][,x])~omics[[i]][,a])$p.value)
		q=p.adjust(p, method="fdr")
		est=sapply(14:ncol(omics[[i]]), function(x) 
				-wilcox.test(as.numeric(omics[[i]][,x])~omics[[i]][,a], conf.int=T)$estimate)
		df=data.frame(var=colnames(omics[[i]])[14:ncol(omics[[i]])], p=p, q=q, est=est,
			gene=sapply(colnames(omics[[i]])[14:ncol(omics[[i]])], function(x) gene.ann[x,]$gene))
		write.table(df, file=paste0(out.path,a,".",names(omics)[i],".txt"),sep="\t", quote=F)
		}
	}

	save.image(paste0(out.path,"omics.Rdata"))
