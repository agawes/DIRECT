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
library(tidyverse)

options(width=150)
setwd("~/Archetypes_Age/")
out.path = "~/Archetypes_Age/omics_vs_phenotypes/"
load("~/Archetypes_Age/omics/omics.Rdata")
names(omics)=c("Biocrates","Metabolon","AB_array","Myriad","Olink","Transcript")
phe=read.table("~/Archetypes_Age/WP2_2.archetype.pheno_matrix.200919.txt",sep="\t",h=T)

### read in the top omics
setwd("~/Archetypes_Age/omics/arch_0.6/")

top_files=list.files(".",pattern="^top_*")
top_list=list()
for (i in top_files){
  top_list[[i]]=scan(i,what="character")
}

## create a top omics matrix:
names(top_list)=c("Transcript","Myriad","Olink","AB_array","Biocrates","Metabolon")
subjects=rownames(phe)
complete_subjects=intersect(subjects[which(subjects %in% rownames(omics[["AB_array"]]))], 
	subjects[which(subjects %in% rownames(omics[["Metabolon"]]))])
m=list()
for (om in names(top_list)){
	m[[om]]=omics[[om]][complete_subjects,top_list[[om]]]
}

## rename some omics
## transcriptomics - ENSG to gene symbol
colnames(m[["Transcript"]]) = sapply(gsub("\\..+","",colnames(m[["Transcript"]])), function(x)
	grch37[which(grch37$ensgene==x),]$symbol)
## prot - AB to ENSG to gene symbol
filepath.proteomics = "/home/Teams/teamVIP/Analysis/MultiOmics/Data/Proteomics/Proteomic_WP22_DIRECT_BiolTechCenter_Release2_3Aug17_02032018_ResidualsInverseNormalised.txt"
prot <- read.delim(filepath.proteomics, header=TRUE, check.names=FALSE)
prot.info <- prot[,c("ENSEMBL_ID","AB")]
prot.info$gene=sapply(strsplit(prot.info$ENSEMBL_ID, split=";"), 
		function(x) paste(sapply(x, function(y) grch37[which(grch37$ensgene == y),]$symbol[1]),collapse=","))

colnames(m[["AB_array"]]) = sapply(colnames(m[["AB_array"]]), function(x)
	prot.info[which(prot.info$AB==x),]$gene)
colnames(m[["AB_array"]])[grepl("CXCL1,CXCL3,CXCL2",colnames(m[["AB_array"]]))] = "CXCL1_2_3"

colnames(m[["Myriad"]])[grepl("Ferritin..FRTN.",colnames(m[["Myriad"]]))] = "Ferritin"
colnames(m[["Myriad"]])[grepl("Interleukin.18..IL.18.",colnames(m[["Myriad"]]))] = "IL-18"
colnames(m[["Myriad"]])[grepl("Interleukin.8..IL.8.",colnames(m[["Myriad"]]))] = "IL-8"
colnames(m[["Myriad"]])[grepl("Macrophage.Inflammatory.Protein.1.beta..MIP.1.beta.",colnames(m[["Myriad"]]))] = "MIP-1-beta"
colnames(m[["Myriad"]])[grepl("Monocyte.Chemotactic.Protein.1..MCP.1.",colnames(m[["Myriad"]]))] = "MCP-1"
colnames(m[["Myriad"]])[grepl("Plasminogen.Activator.Inhibitor.1..PAI.1.",colnames(m[["Myriad"]]))] = "PAI-1"
colnames(m[["Myriad"]])[grepl("Pulmonary.and.Activation.Regulated.Chemokine..PARC.",colnames(m[["Myriad"]]))] = "PARC"
colnames(m[["Myriad"]])[grepl("Tissue.Inhibitor.of.Metalloproteinases.1..TIMP.1.",colnames(m[["Myriad"]]))] = "TIMP-1"
colnames(m[["Myriad"]])[grepl("Tumor.necrosis.factor.receptor.2..TNFR2.",colnames(m[["Myriad"]]))] = "TNFR2"
colnames(m[["Myriad"]])[grepl("Vascular.Cell.Adhesion.Molecule.1..VCAM.1.",colnames(m[["Myriad"]]))] = "VCAM1"

colnames(m[["Metabolon"]])[grepl("Isobar..glucose..fructose..mannose..galactose..allose..altrose..etc.",colnames(m[["Metabolon"]]))] = "Monosaccharides"


top_matrix=do.call(cbind, m)
top_matrix=top_matrix[,-ncol(top_matrix)]


## now run regressions for all top omics against all phenotypes
om_phe_df=matrix(,nrow=ncol(top_matrix), ncol=ncol(phe))
for (o in 1:ncol(top_matrix)){
	for (p in 1:ncol(phe)){

		pval=summary(lm(phe[complete_subjects,p]~top_matrix[,o]))$coefficients[2,4]
		beta=summary(lm(phe[complete_subjects,p]~top_matrix[,o]))$coefficients[2,1]
		om_phe_df[o,p]=-log10(pval)*sign(beta)

	}
}
rownames(om_phe_df)=colnames(top_matrix)
colnames(om_phe_df)=colnames(phe)
om_phe_df=om_phe_df[,which(colnames(om_phe_df) != "age")]  ## rm age , all omics are age-residualized

df=data.frame(om_phe_df)
df$omic=rownames(om_phe_df)
## for plotting - change all values <-25 to -25, and >25 to 25
#exclude=rownames(df)[c(grep("Transcript",rownames(df)), grep("Metabolon",rownames(df)))]
exclude=rownames(df)[c( grep("Metabolon",rownames(df)))]

df_plot=df[which(!(rownames(df) %in% exclude)),]
df_plot[df_plot > 25] <- 25
df_plot[df_plot < (-25)] <- (-25)
df_plot$omic=rownames(df_plot)

## order the omics the same as in heatmap:
order=scan("omics.order.txt",what="character")
order=gsub("^OLINK","Olink",order)
order=gsub("^metab.targ","Biocrates",order)
order=gsub("^prot","AB_array",order)
order=gsub("^myriad","Myriad",order)
order=gsub("^transcript","Transcript",order)

d=dist(om_phe_df[which(!(rownames(om_phe_df) %in% exclude)),])
order=rev(hclust(d)$labels[hclust(d)$order])
write.table(rev(order), file="omics.order.txt",col.names=F,row.names=F,quote=F)

df_plot$omic=factor(df_plot$omic, levels=order)

## order phenotypes:
phe_levels=rev(hclust(dist(t(om_phe_df)))$labels[hclust(dist(t(om_phe_df)))$order])
tmp=df_plot %>%   gather(pheno, value, BMI:Clins)
tmp$pheno=factor(tmp$pheno, levels=phe_levels)

pdf("omics.bubble_chart.pdf", height=11, width=7)
#heatmap(data.matrix(top_matrix)) ## add colSidecolors - archetype
## PCA plot on top_matrix
## tSNE on top matrix?
#heatmap(om_phe_df) 
## bubble chart

  ggplot(tmp, aes(pheno, omic, size = abs(value), col=value)) + 
  geom_point() + scale_colour_gradient2(low = "darkblue", high = "darkred") + theme_bw() +
  xlab("Phenotype") + labs(size="-log10(q)", col="-log10(q)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y=element_blank(),
  	axis.title.y=element_blank()) 

dev.off()