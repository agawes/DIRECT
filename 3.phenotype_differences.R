module load gcc/6.2.0
module load intel/redist/2017.2.174
module load intel/perflibs/64/2017_update2
module load R/3.4.0-ICC-MKL

R

setwd("~/archetypes/")
load("WP2.2_archetypes.Rdata")

library(ggplot2)
library(data.table)
library(GenABEL)    ##  rntransform

### what are the differences in phenotypes between archetypes
pheno_heat_input=data.frame()
for (a in 1:4){
  df=data.frame(archetype=rep(a,ncol(m)), phenotype=colnames(m), 
                p=sapply(1:ncol(m), function(x) summary(lm(coef[,a] ~ m[,x]))$coefficients[8]),
                beta=sapply(1:ncol(m), function(x) summary(lm(coef[,a] ~ m[,x]))$coefficients[2]) )
  pheno_heat_input=rbind(pheno_heat_input, df)
}

pheno_heat_input$logP=-log10(pheno_heat_input$p) * pheno_heat_input$beta/abs(pheno_heat_input$beta)
pheno_heat_input$logP_plot=pheno_heat_input$logP
pheno_heat_input$logP_plot[which(pheno_heat_input$logP_plot<(-50))]=-50
pheno_heat_input$logP_plot[which(pheno_heat_input$logP_plot>50)]=50

pheno_heat_input$star=cut(pheno_heat_input$p, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", "")) 

pheno_heat_input$phenotype=factor(pheno_heat_input$phenotype, levels=rev(c("X2.h.OGIS","Matsuda","Stumvoll",
	"Clinsb","Clins","fasting.HDL","fasting.Chol","fasting.LDL","PFR1","BMI","BSA","WHR","fasting.UCreatinine",
	"fasting.Creatinine","fasting.UCPCR","fasting.UCpep","fasting.Insulin","mean.ins","fasting.Cpep","basal.isr",
	"mmtt.120.Insulin","total.isr","fasting.TG","fasting.ALT","fasting.AST","rate.sens","glu.sens",
	"fasting.Glucose", "mmtt.120.Glucose","mean.glu","fasting.HbA1c")))

pdf("Fig1B.pheno_heat.by_archetype.pdf",width=4,height=10)
p=ggplot(aes(x=archetype, y=phenotype, fill=logP_plot), data=pheno_heat_input)
p + geom_tile() + scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C") + 
  geom_text(aes(label=star), color="black", size=3) +
  labs(y=NULL, x="Archetypes", fill="-log10P") + theme_bw()
dev.off()

###### differences in phenotypes not included in clustering
##### phenotypes not measured in whole cohort (MRI, GLP, diet, physical activity)
MRI=fread("~/pheno_tables/DIRECT_MRI_RESULTS_Master.csv")
MRI=data.frame(MRI)
MRI$LIVER.FAT....=as.numeric(MRI$LIVER.FAT....)
MRI$LIVER.IRON..mg.g.tissue.=as.numeric(MRI$LIVER.IRON..mg.g.tissue.)
MRI$PANCREAS.FAT....=as.numeric(MRI$PANCREAS.FAT....)
MRI$PANCREAS.IRON..mg.g.tissue.=as.numeric(MRI$PANCREAS.IRON..mg.g.tissue.)

MRI= subset(MRI, grepl("WP2.2", MRI$FOLDER.GROUP))
MRI= subset(MRI, !grepl("18", MRI$FOLDER.GROUP))
# to keep only baseline 
names(MRI)[names(MRI) == "ID"] <- "studyid"
MRI$studyid <- unlist(lapply(strsplit(as.character(MRI$studyid), split = "_", fixed = TRUE), "[", 1))
MRI <- MRI[,-c(1:2)]
rownames(MRI)=MRI$studyid
MRI=merge(MRI,archetypes,by="row.names")
MRI=MRI[,-1]

### rntransform so the estimates are not crazy
rnt_list=colnames(MRI)[2:12]
MRI_rnt=apply(data.matrix(MRI[,rnt_list]),2,rntransform)

MRI_heat_input=data.frame()
for (a in 1:4){
  df=data.frame(archetype=rep(a,ncol(MRI_rnt)), phenotype=colnames(MRI_rnt), 
                p=sapply(1:ncol(MRI_rnt), function(x) summary(lm(MRI[,12+a] ~ MRI_rnt[,x]))$coefficients[8]),
                beta=sapply(1:ncol(MRI_rnt), function(x) summary(lm(MRI[,12+a] ~ MRI_rnt[,x]))$coefficients[2]) )
  MRI_heat_input=rbind(MRI_heat_input, df)
}

MRI_heat_input$logP=-log10(MRI_heat_input$p) * MRI_heat_input$beta/abs(MRI_heat_input$beta)
MRI_heat_input$star=cut(MRI_heat_input$p, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", "")) 



Rob_matrix=read.table("/home/Teams/teamVIP/Analysis/MultiOmics/ClinicalVariables/Matrices/_9/merged_pheno_matrix_short_raw_2-2_direct_19_04_2018_9.txt",h=T,sep="\t")
rownames(Rob_matrix)=Rob_matrix$studyid

Rob_pheno=merge(Rob_matrix,archetypes,by="row.names")
Rob_pheno=Rob_pheno[,-1]
Rob_phenos=colnames(Rob_pheno)[c(23,24,48:75)]
Rob_pheno_rnt=apply(data.matrix(Rob_pheno[,Rob_phenos]),2,rntransform)

Rob_heat_input=data.frame()
for (a in 1:4){
  df=data.frame(archetype=rep(a,ncol(Rob_pheno_rnt)), phenotype=colnames(Rob_pheno_rnt), 
                p=sapply(1:ncol(Rob_pheno_rnt), function(x) summary(lm(Rob_pheno[,83+a] ~ Rob_pheno_rnt[,x]))$coefficients[8]),
                beta=sapply(1:ncol(Rob_pheno_rnt), function(x) summary(lm(Rob_pheno[,83+a] ~ Rob_pheno_rnt[,x]))$coefficients[2]) )
  Rob_heat_input=rbind(Rob_heat_input, df)
}

Rob_heat_input$logP=-log10(Rob_heat_input$p) * Rob_heat_input$beta/abs(Rob_heat_input$beta)
Rob_heat_input$star=cut(Rob_heat_input$p, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", "")) 


### Robert's matrix - long - contains the more detailed diet info
Rob_matrix_l=read.table("/home/Teams/teamVIP/Analysis/MultiOmics/ClinicalVariables/Matrices/_9/merged_pheno_matrix_full_raw_2-2_direct_19_04_2018_9.txt",h=T,sep="\t")
## keep only the diet variables, and studyid
Rob_matrix_l=Rob_matrix_l[,c(1,296:345)]
rownames(Rob_matrix_l) = Rob_matrix_l$studyid

Rob_diet=merge(Rob_matrix_l,archetypes,by="row.names")
Rob_diet=Rob_diet[,-1]
Rob_diet_phenos=colnames(Rob_diet)[c(2:51)]
diet_rnt=apply(data.matrix(Rob_diet[,Rob_diet_phenos]),2,rntransform)

diet_heat_input=data.frame()
for (a in 1:4){
  df=data.frame(archetype=rep(a,ncol(diet_rnt)), phenotype=colnames(diet_rnt), 
                p=sapply(1:ncol(diet_rnt), function(x) summary(lm(Rob_diet[,51+a] ~ diet_rnt[,x]))$coefficients[8]),
                beta=sapply(1:ncol(diet_rnt), function(x) summary(lm(Rob_diet[,51+a] ~ diet_rnt[,x]))$coefficients[2]) )
  diet_heat_input=rbind(diet_heat_input, df)
}

diet_heat_input$logP=-log10(diet_heat_input$p) * diet_heat_input$beta/abs(diet_heat_input$beta)
diet_heat_input$star=cut(diet_heat_input$p, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", "")) 


#### GLP, proinsulin data from Adem
glp1=read.csv("~/pheno_tables/from_Adem/WP2.2_GLP1_22_feb.csv",h=T)
glucagon=read.table("~/pheno_tables/from_Adem/WP2.2_glucagon_22_feb.csv",h=T,sep="\t")
#gluc120=read.csv("~/pheno_tables/from_Adem/WP22_120MIN_GLUCAGON.csv",h=T) - too many repeated 
proins=read.csv("~/pheno_tables/from_Adem/wp22_proinsulin.csv",h=T)

Adem_df=merge(glp1,glucagon,by="studyid",all=T)
#Adem_df=merge(Adem_df,gluc120, by="studyid",all=T)
Adem_df=merge(Adem_df,proins, by="studyid",all=T)
rownames(Adem_df) = Adem_df$studyid

Adem_df=merge(Adem_df,archetypes,by="row.names")
Adem_df=Adem_df[,-1]
phenos_Adem=colnames(Adem_df)[5:12]
Adem_rnt=apply(data.matrix(Adem_df[,phenos_Adem]),2,rntransform)

Adem_heat_input=data.frame()
for (a in 1:4){
  df=data.frame(archetype=rep(a,ncol(Adem_rnt)), phenotype=colnames(Adem_rnt), 
                p=sapply(1:ncol(Adem_rnt), function(x) summary(lm(Adem_df[,12+a] ~ Adem_rnt[,x]))$coefficients[8]),
                beta=sapply(1:ncol(Adem_rnt), function(x) summary(lm(Adem_df[,12+a] ~ Adem_rnt[,x]))$coefficients[2]) )
  Adem_heat_input=rbind(Adem_heat_input, df)
}

Adem_heat_input$logP=-log10(Adem_heat_input$p) * Adem_heat_input$beta/abs(Adem_heat_input$beta)
Adem_heat_input$star=cut(Adem_heat_input$p, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", "")) 


plot_phenos=c(as.character(unique(MRI_heat_input[MRI_heat_input$p<0.01,]$phenotype))[1:4],
	as.character(unique(Rob_heat_input[Rob_heat_input$p<0.01,]$phenotype))[c(8,14,15)],
	as.character(unique(Adem_heat_input[Adem_heat_input$p<0.01,]$phenotype))
 )
plot_MRI=as.character(unique(MRI_heat_input[MRI_heat_input$p<0.01,]$phenotype))[c(1:4)]
plot_Rob=as.character(unique(Rob_heat_input[Rob_heat_input$p<0.01,]$phenotype))[c(8,14,15)]
plot_Adem=as.character(unique(Adem_heat_input[Adem_heat_input$p<0.01,]$phenotype))

extra_heat_input=rbind(MRI_heat_input[MRI_heat_input$phenotype %in% plot_MRI,],
	Rob_heat_input[Rob_heat_input$phenotype %in% plot_Rob,])
extra_heat_input=rbind(extra_heat_input, Adem_heat_input[Adem_heat_input$phenotype %in% plot_Adem,])
extra_heat_input$phenotype=factor(as.character(extra_heat_input$phenotype))
levels(extra_heat_input$phenotype)=c("BP_diastolic","BP_systolic","Glucagon_0min","Glucagon_60min",
	"Liver_fat","Subcutaneuous_abd_fat","Trunk_fat","Visceral_fat","Physical_activity",
	"Proinsulin_60min","GLP1_0min","GLP1_60min")
extra_heat_input$phenotype=factor(extra_heat_input$phenotype, levels=c("GLP1_60min","GLP1_0min","Glucagon_60min","Glucagon_0min",
	"Proinsulin_60min","Subcutaneuous_abd_fat","Trunk_fat","Visceral_fat","Liver_iron","Liver_fat","BP_systolic","BP_diastolic",
	"Physical_activity"))

pdf("Fig1C.Extra_pheno_heat.by_archetype.pdf",width=4,height=4)
p <- ggplot(aes(x=archetype, y=phenotype, fill=logP), data=extra_heat_input)
p + geom_tile() + scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C") + 
  geom_text(aes(label=star), color="black", size=3) + 
  labs(y=NULL, x="Archetypes", fill="-log10P") + theme_bw() 
dev.off()
