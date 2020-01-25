module load gcc/6.2.0
module load intel/redist/2017.2.174
module load intel/perflibs/64/2017_update2
module load R/3.4.0-ICC-MKL

R

library(ggplot2)
library(data.table)
library(GenABEL)    ##  rntransform
options(width=150)

setwd("~/Archetypes_Age/")
phe=read.table("WP2_2.archetype.pheno_matrix.200919.txt",sep="\t",h=T)
archetypes=read.table("archetypes.wp2.2.txt",h=T,sep="\t")
archetypes=cbind(archetypes, model.matrix(~arch_0.6+0, data=archetypes))

### run Wilcox.test for each phenotype and each group against others
rownames(archetypes)=archetypes$studyid

arch_pheno=merge(archetypes,phe,by="row.names")
arch_pheno=arch_pheno[,-1]
arch=colnames(arch_pheno)[7:11]
phenos=colnames(phe)

pheno_heat_input = data.frame()
for (a in arch){
	for (p in phenos){
		df=data.frame(archetype=a, phenotype=p, p=wilcox.test(as.numeric(arch_pheno[,p])~arch_pheno[,a])$p.value,
			est=-wilcox.test(as.numeric(arch_pheno[,p])~arch_pheno[,a], conf.int=T)$estimate)
		pheno_heat_input=rbind(pheno_heat_input,df)
	}
}
pheno_heat_input$archetype=gsub("arch_0.6","",pheno_heat_input$archetype)

m=matrix(, nrow=5,ncol=length(phenos))
for (i in 1:5){
	m[i,]=sapply(phenos, function(x) 
		-log10(wilcox.test(as.numeric(arch_pheno[,x])~arch_pheno[,arch[i]])$p.value) * 
		sign(-wilcox.test(as.numeric(arch_pheno[,x])~arch_pheno[,arch[i]], conf.int=T)$estimate) )
}
rownames(m)=c("A","B","C","D","mix")
colnames(m)=phenos

d=dist(t(m))
h=hclust(d, method="ward.D2")
phe_ordered=h$labels[h$order]
pheno_heat_input$phenotype=factor(pheno_heat_input$phenotype, 
	levels=rev(phe_ordered[c(1,4,5,2,3,19,20,18,21,24,23,6,7,30:32,28,29,25,8:13, 22,26,27,14:17)]))


#pheno_heat_input$phenotype=factor(pheno_heat_input$phenotype, levels=rev(c("Matsuda","Stumvoll","X2.h.OGIS","Clinsb","Clins","fasting.HDL",
#	"PFR1","age","BMI","BSA","WHR","fasting.UCreatinine","fasting.Insulin","mean.ins","fasting.Cpep","basal.isr",
#	"mmtt.120.Insulin","total.isr","fasting.UCPCR","fasting.UCpep","rate.sens","fasting.TG","fasting.ALT","fasting.AST","fasting.Chol","fasting.LDL","fasting.Glucose", "mmtt.120.Glucose","mean.glu","fasting.HbA1c",
#	"fasting.Creatinine","glu.sens")))

pheno_heat_input$star=cut(pheno_heat_input$p, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", "")) 
pheno_heat_input$archetype=factor(pheno_heat_input$archetype, levels=c("A","B","C","D","mix"))

pdf("Fig1B.pheno_heat.by_group.incl_gray.pdf",width=4,height=10)
p <- ggplot(aes(x=archetype, y=phenotype, fill=-log10(p)*sign(est)), data=pheno_heat_input)
p + geom_tile() + scale_fill_gradient2(low="darkblue", mid="white", high="darkred") + 
  geom_text(aes(label=star), color="black", size=3) + 
  labs(y=NULL, x=NULL, fill="-log10P") + theme_bw() #+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
dev.off()

## save table & kruskal beween subgroups 
p=levels(pheno_heat_input$phenotype)
phenos_df=data.frame(phenotype=p,
	Kruskal=sapply(p, function(x) 
		kruskal.test(as.numeric(arch_pheno[,x])~arch_pheno$arch_0.6)$p.value),
	A_p=as.numeric(sapply(p, function(x) subset(pheno_heat_input, archetype=="A" & phenotype==x, p))),
	A_est=as.numeric(sapply(p, function(x) subset(pheno_heat_input, archetype=="A" & phenotype==x, est))),
	B_p=as.numeric(sapply(p, function(x) subset(pheno_heat_input, archetype=="B" & phenotype==x, p))),
	B_est=as.numeric(sapply(p, function(x) subset(pheno_heat_input, archetype=="B" & phenotype==x, est))),
	C_p=as.numeric(sapply(p, function(x) subset(pheno_heat_input, archetype=="C" & phenotype==x, p))),
	C_est=as.numeric(sapply(p, function(x) subset(pheno_heat_input, archetype=="C" & phenotype==x, est))),
	D_p=as.numeric(sapply(p, function(x) subset(pheno_heat_input, archetype=="D" & phenotype==x, p))),
	D_est=as.numeric(sapply(p, function(x) subset(pheno_heat_input, archetype=="D" & phenotype==x, est))),
	mix_p=as.numeric(sapply(p, function(x) subset(pheno_heat_input, archetype=="mix" & phenotype==x, p))),
	mix_est=as.numeric(sapply(p, function(x) subset(pheno_heat_input, archetype=="mix" & phenotype==x, est)))
)
write.table(phenos_df, file="Sup.Table1.phenotype_differences.txt",sep="\t",quote=F)
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
for (a in arch){
	for (p in colnames(MRI_rnt)){
		df=data.frame(archetype=a, phenotype=p, p=wilcox.test(as.numeric(MRI_rnt[,p])~MRI[,a])$p.value,
			est=-wilcox.test(as.numeric(MRI_rnt[,p])~MRI[,a], conf.int=T)$estimate)
		MRI_heat_input=rbind(MRI_heat_input,df)
	}
}
MRI_heat_input$logP=-log10(MRI_heat_input$p) * MRI_heat_input$est/abs(MRI_heat_input$est)
MRI_heat_input$star=cut(MRI_heat_input$p, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", "")) 

p=levels(MRI_heat_input$phenotype)
levels(MRI_heat_input$archetype)=gsub("arch_0.6","",levels(MRI_heat_input$archetype))
MRI_df=data.frame(phenotype=p,
	Kruskal=sapply(p, function(x) 
		kruskal.test(as.numeric(MRI_rnt[,x])~MRI$arch_0.6)$p.value),
	A_p=as.numeric(sapply(p, function(x) subset(MRI_heat_input, archetype=="A" & phenotype==x, p))),
	A_est=as.numeric(sapply(p, function(x) subset(MRI_heat_input, archetype=="A" & phenotype==x, est))),
	B_p=as.numeric(sapply(p, function(x) subset(MRI_heat_input, archetype=="B" & phenotype==x, p))),
	B_est=as.numeric(sapply(p, function(x) subset(MRI_heat_input, archetype=="B" & phenotype==x, est))),
	C_p=as.numeric(sapply(p, function(x) subset(MRI_heat_input, archetype=="C" & phenotype==x, p))),
	C_est=as.numeric(sapply(p, function(x) subset(MRI_heat_input, archetype=="C" & phenotype==x, est))),
	D_p=as.numeric(sapply(p, function(x) subset(MRI_heat_input, archetype=="D" & phenotype==x, p))),
	D_est=as.numeric(sapply(p, function(x) subset(MRI_heat_input, archetype=="D" & phenotype==x, est))),
	mix_p=as.numeric(sapply(p, function(x) subset(MRI_heat_input, archetype=="mix" & phenotype==x, p))),
	mix_est=as.numeric(sapply(p, function(x) subset(MRI_heat_input, archetype=="mix" & phenotype==x, est)))
)

Rob_matrix=read.table("/home/Teams/teamVIP/Analysis/MultiOmics/ClinicalVariables/Matrices/_9/merged_pheno_matrix_short_raw_2-2_direct_19_04_2018_9.txt",h=T,sep="\t")
rownames(Rob_matrix)=Rob_matrix$studyid

Rob_pheno=merge(Rob_matrix,archetypes,by="row.names")
Rob_pheno=Rob_pheno[,-1]
Rob_phenos=colnames(Rob_pheno)[c(23,24,48:75)]
Rob_pheno_rnt=apply(data.matrix(Rob_pheno[,Rob_phenos]),2,rntransform)

Rob_heat_input=data.frame()
for (a in arch){
	for (p in colnames(Rob_pheno_rnt)){
		df=data.frame(archetype=a, phenotype=p, p=wilcox.test(as.numeric(Rob_pheno_rnt[,p])~Rob_pheno[,a])$p.value,
			est=-wilcox.test(as.numeric(Rob_pheno_rnt[,p])~Rob_pheno[,a], conf.int=T)$estimate)
		Rob_heat_input=rbind(Rob_heat_input,df)
	}
}

Rob_heat_input$logP=-log10(Rob_heat_input$p) * Rob_heat_input$est/abs(Rob_heat_input$est)
Rob_heat_input$star=cut(Rob_heat_input$p, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", "")) 

p=levels(Rob_heat_input$phenotype)
levels(Rob_heat_input$archetype)=gsub("arch_0.6","",levels(Rob_heat_input$archetype))
Rob_df=data.frame(phenotype=p,
	Kruskal=sapply(p, function(x) 
		kruskal.test(as.numeric(Rob_pheno_rnt[,x])~Rob_pheno$arch_0.6)$p.value),
	A_p=as.numeric(sapply(p, function(x) subset(Rob_heat_input, archetype=="A" & phenotype==x, p))),
	A_est=as.numeric(sapply(p, function(x) subset(Rob_heat_input, archetype=="A" & phenotype==x, est))),
	B_p=as.numeric(sapply(p, function(x) subset(Rob_heat_input, archetype=="B" & phenotype==x, p))),
	B_est=as.numeric(sapply(p, function(x) subset(Rob_heat_input, archetype=="B" & phenotype==x, est))),
	C_p=as.numeric(sapply(p, function(x) subset(Rob_heat_input, archetype=="C" & phenotype==x, p))),
	C_est=as.numeric(sapply(p, function(x) subset(Rob_heat_input, archetype=="C" & phenotype==x, est))),
	D_p=as.numeric(sapply(p, function(x) subset(Rob_heat_input, archetype=="D" & phenotype==x, p))),
	D_est=as.numeric(sapply(p, function(x) subset(Rob_heat_input, archetype=="D" & phenotype==x, est))),
	mix_p=as.numeric(sapply(p, function(x) subset(Rob_heat_input, archetype=="mix" & phenotype==x, p))),
	mix_est=as.numeric(sapply(p, function(x) subset(Rob_heat_input, archetype=="mix" & phenotype==x, est)))
)
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
for (a in arch){
	for (p in colnames(diet_rnt)){
		df=data.frame(archetype=a, phenotype=p, p=wilcox.test(as.numeric(diet_rnt[,p])~Rob_diet[,a])$p.value,
			est=-wilcox.test(as.numeric(diet_rnt[,p])~Rob_diet[,a], conf.int=T)$estimate)
		diet_heat_input=rbind(diet_heat_input,df)
	}
}
diet_heat_input$logP=-log10(diet_heat_input$p) * diet_heat_input$est/abs(diet_heat_input$est)
diet_heat_input$star=cut(diet_heat_input$p, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", "")) 

p=levels(diet_heat_input$phenotype)
levels(diet_heat_input$archetype)=gsub("arch_0.6","",levels(diet_heat_input$archetype))
diet_df=data.frame(phenotype=p,
	Kruskal=sapply(p, function(x) 
		kruskal.test(as.numeric(diet_rnt[,x])~Rob_diet$arch_0.6)$p.value),
	A_p=as.numeric(sapply(p, function(x) subset(diet_heat_input, archetype=="A" & phenotype==x, p))),
	A_est=as.numeric(sapply(p, function(x) subset(diet_heat_input, archetype=="A" & phenotype==x, est))),
	B_p=as.numeric(sapply(p, function(x) subset(diet_heat_input, archetype=="B" & phenotype==x, p))),
	B_est=as.numeric(sapply(p, function(x) subset(diet_heat_input, archetype=="B" & phenotype==x, est))),
	C_p=as.numeric(sapply(p, function(x) subset(diet_heat_input, archetype=="C" & phenotype==x, p))),
	C_est=as.numeric(sapply(p, function(x) subset(diet_heat_input, archetype=="C" & phenotype==x, est))),
	D_p=as.numeric(sapply(p, function(x) subset(diet_heat_input, archetype=="D" & phenotype==x, p))),
	D_est=as.numeric(sapply(p, function(x) subset(diet_heat_input, archetype=="D" & phenotype==x, est))),
	mix_p=as.numeric(sapply(p, function(x) subset(diet_heat_input, archetype=="mix" & phenotype==x, p))),
	mix_est=as.numeric(sapply(p, function(x) subset(diet_heat_input, archetype=="mix" & phenotype==x, est)))
)
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
for (a in arch){
	for (p in colnames(Adem_rnt)){
		df=data.frame(archetype=a, phenotype=p, p=wilcox.test(as.numeric(Adem_rnt[,p])~Adem_df[,a])$p.value,
			est=-wilcox.test(as.numeric(Adem_rnt[,p])~Adem_df[,a], conf.int=T)$estimate)
		Adem_heat_input=rbind(Adem_heat_input,df)
	}
}

Adem_heat_input$logP=-log10(Adem_heat_input$p) * Adem_heat_input$est/abs(Adem_heat_input$est)
Adem_heat_input$star=cut(Adem_heat_input$p, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", "")) 

p=levels(Adem_heat_input$phenotype)
levels(Adem_heat_input$archetype)=gsub("arch_0.6","",levels(Adem_heat_input$archetype))
GLP_df=data.frame(phenotype=p,
	Kruskal=sapply(p, function(x) 
		kruskal.test(as.numeric(Adem_rnt[,x])~Adem_df$arch_0.6)$p.value),
	A_p=as.numeric(sapply(p, function(x) subset(Adem_heat_input, archetype=="A" & phenotype==x, p))),
	A_est=as.numeric(sapply(p, function(x) subset(Adem_heat_input, archetype=="A" & phenotype==x, est))),
	B_p=as.numeric(sapply(p, function(x) subset(Adem_heat_input, archetype=="B" & phenotype==x, p))),
	B_est=as.numeric(sapply(p, function(x) subset(Adem_heat_input, archetype=="B" & phenotype==x, est))),
	C_p=as.numeric(sapply(p, function(x) subset(Adem_heat_input, archetype=="C" & phenotype==x, p))),
	C_est=as.numeric(sapply(p, function(x) subset(Adem_heat_input, archetype=="C" & phenotype==x, est))),
	D_p=as.numeric(sapply(p, function(x) subset(Adem_heat_input, archetype=="D" & phenotype==x, p))),
	D_est=as.numeric(sapply(p, function(x) subset(Adem_heat_input, archetype=="D" & phenotype==x, est))),
	mix_p=as.numeric(sapply(p, function(x) subset(Adem_heat_input, archetype=="mix" & phenotype==x, p))),
	mix_est=as.numeric(sapply(p, function(x) subset(Adem_heat_input, archetype=="mix" & phenotype==x, est)))
)

library(dplyr)
extra_df=bind_rows(MRI_df, Rob_df, diet_df, GLP_df)
write.table(extra_df, file="Sup.Table2.extra_phenotypes.txt",sep="\t",quote=F)

plot_MRI=as.character(unique(MRI_heat_input[MRI_heat_input$p<0.05,]$phenotype))[c(1:5,8:9)]
plot_Rob=as.character(unique(Rob_heat_input[Rob_heat_input$p<0.05,]$phenotype))[c(9,17)]
plot_Adem=as.character(unique(Adem_heat_input[Adem_heat_input$p<0.05,]$phenotype))

extra_heat_input=rbind(MRI_heat_input[MRI_heat_input$phenotype %in% plot_MRI,],
	Rob_heat_input[Rob_heat_input$phenotype %in% plot_Rob,])
extra_heat_input=rbind(extra_heat_input, Adem_heat_input[Adem_heat_input$phenotype %in% plot_Adem,])
extra_heat_input$phenotype=factor(as.character(extra_heat_input$phenotype))
levels(extra_heat_input$phenotype)=c("BP_diastolic","Glucagon_0min","Glucagon_60min",
	"Liver_fat","Liver_iron","Pancreas_curvature","Pancreas_volume",
	"Subcutaneuous_abd_fat","Trunk_fat","Visceral_fat","Physical_activity",
	"Active_GLP1_0min","Proinsulin_60min","Total_GLP1_0min")
extra_heat_input$phenotype=factor(extra_heat_input$phenotype, levels=c(
	"BP_diastolic","Liver_iron","Liver_fat","Trunk_fat","Subcutaneuous_abd_fat","Visceral_fat",
	"Proinsulin_60min","Total_GLP1_0min",
	"Active_GLP1_0min","Glucagon_60min","Glucagon_0min",
	"Pancreas_volume","Pancreas_curvature","Physical_activity"))

levels(extra_heat_input$archetype)=gsub("arch_0.6","",levels(extra_heat_input$archetype))

pdf("Fig1C.Extra_pheno_heat.by_archetype.pdf",width=4,height=4)
p <- ggplot(aes(x=archetype, y=phenotype, fill=logP), data=extra_heat_input)
p + geom_tile() + scale_fill_gradient2(low="darkblue", mid="white", high="darkred") + 
  geom_text(aes(label=star), color="black", size=3) + 
  labs(y=NULL, x="Archetypes", fill="-log10P") + theme_bw() #+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()
