module load gcc/6.2.0
module load intel/redist/2017.2.174
module load intel/perflibs/64/2017_update2
module load R/3.4.0-ICC-MKL

R
library(reshape)
library(stats)
library(GenABEL) 

CRFs=read.csv(file='/home/Teams/teamVIP/extract/direct_03_09_2018/CRF_Data_WP2.2/All_CRFs_Matrix_WP2.2.csv')
clinical_assays=read.csv(file='/home/Teams/teamVIP/extract/direct_03_09_2018/Clinical_Data/Clinical_assays_Exeter.2.2.csv')
glyc_modeling=read.table('/home/Teams/team09/Andrea/modeling/WP2.2/WP22-modeling-parameters-2018-03-12.txt',sep="\t",h=T)
#medication=read.table("/home/Teams/teamVIP/extract/glucose_lowering/drugs_as_covariates/Cleanto18m2_2tab_delineated12_7_16.raw",sep="\t",h=T)
#Robo_matrix=read.table("/home/Teams/teamVIP/Analysis/MultiOmics/ClinicalVariables/Matrices/_8/merged_pheno_matrix_short_raw_2-2_direct_03_11_2017_8.txt",h=T,sep="\t")
Robo_matrix=read.table("/home/Teams/teamVIP/Analysis/MultiOmics/ClinicalVariables/Matrices/_9/merged_pheno_matrix_short_raw_2-2_direct_19_04_2018_9.txt",h=T,sep="\t")

### baseline (m0)
### CRF - study_id, age, gender,height, waist, weight, hip
### - derive: center, WHR
### BMI - from glycemic modeling
CRF_0m =  subset(CRFs,,select=c("StudyID","CRF_2_100_eligibility.age_at_visit","CRF_2_100_eligibility.sex","CRF_2_102_general.subject_height","CRF_2_102_general.subject_waist","CRF_2_102_general.subject_weight","CRF_2_102_general.subject_hip"))
CRF_0m$center=factor(as.numeric(gsub("([0-9]+).*", "\\1",CRF_0m$StudyID)))
CRF_0m$WHR=CRF_0m$CRF_2_102_general.subject_waist/CRF_0m$CRF_2_102_general.subject_hip
CRF_0m<-rename(CRF_0m,c(StudyID="studyid"))

### clinical assays
### Hba1c, visitnum=2, samplenum=3
#visitid: 1=fasting, 2=MMTT
#if visitid=1 & samplenum: 3=HbA1c, 5=serum 6=plasma 8=urines
#if visitid=2 & samplenum: 1=0min, 2=30min, 3=60min 4=90min 5=120min for each Glucose/Insulin/Cpeptide
#use testtitle to specify "HbA1c" or: 
# serum: "GAD" "GAD autoantibodies" "IA-2" "IA2 autoantibodies"
# plasma: "ALT" "AST" "Cholesterol" "Glucose" "C-peptide" "Creatinine" "HDL" "Insulin" "LDL" "Triglycerides"
# urine: "UCPCR""Urine C-peptide" "Urine Creatinine"
fasting.HbA1c=subset(clinical_assays,visitid==1&samplenum==3&testtitle=="HbA1c",select = c(studyid,result_value))
fasting.Glucose<-subset(clinical_assays,visitid==1&samplenum==6&testtitle=="Glucose",select = c(studyid,result_value))
fasting.Insulin<-subset(clinical_assays,visitid==1&samplenum==6&testtitle=="Insulin",select = c(studyid,result_value))
mmtt.120.Glucose<-subset(clinical_assays,visitid==2&samplenum==5&testtitle=="Glucose",select = c(studyid,result_value))
mmtt.120.Insulin<-subset(clinical_assays,visitid==2&samplenum==5&testtitle=="Insulin",select = c(studyid,result_value))
fasting.Cpep<-subset(clinical_assays,visitid==1&samplenum==6&testtitle=="C-peptide",select = c(studyid,result_value))
fasting.Creatinine<-subset(clinical_assays,visitid==1&samplenum==6&testtitle=="Creatinine",select = c(studyid,result_value))
fasting.UCPCR<-subset(clinical_assays,visitid==1&samplenum==8&testtitle=="UCPCR",select = c(studyid,result_value))
fasting.UCpep<-subset(clinical_assays,visitid==1&samplenum==8&testtitle=="Urine C-peptide",select = c(studyid,result_value))
fasting.UCreatinine<-subset(clinical_assays,visitid==1&samplenum==8&testtitle=="Urine Creatinine",select = c(studyid,result_value))
fasting.HDL<-subset(clinical_assays,visitid==1&samplenum==6&testtitle=="HDL",select = c(studyid,result_value))
fasting.LDL<-subset(clinical_assays,visitid==1&samplenum==6&testtitle=="LDL",select = c(studyid,result_value))
fasting.TG<-subset(clinical_assays,visitid==1&samplenum==6&testtitle=="Triglycerides",select = c(studyid,result_value))
fasting.ALT<-subset(clinical_assays,visitid==1&samplenum==6&testtitle=="ALT",select = c(studyid,result_value))
fasting.AST<-subset(clinical_assays,visitid==1&samplenum==6&testtitle=="AST",select = c(studyid,result_value))
fasting.Chol<-subset(clinical_assays,visitid==1&samplenum==6&testtitle=="Cholesterol",select = c(studyid,result_value))

fasting.HbA1c<-rename(fasting.HbA1c,c(result_value="fasting.HbA1c"))
fasting.Glucose<-rename(fasting.Glucose,c(result_value="fasting.Glucose"))
fasting.Insulin<-rename(fasting.Insulin,c(result_value="fasting.Insulin"))
mmtt.120.Glucose<-rename(mmtt.120.Glucose,c(result_value="mmtt.120.Glucose"))
mmtt.120.Insulin<-rename(mmtt.120.Insulin,c(result_value="mmtt.120.Insulin"))
fasting.Cpep<-rename(fasting.Cpep,c(result_value="fasting.Cpep"))
fasting.Creatinine<-rename(fasting.Creatinine,c(result_value="fasting.Creatinine"))
fasting.UCPCR<-rename(fasting.UCPCR,c(result_value="fasting.UCPCR"))
fasting.UCpep<-rename(fasting.UCpep,c(result_value="fasting.UCpep"))
fasting.UCreatinine<-rename(fasting.UCreatinine,c(result_value="fasting.UCreatinine"))
fasting.HDL<-rename(fasting.HDL,c(result_value="fasting.HDL"))
fasting.LDL<-rename(fasting.LDL,c(result_value="fasting.LDL"))
fasting.TG<-rename(fasting.TG,c(result_value="fasting.TG"))
fasting.ALT<-rename(fasting.ALT,c(result_value="fasting.ALT"))
fasting.AST<-rename(fasting.AST,c(result_value="fasting.AST"))
fasting.Chol<-rename(fasting.Chol,c(result_value="fasting.Chol"))

bc.2<-merge(fasting.HbA1c,fasting.Glucose, by="studyid", all=TRUE)
bc.2<-merge(bc.2,fasting.Insulin, by="studyid", all=TRUE)
bc.2<-merge(bc.2,mmtt.120.Glucose, by="studyid", all=TRUE)
bc.2<-merge(bc.2,mmtt.120.Insulin, by="studyid", all=TRUE)
bc.2<-merge(bc.2,fasting.Cpep, by="studyid", all=TRUE)
bc.2<-merge(bc.2,fasting.Creatinine, by="studyid", all=TRUE)
bc.2<-merge(bc.2,fasting.UCPCR, by="studyid", all=TRUE)
bc.2<-merge(bc.2,fasting.UCpep, by="studyid", all=TRUE)
bc.2<-merge(bc.2,fasting.UCreatinine, by="studyid", all=TRUE)
bc.2<-merge(bc.2,fasting.HDL, by="studyid", all=TRUE)
bc.2<-merge(bc.2,fasting.LDL, by="studyid", all=TRUE)
bc.2<-merge(bc.2,fasting.TG, by="studyid", all=TRUE)
bc.2<-merge(bc.2,fasting.ALT, by="studyid", all=TRUE)
bc.2<-merge(bc.2,fasting.AST, by="studyid", all=TRUE)
bc.2<-merge(bc.2,fasting.Chol, by="studyid", all=TRUE)

clin_assays_0m=bc.2

rm(fasting.Glucose)
rm(fasting.HbA1c)
rm(fasting.Insulin)
rm(mmtt.120.Glucose)
rm(mmtt.120.Insulin)
rm(fasting.Cpep)
rm(fasting.Creatinine)
rm(fasting.UCPCR)
rm(fasting.UCpep)
rm(fasting.UCreatinine)
rm(fasting.HDL)
rm(fasting.LDL)
rm(fasting.TG)
rm(fasting.ALT)
rm(fasting.AST)
rm(fasting.Chol)

#### glycemic modeling:
### only flag=0, no problems according to Andrea
glyc_mod_0m = subset(glyc_modeling,visit=="M0" & flag==0,select=c("ID","flag","age","BMI","BSA","basal.glu","mean.glu","basal.ins","mean.ins","basal.isr","glu.sens","rate.sens","PFR1","total.isr","X2.h.OGIS","Stumvoll","Matsuda","Clinsb","Clins"))
glyc_mod_0m<-rename(glyc_mod_0m,c(ID="studyid"))

## create a file listing individuals removed from analysis because flag!=0 in glycemic modeling
glyc_mod_0m_exclude = subset(glyc_modeling,visit=="M0" & flag!=0,select=c("ID","flag"))
#write.table(glyc_mod_0m_exclude, file="M0.glycemic_modeling_flag.individuals_excluded.txt", row.names=F, quote=F,sep="\t")

#### medication:
#med_0m = subset(medication, Metformin0m==1,select=c("studyid","Metformin0m"))

### filter out individuals flagged by Ana/Robert:

include=subset(Robo_matrix,Pass_InclusionExclusion=="Yes", select=c("studyid","Pass_InclusionExclusion"))

## create a file listing individuals removed from analysis because they were flagged in Robert's matrices
criteria_0m_exclude = subset(Robo_matrix,Pass_InclusionExclusion=="No", select=c("studyid","Pass_InclusionExclusion"))
#write.table(criteria_0m_exclude, file="M0.Pass_InclusionExclusion.individuals_excluded.txt", row.names=F, quote=F,sep="\t")

#### merge into final:
m0<-merge(CRF_0m,clin_assays_0m, by="studyid", all=TRUE)
m0<-merge(m0,glyc_mod_0m, by="studyid", all=TRUE)
#m0<-merge(m0,med_0m, by="studyid", all=TRUE)
m0<-merge(m0,include, by="studyid", all=TRUE)

## create a file listing individuals removed at this stage because they were not in the "Pass_InclusionExclusion" list"
missing_Pass_InclusionExclusion=subset(m0,is.na(Pass_InclusionExclusion))
#write.table(missing_Pass_InclusionExclusion, file="M0.miss_Pass_InclusionExclusion.individuals_excluded.txt", row.names=F, quote=F,sep="\t")


m0=subset(m0,Pass_InclusionExclusion=="Yes")
m0 = m0[,which(!(colnames(m0) %in% c("flag","Pass_InclusionExclusion","CRF_2_100_eligibility.age_at_visit",
                                    "CRF_2_102_general.subject_height","CRF_2_102_general.subject_weight",
                                    "CRF_2_102_general.subject_waist","CRF_2_102_general.subject_hip","basal.glu","basal.ins")))]

#Rank normal transformation
nonrnt_list<-c("studyid","CRF_2_100_eligibility.sex","center")
rnt_list<-c("age","BMI","BSA","WHR","fasting.HbA1c","fasting.Glucose","fasting.Insulin","mmtt.120.Glucose","mmtt.120.Insulin",
            "fasting.Cpep","fasting.Creatinine","fasting.UCPCR","fasting.UCpep","fasting.UCreatinine","fasting.HDL","fasting.LDL",
            "fasting.TG","fasting.ALT","fasting.AST","fasting.Chol","mean.glu","mean.ins","basal.isr","glu.sens","rate.sens","PFR1",
            "total.isr","X2.h.OGIS","Stumvoll","Matsuda","Clinsb","Clins")

m0.rnt<-cbind(m0[nonrnt_list],as.data.frame(lapply(m0[rnt_list],rntransform)))

m0.rnt$male=NA
m0.rnt$male[m0.rnt$CRF_2_100_eligibility.sex=="female"]=0
m0.rnt$male[m0.rnt$CRF_2_100_eligibility.sex=="male"]=1

### residualise - for center, age and gender:
nonres_list<-c("studyid","center","male","age")
res_list<-c("BMI","BSA","WHR","fasting.HbA1c","fasting.Glucose","fasting.Insulin","mmtt.120.Glucose","mmtt.120.Insulin",
            "fasting.Cpep","fasting.Creatinine","fasting.UCPCR","fasting.UCpep","fasting.UCreatinine","fasting.HDL","fasting.LDL",
            "fasting.TG","fasting.ALT","fasting.AST","fasting.Chol","mean.glu","mean.ins","basal.isr","glu.sens","rate.sens","PFR1",
            "total.isr","X2.h.OGIS","Stumvoll","Matsuda","Clinsb","Clins")

# residualise<-function(var){
# resid(lm(var~m0.rnt$center+m0.rnt$age+m0.rnt$CRF_2_100_eligibility.sex, na.action=na.exclude))
# }
# m0.rnt.resid<-cbind(m0.rnt[nonres_list],as.data.frame(lapply(m0.rnt[res_list],residualise)))

write.table(m0, file="M0.raw.221118.txt",sep="\t",quote=F)
# write.table(m0.rnt.resid, file="M0.rnt.resid.221118.txt",sep="\t",quote=F)

### add residualization only by sex + center; age is only included as covariate (to compare with Ahlqvist clusters)
residualise<-function(var){
resid(lm(var~m0.rnt$center+m0.rnt$CRF_2_100_eligibility.sex, na.action=na.exclude))
}
m0.rnt.resid<-cbind(m0.rnt[nonres_list],as.data.frame(lapply(m0.rnt[res_list],residualise)))

write.table(m0.rnt.resid, file="M0.rnt.resid.age_as_covariate.190919.txt",sep="\t",quote=F)

# ### add also GLP1 variables:
# glp1=read.csv("~/pheno_tables/from_Adem/WP2.2_GLP1_22_feb.csv",h=T)
# rownames(glp1)=glp1$studyid
# m0=cbind(m0, glp1[as.character(m0$studyid),5:7])

# glucagon=read.table("~/pheno_tables/from_Adem/WP2.2_glucagon_22_feb.csv",h=T,sep="\t")
# rownames(glucagon)=glucagon$studyid
# m0=cbind(m0, glucagon[as.character(m0$studyid),4:5])

# #Rank normal transformation
# nonrnt_list<-c("studyid","CRF_2_100_eligibility.sex","center")
# rnt_list<-c("age","BMI","BSA","WHR","fasting.HbA1c",
# 	"fasting.Glucose","fasting.Insulin","mmtt.120.Glucose","mmtt.120.Insulin","fasting.Cpep",
# 	"fasting.Creatinine","fasting.UCPCR","fasting.UCpep","fasting.UCreatinine","fasting.HDL",
# 	"fasting.LDL","fasting.TG","fasting.ALT","fasting.AST","fasting.Chol","mean.glu","mean.ins",
# 	"basal.isr","glu.sens","rate.sens","PFR1","total.isr","X2.h.OGIS","Stumvoll","Matsuda","Clinsb",
# 	"Clins", "active_glp1_conc_0min","total_glp1_conc_0min","total_glp1_conc_60min", 
# 	"Glucagon_conc_0min_pg_ml_corrected", "Glucagon_conc_60min_pg_ml_corrected")

# m0.rnt<-cbind(m0[nonrnt_list],as.data.frame(lapply(m0[rnt_list],rntransform)))

# m0.rnt$male=NA
# m0.rnt$male[m0.rnt$CRF_2_100_eligibility.sex=="female"]=0
# m0.rnt$male[m0.rnt$CRF_2_100_eligibility.sex=="male"]=1

# ### residualise - for center, and gender:
# nonres_list<-c("studyid","center","male","age")
# res_list<-c("BMI","BSA","WHR","fasting.HbA1c",
# 	"fasting.Glucose","fasting.Insulin","mmtt.120.Glucose","mmtt.120.Insulin","fasting.Cpep",
# 	"fasting.Creatinine","fasting.UCPCR","fasting.UCpep","fasting.UCreatinine","fasting.HDL",
# 	"fasting.LDL","fasting.TG","fasting.ALT","fasting.AST","fasting.Chol","mean.glu","mean.ins",
# 	"basal.isr","glu.sens","rate.sens","PFR1","total.isr","X2.h.OGIS","Stumvoll","Matsuda","Clinsb",
# 	"Clins","active_glp1_conc_0min","total_glp1_conc_0min","total_glp1_conc_60min", 
# 	"Glucagon_conc_0min_pg_ml_corrected", "Glucagon_conc_60min_pg_ml_corrected")

# residualise<-function(var){
# resid(lm(var~m0.rnt$center+m0.rnt$CRF_2_100_eligibility.sex, na.action=na.exclude))
# }
# m0.rnt.resid<-cbind(m0.rnt[nonres_list],as.data.frame(lapply(m0.rnt[res_list],residualise)))

# write.table(m0.rnt.resid, file="m0.rnt.resid.glp_glucagon.age_as_covariate.200919.txt",sep="\t",quote=F)

