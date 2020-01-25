### create the phenotype matrices at baseline (0m) and 18months
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

### filter out individuals flagged by Ana/Robert:

include=subset(Robo_matrix,, select=c("studyid","Pass_InclusionExclusion"))

#### merge into final:
m0<-merge(CRF_0m,clin_assays_0m, by="studyid", all=TRUE)
m0<-merge(m0,glyc_mod_0m, by="studyid", all=TRUE)
m0<-merge(m0,include, by="studyid", all=TRUE)

m0=subset(m0,Pass_InclusionExclusion=="Yes")
m0 = m0[,which(!(colnames(m0) %in% c("flag","Pass_InclusionExclusion","CRF_2_100_eligibility.age_at_visit",
                                    "CRF_2_102_general.subject_height","CRF_2_102_general.subject_weight",
                                    "CRF_2_102_general.subject_waist","CRF_2_102_general.subject_hip","basal.glu","basal.ins")))]

### first follow-up (18m)
### CRF - study_id, age, gender,height, waist, weight, hip
### - derive: center, WHR
### get age and BMI from glycemic modeling!

CRF_18m =  subset(CRFs,,select=c("StudyID","CRF_2_100_eligibility.sex","CRF_2_108_general.subject_height","CRF_2_108_general.subject_waist","CRF_2_108_general.subject_weight","CRF_2_108_general.subject_hip"))
CRF_18m$center=factor(as.numeric(gsub("([0-9]+).*", "\\1",CRF_18m$StudyID)))
CRF_18m$WHR=CRF_18m$CRF_2_108_general.subject_waist/CRF_18m$CRF_2_108_general.subject_hip
CRF_18m<-rename(CRF_18m,c(StudyID="studyid"))

### clinical assays - 18 mon
fasting.HbA1c=subset(clinical_assays,visitid==5&samplenum==3&testtitle=="HbA1c",select = c(studyid,result_value))
fasting.Glucose<-subset(clinical_assays,visitid==5&samplenum==4&testtitle=="Glucose",select = c(studyid,result_value))
fasting.Insulin<-subset(clinical_assays,visitid==5&samplenum==4&testtitle=="Insulin",select = c(studyid,result_value))
mmtt.120.Glucose<-subset(clinical_assays,visitid==6&samplenum==5&testtitle=="Glucose",select = c(studyid,result_value))
mmtt.120.Insulin<-subset(clinical_assays,visitid==6&samplenum==5&testtitle=="Insulin",select = c(studyid,result_value))
fasting.Cpep<-subset(clinical_assays,visitid==5&samplenum==4&testtitle=="C-peptide",select = c(studyid,result_value))
fasting.Creatinine<-subset(clinical_assays,visitid==5&samplenum==4&testtitle=="Creatinine",select = c(studyid,result_value))
fasting.UCPCR<-subset(clinical_assays,visitid==5&samplenum==6&testtitle=="UCPCR",select = c(studyid,result_value))
fasting.UCpep<-subset(clinical_assays,visitid==5&samplenum==6&testtitle=="Urine C-peptide",select = c(studyid,result_value))
fasting.UCreatinine<-subset(clinical_assays,visitid==5&samplenum==6&testtitle=="Urine Creatinine",select = c(studyid,result_value))
fasting.HDL<-subset(clinical_assays,visitid==5&samplenum==4&testtitle=="HDL",select = c(studyid,result_value))
fasting.LDL<-subset(clinical_assays,visitid==5&samplenum==4&testtitle=="LDL",select = c(studyid,result_value))
fasting.TG<-subset(clinical_assays,visitid==5&samplenum==4&testtitle=="Triglycerides",select = c(studyid,result_value))
fasting.ALT<-subset(clinical_assays,visitid==5&samplenum==4&testtitle=="ALT",select = c(studyid,result_value))
fasting.AST<-subset(clinical_assays,visitid==5&samplenum==4&testtitle=="AST",select = c(studyid,result_value))
fasting.Chol<-subset(clinical_assays,visitid==5&samplenum==4&testtitle=="Cholesterol",select = c(studyid,result_value))

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

clin_assays_18m=bc.2

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
glyc_mod_18m = subset(glyc_modeling,visit=="M18" & flag==0,select=c("ID","flag","age","BMI","BSA","basal.glu","mean.glu","basal.ins","mean.ins","basal.isr","glu.sens","rate.sens","PFR1","total.isr","X2.h.OGIS","Stumvoll","Matsuda","Clinsb","Clins"))
glyc_mod_18m<-rename(glyc_mod_18m,c(ID="studyid"))

include=subset(Robo_matrix,, select=c("studyid","Pass_InclusionExclusion"))

#### merge into final:
m18<-merge(CRF_18m,clin_assays_18m, by="studyid", all=TRUE)
m18<-merge(m18,glyc_mod_18m, by="studyid", all=TRUE)
m18<-merge(m18,include, by="studyid", all=TRUE)

m18=subset(m18,Pass_InclusionExclusion=="Yes")
m18 = m18[,which(!(colnames(m18) %in% c("flag","Pass_InclusionExclusion",
                                    "CRF_2_108_general.subject_height","CRF_2_108_general.subject_weight",
                                    "CRF_2_108_general.subject_waist","CRF_2_108_general.subject_hip","basal.glu","basal.ins")))]

#Rank normal transformation
nonrnt_list<-c("studyid","CRF_2_100_eligibility.sex","center")
rnt_list<-c("age","BMI","BSA","WHR","fasting.HbA1c","fasting.Glucose",
	"fasting.Insulin","mmtt.120.Glucose","mmtt.120.Insulin","fasting.Cpep","fasting.Creatinine",
	"fasting.UCPCR","fasting.UCpep","fasting.UCreatinine","fasting.HDL","fasting.LDL","fasting.TG",
	"fasting.ALT","fasting.AST","fasting.Chol","mean.glu","mean.ins","basal.isr","glu.sens","rate.sens",
	"PFR1","total.isr","X2.h.OGIS","Stumvoll","Matsuda","Clinsb","Clins")

m18.rnt<-cbind(m18[nonrnt_list],as.data.frame(lapply(m18[rnt_list],rntransform)))

m18.rnt$male=NA
m18.rnt$male[m18.rnt$CRF_2_100_eligibility.sex=="female"]=0
m18.rnt$male[m18.rnt$CRF_2_100_eligibility.sex=="male"]=1

### residualise - for center, and gender:
nonres_list<-c("studyid","center","male","age")
res_list<-c("BMI","BSA","WHR","fasting.HbA1c","fasting.Glucose","fasting.Insulin","mmtt.120.Glucose",
	"mmtt.120.Insulin","fasting.Cpep","fasting.Creatinine","fasting.UCPCR","fasting.UCpep",
	"fasting.UCreatinine","fasting.HDL","fasting.LDL","fasting.TG","fasting.ALT","fasting.AST",
	"fasting.Chol","mean.glu","mean.ins","basal.isr","glu.sens","rate.sens","PFR1","total.isr",
	"X2.h.OGIS","Stumvoll","Matsuda","Clinsb","Clins")

residualise<-function(var){
resid(lm(var~m18.rnt$center+m18.rnt$male, na.action=na.exclude))
}
m18.rnt.resid<-cbind(m18.rnt[nonres_list],as.data.frame(lapply(m18.rnt[res_list],residualise)))

write.table(m18.rnt.resid, file="~/Archetypes_Age/m18.rnt.resid.011019.txt",sep="\t",quote=F)

### 36 months follow-up visit
CRF_36m =  subset(CRFs,,select=c("StudyID","CRF_2_100_eligibility.sex","CRF_2_301_visit.subject_height","CRF_2_301_visit.subject_waist","CRF_2_301_visit.subject_weight","CRF_2_301_visit.subject_hip"))
CRF_36m$center=factor(as.numeric(gsub("([0-9]+).*", "\\1",CRF_36m$StudyID)))
CRF_36m$WHR=CRF_36m$CRF_2_301_visit.subject_waist/CRF_36m$CRF_2_301_visit.subject_hip
CRF_36m<-rename(CRF_36m,c(StudyID="studyid"))


### clinical assays - 36 month - 04/09/2018 - visitid=0
fasting.HbA1c=subset(clinical_assays,visitid==0&samplenum==3&testtitle=="HbA1c",select = c(studyid,result_value))
fasting.Glucose<-subset(clinical_assays,visitid==0&samplenum==6&testtitle=="Glucose",select = c(studyid,result_value))
fasting.Insulin<-subset(clinical_assays,visitid==0&samplenum==6&testtitle=="Insulin",select = c(studyid,result_value))
mmtt.120.Glucose<-subset(clinical_assays,visitid==0&samplenum==30&testtitle=="Glucose",select = c(studyid,result_value))
mmtt.120.Insulin<-subset(clinical_assays,visitid==0&samplenum==30&testtitle=="Insulin",select = c(studyid,result_value))
fasting.Cpep<-subset(clinical_assays,visitid==0&samplenum==6&testtitle=="C-peptide",select = c(studyid,result_value))
fasting.Creatinine<-subset(clinical_assays,visitid==0&samplenum==6&testtitle=="Creatinine",select = c(studyid,result_value))
fasting.UCPCR<-subset(clinical_assays,visitid==0&samplenum==8&testtitle=="UCPCR",select = c(studyid,result_value))
fasting.UCpep<-subset(clinical_assays,visitid==0&samplenum==8&testtitle=="Urine C-peptide",select = c(studyid,result_value))
fasting.UCreatinine<-subset(clinical_assays,visitid==0&samplenum==8&testtitle=="Urine Creatinine",select = c(studyid,result_value))
fasting.HDL<-subset(clinical_assays,visitid==0&samplenum==6&testtitle=="HDL",select = c(studyid,result_value))
fasting.LDL<-subset(clinical_assays,visitid==0&samplenum==6&testtitle=="LDL",select = c(studyid,result_value))
fasting.TG<-subset(clinical_assays,visitid==0&samplenum==6&testtitle=="Triglycerides",select = c(studyid,result_value))
fasting.ALT<-subset(clinical_assays,visitid==0&samplenum==6&testtitle=="ALT",select = c(studyid,result_value))
fasting.AST<-subset(clinical_assays,visitid==0&samplenum==6&testtitle=="AST",select = c(studyid,result_value))
fasting.Chol<-subset(clinical_assays,visitid==0&samplenum==6&testtitle=="Cholesterol",select = c(studyid,result_value))

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

clin_assays_36m=bc.2

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
glyc_mod_36m = subset(glyc_modeling,visit=="M36" & flag==0,select=c("ID","flag","age","BMI","BSA","basal.glu","mean.glu","basal.ins","mean.ins","basal.isr","glu.sens","rate.sens","PFR1","total.isr","X2.h.OGIS","Stumvoll","Matsuda","Clinsb","Clins"))
glyc_mod_36m<-rename(glyc_mod_36m,c(ID="studyid"))

### filter out individuals flagged by Ana/Robert:
include=subset(Robo_matrix,, select=c("studyid","Pass_InclusionExclusion"))

#### merge into final:
m36<-merge(CRF_36m,clin_assays_36m, by="studyid", all=TRUE)
m36<-merge(m36,glyc_mod_36m, by="studyid", all=TRUE)
m36<-merge(m36,include, by="studyid", all=TRUE)

m36=subset(m36,Pass_InclusionExclusion=="Yes")
m36 = m36[,which(!(colnames(m36) %in% c("flag","Pass_InclusionExclusion",
                                    "CRF_2_301_visit.subject_height","CRF_2_301_visit.subject_weight",
                                    "CRF_2_301_visit.subject_waist","CRF_2_301_visit.subject_hip",
                                    "basal.glu","basal.ins")))]
#Rank normal transformation
nonrnt_list<-c("studyid","CRF_2_100_eligibility.sex","center")
rnt_list<-c("age","BMI","BSA","WHR","fasting.HbA1c","fasting.Glucose",
	"fasting.Insulin","mmtt.120.Glucose","mmtt.120.Insulin","fasting.Cpep","fasting.Creatinine",
	"fasting.UCPCR","fasting.UCpep","fasting.UCreatinine","fasting.HDL","fasting.LDL","fasting.TG",
	"fasting.ALT","fasting.AST","fasting.Chol","mean.glu","mean.ins","basal.isr","glu.sens","rate.sens",
	"PFR1","total.isr","X2.h.OGIS","Stumvoll","Matsuda","Clinsb","Clins")

m36.rnt<-cbind(m36[nonrnt_list],as.data.frame(lapply(m36[rnt_list],rntransform)))

m36.rnt$male=NA
m36.rnt$male[m36.rnt$CRF_2_100_eligibility.sex=="female"]=0
m36.rnt$male[m36.rnt$CRF_2_100_eligibility.sex=="male"]=1

### residualise - for center, and gender:
nonres_list<-c("studyid","center","male","age")
res_list<-c("BMI","BSA","WHR","fasting.HbA1c","fasting.Glucose","fasting.Insulin","mmtt.120.Glucose",
	"mmtt.120.Insulin","fasting.Cpep","fasting.Creatinine","fasting.UCPCR","fasting.UCpep",
	"fasting.UCreatinine","fasting.HDL","fasting.LDL","fasting.TG","fasting.ALT","fasting.AST",
	"fasting.Chol","mean.glu","mean.ins","basal.isr","glu.sens","rate.sens","PFR1","total.isr",
	"X2.h.OGIS","Stumvoll","Matsuda","Clinsb","Clins")

residualise<-function(var){
resid(lm(var~m36.rnt$center+m36.rnt$male, na.action=na.exclude))
}
m36.rnt.resid<-cbind(m36.rnt[nonres_list],as.data.frame(lapply(m36.rnt[res_list],residualise)))

write.table(m36.rnt.resid, file="~/Archetypes_Age/m36.rnt.resid.011019.txt",sep="\t",quote=F)


## rbind pheno matrices for M0, M18 and M36
### rnt and residualized together

rownames(m0) = paste0(m0$studyid,"_M0")
rownames(m18) = paste0(m18$studyid,"_M18")
rownames(m36) = paste0(m36$studyid,"_M36")

m0_18 = rbind(m0,m18)
m0_18_36=rbind(m0_18,m36)

write.table(m0_18_36,file="~/Archetypes_Age/m0_18_36.raw.pheno_matrix.txt",sep="\t", quote=F)

#Rank normal transformation
nonrnt_list<-c("studyid","CRF_2_100_eligibility.sex","center")
rnt_list<-c("age","BMI","BSA","WHR","fasting.HbA1c","fasting.Glucose",
	"fasting.Insulin","mmtt.120.Glucose","mmtt.120.Insulin","fasting.Cpep","fasting.Creatinine",
	"fasting.UCPCR","fasting.UCpep","fasting.UCreatinine","fasting.HDL","fasting.LDL","fasting.TG",
	"fasting.ALT","fasting.AST","fasting.Chol","mean.glu","mean.ins","basal.isr","glu.sens","rate.sens",
	"PFR1","total.isr","X2.h.OGIS","Stumvoll","Matsuda","Clinsb","Clins")

m0_18_36.rnt<-cbind(m0_18_36[nonrnt_list],as.data.frame(lapply(m0_18_36[rnt_list],rntransform)))

m0_18_36.rnt$male=NA
m0_18_36.rnt$male[m0_18_36.rnt$CRF_2_100_eligibility.sex=="female"]=0
m0_18_36.rnt$male[m0_18_36.rnt$CRF_2_100_eligibility.sex=="male"]=1

### residualise - for center, and gender:
nonres_list<-c("studyid","center","male","age")
res_list<-c("BMI","BSA","WHR","fasting.HbA1c","fasting.Glucose","fasting.Insulin","mmtt.120.Glucose",
	"mmtt.120.Insulin","fasting.Cpep","fasting.Creatinine","fasting.UCPCR","fasting.UCpep",
	"fasting.UCreatinine","fasting.HDL","fasting.LDL","fasting.TG","fasting.ALT","fasting.AST",
	"fasting.Chol","mean.glu","mean.ins","basal.isr","glu.sens","rate.sens","PFR1","total.isr",
	"X2.h.OGIS","Stumvoll","Matsuda","Clinsb","Clins")

residualise<-function(var){
resid(lm(var~m0_18_36.rnt$center+m0_18_36.rnt$male, na.action=na.exclude))
}
m0_18_36.rnt.resid<-cbind(m0_18_36.rnt[nonres_list],as.data.frame(lapply(m0_18_36.rnt[res_list],residualise)))

write.table(m0_18_36.rnt.resid, file="~/Archetypes_Age/m0_18_36.rnt.resid.011019.txt",sep="\t",quote=F)

### filter to the final set of individuals agreed with Caroline:
tmp=read.table("~/FINAL_baseline_analysis_2019/FINAL.pheno_matrix.wp2_2.m0_m18_m36.txt",h=T,sep="\t")
write.table(m0_18_36.rnt.resid[rownames(tmp),], 
	file="~/Archetypes_Age/WP2_2.m0_18_36.rnt_resid_flt.011019.txt",sep="\t",quote=F)
