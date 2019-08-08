module load gcc/6.2.0
module load intel/redist/2017.2.174
module load intel/perflibs/64/2017_update2
module load R/3.4.0-ICC-MKL

R

options(width=150)

library(dplyr)
library(annotables)
library(gage)
library(reshape)

setwd("~/archetypes/omics/")

immune_dir = "~/immune_cells/"
files<-list.files(path=immune_dir,recursive=TRUE)
immune_cells = list()
files
for(f in 1:length(files)){
  name = gsub(".txt","",files[f])
  df = read.table(paste0(immune_dir,files[f]),header=T)
  immune_cells[[name]] = df
}

a1=read.table("Arch1.transcript.txt")
a2=read.table("Arch2.transcript.txt")
a3=read.table("Arch3.transcript.txt")
a4=read.table("Arch4.transcript.txt")

ens=sapply(strsplit(rownames(a1),split="\\."), function(x) x[1])
HGNC=sapply(ens, function(x) grch37[which(grch37$ensgene==x),]$symbol[1])
a1$ann=HGNC
a2$ann=HGNC
a3$ann=HGNC
a4$ann=HGNC

group_DEGs=list(
  a1_up=as.character(a1[which(a1$q<0.05 & a1$beta>0),]$ann),
  a1_down=as.character(a1[which(a1$q<0.05 & a1$beta<0),]$ann),
  a2_up=as.character(a2[which(a2$q<0.05 & a2$beta>0),]$ann),
  a2_down=as.character(a2[which(a2$q<0.05 & a2$beta<0),]$ann),
  a3_up=as.character(a3[which(a3$q<0.05 & a3$beta>0),]$ann),
  a3_down=as.character(a3[which(a3$q<0.05 & a3$beta<0),]$ann),
  a4_up=as.character(a4[which(a4$q<0.05 & a4$beta>0),]$ann),
  a4_down=as.character(a4[which(a4$q<0.05 & a4$beta<0),]$ann)
)


### define gage function - returns just the enrichment p-value
gage_tab<-function(deseq2.res, genesets){
  deseq2.fc=deseq2.res$logFC
  names(deseq2.fc)=deseq2.res$Gene.symbol
  input = deseq2.fc
  input<-input[!is.na(names(input))]
  gage_res <- gage(input, gsets = genesets, ref = NULL, samp = NULL, same.dir = T, set.size = c(5, 3500), full.table=T)
  return(gage_res$greater[,4])
}

modules_immune = data.frame(gage_tab(immune_cells[[1]], group_DEGs))
colnames(modules_immune) = names(immune_cells)[1]
for (i in 2:length(immune_cells)){
  res=data.frame(gage_tab(immune_cells[[i]], group_DEGs))
  modules_immune = merge(modules_immune,res,by="row.names")
  rownames(modules_immune) = modules_immune$Row.names
  modules_immune = modules_immune[,-1]
  colnames(modules_immune)[i] = names(immune_cells)[i]
}
modules_immune=modules_immune[,c(2:16,1,17:18)]

write.table(modules_immune, file="immune_cell_enrichment.archetype_DEGs.tab",sep="\t",quote=F)


m=matrix(,nrow=4,ncol=ncol(modules_immune))
rownames(m)=c("Arch1","Arch2","Arch3","Arch4")
colnames(m)=colnames(modules_immune)
for (j in 1:ncol(modules_immune)){
  for (i in 1:4){
    if(modules_immune[2*i-1,j]<modules_immune[2*i,j]){
      m[i,j]=log10(modules_immune[2*i-1,j])
    }else{m[i,j]=-log10(modules_immune[2*i,j])}
  }
}
m=m[,-grep("stim",colnames(m))]
colnames(m)=c("B cells","Basophils","cm T cells","DC","em T cells","Eosinophils","immDC","Macrophages","Mast cells","Neutrophils","NK cells","Th1","Th2")
m=data.frame(m)
m$group=rownames(m)
mm=melt(m)
mm$star=cut(-abs(mm$value), breaks=c(-Inf, -3, -2, -1.3, Inf), label=c("***", "**", "*", ""))
mm$variable=factor(mm$variable, levels=rev(levels(mm$variable)))

pdf("immune_cell_enrichment.q0.05.pdf",width=4, height=5)
p <- ggplot(aes(x=group, y=variable, fill=value), data=mm)
p + geom_tile() + scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C") + 
  geom_text(aes(label=star), color="black", size=5) + 
  labs(y=NULL, x=NULL, fill="-log10(p)") + theme_bw() 
dev.off()
