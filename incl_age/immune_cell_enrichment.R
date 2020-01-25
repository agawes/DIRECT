setwd("~/Desktop/DIRECT/ARCHETYPES_CLEAN_FINAL/omics")

library(dplyr)
library(annotables)
library(gage)

immune_dir = "../../immune_cells/"
files<-list.files(path=immune_dir,recursive=TRUE)
immune_cells = list()
files
for(f in 1:length(files)){
    name = gsub(".txt","",files[f])
    df = read.table(paste0(immune_dir,files[f]),header=T)
    immune_cells[[name]] = df
}

trans_files=list.files(path=".", pattern=".transcript.txt", recursive=T)
geneset_trans_results=list()
for(f in 1:length(trans_files)){
  name = gsub(".transcript.txt","",trans_files[f])
  df = read.table(trans_files[f],header=T)
  geneset_trans_results[[paste0(name,"_up")]] = as.character(df[which(df$q<0.05 & df$est >0),]$gene)
  geneset_trans_results[[paste0(name,"_down")]] = as.character(df[which(df$q<0.05 & df$est <0),]$gene)
  
}

### define gage function - returns just the enrichment p-value
gage_tab<-function(deseq2.res, genesets){
  deseq2.fc=deseq2.res$logFC
  names(deseq2.fc)=deseq2.res$Gene.symbol
  input = deseq2.fc
  input<-input[!is.na(names(input))]
  gage_res <- gage(input, gsets = genesets, ref = NULL, samp = NULL, same.dir = T, set.size = c(5, 3500), full.table=T)
  return(gage_res$greater[,4])
}

modules_immune = data.frame(gage_tab(immune_cells[[1]], geneset_trans_results))
colnames(modules_immune) = names(immune_cells)[1]
for (i in 2:length(immune_cells)){
	res=data.frame(gage_tab(immune_cells[[i]], geneset_trans_results))
	modules_immune = merge(modules_immune,res,by="row.names")
	rownames(modules_immune) = modules_immune$Row.names
	modules_immune = modules_immune[,-1]
	colnames(modules_immune)[i] = names(immune_cells)[i]
}

write.table(modules_immune, file="immune_cell_enrichment.wp2.2_archetypes_transcripts.tab",sep="\t",quote=F)
