module load gcc/6.2.0
module load intel/redist/2017.2.174
module load intel/perflibs/64/2017_update2
module load R/3.4.0-ICC-MKL

R

library(ggplot2)
library(annotables)

setwd("~/archetypes/omics/")

top_files=list.files(".",pattern="top_*")
top_list=list()
for (i in top_files){
  top_list[[i]]=scan(i,what="character")
}

heat_input=data.frame()

files=list.files(".",pattern="\\.prot.txt")
matrix=matrix(,,ncol=length(top_list$top_prot))
for (i in files){
  group=gsub("\\..+","",i)
  tmp=read.table(i, sep="\t")
  tmp1=tmp[top_list$top_prot,]
  df=data.frame(group=rep(group,length(top_list$top_prot)), omics_var=top_list$top_prot, beta=tmp1$beta,
                q=tmp1$q, q_star=cut(tmp1$q, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", "")), 
                omics_group=rep("prot",length(top_list$top_prot)), label=tmp1$var)
  heat_input=rbind(heat_input,df)  
  matrix=rbind(matrix,tmp1$beta)
  colnames(matrix)=tmp1$var
}
matrix=matrix[-1,]
prot_labels=hclust(dist(t(matrix)))$labels[hclust(dist(t(matrix)))$order]


files=list.files(".",pattern="myriad.txt")
matrix=matrix(,,ncol=length(top_list$top_myriad))

for (i in files){
  group=gsub("\\..+","",i)
  tmp=read.table(i, sep="\t")
  tmp1=tmp[top_list$top_myriad,]
  df=data.frame(group=rep(group,length(top_list$top_myriad)), omics_var=top_list$top_myriad, beta=tmp1$beta,
                q=tmp1$q, q_star=cut(tmp1$q, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", "")), 
                omics_group=rep("myriad",length(top_list$top_myriad)), label=tmp1$var)
  heat_input=rbind(heat_input,df)  
  matrix=rbind(matrix,tmp1$beta)
  colnames(matrix)=tmp1$var
}

matrix=matrix[-1,]
myriad_labels=hclust(dist(t(matrix)))$labels[hclust(dist(t(matrix)))$order]

files=list.files(".",pattern="OLINK.txt")
matrix=matrix(,,ncol=length(top_list$top_olink))

for (i in files){
  group=gsub("\\..+","",i)
  tmp=read.table(i, sep="\t")
  tmp1=tmp[top_list$top_olink,]
  df=data.frame(group=rep(group,length(top_list$top_olink)), omics_var=top_list$top_olink, beta=tmp1$beta,
                q=tmp1$q, q_star=cut(tmp1$q, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", "")), 
                omics_group=rep("olink",length(top_list$top_olink)), label=tmp1$var)
  heat_input=rbind(heat_input,df)  
  matrix=rbind(matrix,tmp1$beta)
  colnames(matrix)=tmp1$var
}

matrix=matrix[-1,]
olink_labels=hclust(dist(t(matrix)))$labels[hclust(dist(t(matrix)))$order]


files=list.files(".",pattern="\\.targ.txt")
matrix=matrix(,,ncol=length(top_list$top_targ_met))
for (i in files){
  group=gsub("\\..+","",i)
  tmp=read.table(i, sep="\t")
  tmp1=tmp[top_list$top_targ_met,]
  df=data.frame(group=rep(group,length(top_list$top_targ_met)), omics_var=top_list$top_targ_met, beta=tmp1$beta,
                q=tmp1$q, q_star=cut(tmp1$q, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", "")), 
                omics_group=rep("metab",length(top_list$top_targ_met)), label=tmp1$var)
  heat_input=rbind(heat_input,df)  
  matrix=rbind(matrix,tmp1$beta)
  colnames(matrix)=tmp1$var
}

matrix=matrix[-1,]
metab_labels=hclust(dist(t(matrix)))$labels[hclust(dist(t(matrix)))$order]

files=list.files(".",pattern="\\.untarg.txt")
matrix=matrix(,,ncol=length(top_list$top_untarg_met))
for (i in files){
  group=gsub("\\..+","",i)
  tmp=read.table(i, sep="\t")
  tmp1=tmp[top_list$top_untarg_met,]
  df=data.frame(group=rep(group,length(top_list$top_untarg_met)), omics_var=top_list$top_untarg_met, beta=tmp1$beta,
                q=tmp1$q, q_star=cut(tmp1$q, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", "")), 
                omics_group=rep("unt_met",length(top_list$top_untarg_met)), label=tmp1$var)
  heat_input=rbind(heat_input,df)  
  matrix=rbind(matrix,tmp1$beta)
  colnames(matrix)=tmp1$var
}

matrix=matrix[-1,]
unt_metab_labels=hclust(dist(t(matrix)))$labels[hclust(dist(t(matrix)))$order]

files=list.files(".",pattern="\\.transcript.txt")
matrix=matrix(,,ncol=length(top_list$top_genes))

for (i in files){
  group=gsub("\\..+","",i)
  tmp=read.table(i, sep="\t")
  tmp1=tmp[top_list$top_genes,]
  df=data.frame(group=rep(group,length(top_list$top_genes)), omics_var=top_list$top_genes, beta=tmp1$beta,
                q=tmp1$q, q_star=cut(tmp1$q, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", "")), 
                omics_group=rep("transcript",length(top_list$top_genes)), label=grch37[tmp1$var,]$symbol)
  heat_input=rbind(heat_input,df)  
  matrix=rbind(matrix,tmp1$beta)
  colnames(matrix)=grch37[tmp1$var,]$symbol
}

matrix=matrix[-1,]
trans_labels=hclust(dist(t(matrix)))$labels[hclust(dist(t(matrix)))$order]

heat_input$label=factor(as.character(heat_input$label))
## rename some labels
heat_input$label=factor(heat_input$label, 
  levels=c(prot_labels,myriad_labels,olink_labels, metab_labels,unt_metab_labels,trans_labels))
levels(heat_input$label)[grep("CXCL1,CXCL3,CXCL2",levels(heat_input$label))]="CXCL1..2..3"
levels(heat_input$label)[grep("Ferritin..FRTN.",levels(heat_input$label))]="Ferritin"
levels(heat_input$label)[grep("Interleukin.18..IL.18.",levels(heat_input$label))]="IL-18"
levels(heat_input$label)[grep("Interleukin.8..IL.8.",levels(heat_input$label))]="IL-8"
levels(heat_input$label)[grep("Isobar..glucose..fructose..mannose..galactose..allose..altrose..etc.",levels(heat_input$label))]="monosaccharides"
levels(heat_input$label)[grep("Plasminogen.Activator.Inhibitor.1..PAI.1.",levels(heat_input$label))]="PAI-1"
levels(heat_input$label)[grep("Tissue.Inhibitor.of.Metalloproteinases.1..TIMP.1.",levels(heat_input$label))]="TIMP-1"
levels(heat_input$label)[grep("Tumor.necrosis.factor.receptor.2..TNFR2.",levels(heat_input$label))]="TNFR2"
levels(heat_input$label)[grep("Macrophage.Inflammatory.Protein.1.beta..MIP.1.beta.",levels(heat_input$label))]="MIP-1-beta"
levels(heat_input$label)[grep("X1.5.anhydroglucitol..1.5.AG.",levels(heat_input$label))]="X1.5.anhydroglucitol"
levels(heat_input$label)[grep("Monocyte.Chemotactic.Protein.1..MCP.1.",levels(heat_input$label))]="MCP-1"
levels(heat_input$label)[grep("Vascular.Cell.Adhesion.Molecule.1..VCAM.1.",levels(heat_input$label))]="VCAM-1"
levels(heat_input$label)[grep("Pulmonary.and.Activation.Regulated.Chemokine..PARC.",levels(heat_input$label))]="PARC"
levels(heat_input$label)[grep("glycerophosphorylcholine..GPC.",levels(heat_input$label))]="glycerophosphorylcholine"

 

 

heat_input$group=factor(heat_input$group, levels=c("Arch4","Arch3","Arch2","Arch1"))
pdf("omics_heatmap.by_group.h1.pdf",width=20,height=4)
p <- ggplot(aes(x=label, y=group, fill=beta), data=heat_input)
fig12 <- p + geom_tile() + scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C") + 
  geom_text(aes(label=q_star), color="black", size=3) + 
  labs(y=NULL, x=NULL, fill="estimate") + theme_bw() + theme(axis.text.x=element_text(angle = -45, hjust = 0)) +
  geom_vline(xintercept=29.5, size=1.5, color="grey50") + 
  geom_vline(xintercept=41.5, size=1.5, color="grey50") + 
  geom_vline(xintercept=69.5, size=1.5, color="grey50") + 
  geom_vline(xintercept=69.5, size=1.5, color="grey50") + 
  geom_vline(xintercept=120.5, size=1.5, color="grey50") 
fig12
dev.off()

###### separate proteomics heatmap:
prot_input=subset(heat_input, grepl("prot|myriad|olink",omics_group))

pdf("proteomics_heatmap.pdf",width=12,height=4)
p <- ggplot(aes(x=label, y=group, fill=beta), data=prot_input)
fig <- p + geom_tile() + scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C") + 
  geom_text(aes(label=q_star), color="black", size=3) + 
  labs(y=NULL, x=NULL, fill="beta") + theme_bw() + theme(axis.text.x=element_text(angle = -45, hjust = 0)) +
  geom_vline(xintercept=29.5, size=1.5, color="grey50") + 
  geom_vline(xintercept=41.5, size=1.5, color="grey50") 
fig
dev.off()

###### separate metabolomics heatmap:
met_input=subset(heat_input, grepl("met",omics_group))
## remove all unannotated metabolites
met_input=subset(met_input, !grepl("X\\.\\.\\.",label))
levels(met_input$label)=gsub("^X","",levels(met_input$label))
pdf("metabolomics_heatmap.pdf",width=10,height=4)
p <- ggplot(aes(x=label, y=group, fill=beta), data=met_input)
fig <- p + geom_tile() + scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C") + 
  geom_text(aes(label=q_star), color="black", size=3) + 
  labs(y=NULL, x=NULL, fill="beta") + theme_bw() + theme(axis.text.x=element_text(angle = -45, hjust = 0)) +
  geom_vline(xintercept=22.5, size=1.5, color="grey50")  
fig
dev.off()

###### separate gene expression heatmap:
genes_input=subset(heat_input, grepl("transcript",omics_group))

pdf("transcriptomics_heatmap.pdf",width=5,height=4)
p <- ggplot(aes(x=label, y=group, fill=beta), data=genes_input)
fig <- p + geom_tile() + scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C") + 
  geom_text(aes(label=q_star), color="black", size=3) + 
  labs(y=NULL, x=NULL, fill="beta") + theme_bw() + theme(axis.text.x=element_text(angle = -45, hjust = 0)) 
fig
dev.off()

