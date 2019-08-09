module load gcc/6.2.0
module load intel/redist/2017.2.174
module load intel/perflibs/64/2017_update2
module load R/3.4.0-ICC-MKL

R

library("archetypes")
library(fossil)

setwd("~/archetypes/")

phe=read.table("../FINAL_baseline_analysis_2019/M0_FINAL.pheno_matrix.wp2_2.flt_rnt_resid.txt",sep="\t",h=T)
m=phe[,c(5:35)]


#### decide how many archetypes to go for - using scree plot
### let's run this for 10 seeds, and see if we always get similar answers
seeds=c("1234","777","2019","98765","13579","24680","13","55555","31415","111")
test_archetypes=c(1:10)

aa=list()
for (s in seeds){
  print(paste0("Seed: ",s))
  set.seed(s)
  aa[[s]]<-stepArchetypes(m, k=test_archetypes, nrep=nrep, method=robustArchetypes)
}


## check how different are the best archetype solutions across the different seed values

for (n_archetypes in test_archetypes){
  best_archetypes=lapply(aa, function(x) bestModel(x[[n_archetypes]]))  
  coefs<-lapply(best_archetypes, function(x) coef(x))
  
  ### at different cut-offs of archetype memberships assign best archetype
  
  adj_rand_df=data.frame()
  for (cut in seq(from=0.01, to=0.99, by=0.01)){
    ## get the dominant archetype - using cut as cut-off for calling
    arch_at_cut=lapply(coefs, function(x) 
      sapply(1:nrow(x), function(y) ifelse(max(x[y,])>=cut, which(x[y,]==max(x[y,])), NA))
    )
    df=data.frame(cut=rep(cut,choose(10,2)), adj_rand= unlist(sapply(1:9, function(sol1) 
      sapply((sol1+1):10, function(sol2) adj.rand.index(arch_at_cut[[sol1]],arch_at_cut[[sol2]]) )))  )
    adj_rand_df=rbind(adj_rand_df, df)
  }
  
  adj_rand_df$cut=factor(adj_rand_df$cut)
  pdf(paste0(n_archetypes,"_archetypes.stability_by_threshold.pdf"))
  p=ggplot(adj_rand_df, aes(x=cut, y=adj_rand)) + geom_boxplot() + theme_bw() +
    scale_x_discrete(labels=c(rep("",9),"0.1",rep("",9),"0.2",rep("",9),"0.3",rep("",9),"0.4",
                              rep("",9),"0.5",rep("",9),"0.6",rep("",9),"0.7",rep("",9),"0.8",
                              rep("",9),"0.9",rep("",9))) +
    xlab("Archetype membership threshold") +  ylab("Adjusted Rand index") + ylim(0,1) +
    ggtitle(paste0(n_archetypes," archetypes"))
  print(p)
  dev.off()
}