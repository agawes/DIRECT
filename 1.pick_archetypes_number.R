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

### let's run this for 10 different seeds, and see if we always get similar answers
seeds=sample(1:1000000, 10)

nrep=10
max_k=10
aa=list()
for (s in seeds){
  print(paste0("Seed: ",s))
  set.seed(s)
  ### which method to use: archetypes, wighted Archetypes, or robustArchetypes
  aa[[s]]<-stepArchetypes(m, k=1:max_k, nrep=nrep, method=robustArchetypes)
}

pdf("pick_archetypes_number.screeplot.pdf")
sapply(seeds, function(s) print(screeplot(aa[[s]], main=s)))
dev.off()

  ### evaluate stability of the clusters:
 ind=list()
 adj.ind=list()

for (s in seeds){
	ind[[s]]=list()
	adj.ind[[s]]=list()
  for (i in 1:max_k){ ## for each k
    ind[[s]][[i]]=rep(NA,nrep)
    adj.ind[[s]][[i]]=rep(NA,nrep)
    
    best=as.numeric(which(rss(aa[[s]])[i,]==min(rss(aa[[s]])[i,]))[1]) ## find best solution for each seed, and each k to use as reference
    if (!is.na(best)){
    m_best=max.col(coef(aa[[s]][[i]][[best]]))
    for (j in 1:nrep){ ## for all other models
      print(paste(s,i,j))
      if (best != j){
        m2=max.col(coef(aa[[s]][[i]][[j]]))
        adj.ind[[s]][[i]][j]=adj.rand.index(m_best,m2)   ### find better way of evaluating stability here!
        ind[[s]][[i]][j]=rand.index(m_best,m2)   
        
      }
    }
    }
  }
}
 
pdf("rand_index.pdf")
for (s in seeds){
  plot(1:max_k, sapply(ind[[s]],function(x) mean(x, na.rm=T)), type="l", ylim=c(0,1), xlab="Archetypes",ylab="Rand index", main=s)
  lines(1:max_k, sapply(adj.ind[[s]],function(x) mean(x, na.rm=T)), col="red")
  legend("topright", c("Standard","Adjusted"), col=c("black","red"),lty=1, bty="n")
}
 dev.off()

#### 4 archetypes consistently performs best!!! #####

save.image("1.pick_archetypes_number.Rdata")