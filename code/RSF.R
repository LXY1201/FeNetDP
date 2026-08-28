#----------------------------------------------------------#
# Random survival forest with independent variable hunting #
#----------------------------------------------------------#
rm(list=ls())
library(randomForestSRC)
library(survival)
library(Hmisc)

#1.CPTAC_data
load("data/exp/CPTAC_surv_expr.Rdata")
load("result/Bootstrap/bootGene.Rdata")


#2. Select gene
surv <- CPTAC_surv_tumor[,-c(1:2)]
surv.cli <- CPTAC_surv_tumor[,c(1:2)]

surv.exp <- surv[, colnames(surv) %in%bootGene$gene]
surv <- cbind(surv.cli,surv.exp)

nodesize <- 20
ntree <- 1000
nrep <- 10
nsplit <- 10
conservative <- "high"    #"low"

res.rsf <- rfsrc(Surv(OS.time, OS) ~ ., surv, 
                 nodesize = nodesize, 
                 proximity = T, 
                 tree.err = T, 
                 forest = T, 
                 ntree = ntree,
                 splitrule = "logrankscore", 
                 importance = TRUE)


res.trc <- res.trcoob <- c()
topvars <- list()

set.seed(1234)

for (j in 1:1000) {
  set.seed(j)
  print(paste("trying for", j, "times"))
  vars<-var.select(object = res.rsf,
                   cause = 1,
                   method =   "vh.vimp", 
                   conservative = conservative, 
                   ntree = ntree,
                   nodesize = nodesize,
                   splitrule = "logrankscore", 
                   nsplit = nsplit, 
                   xvar.wt = NULL,
                   refit = T, 
                   fast = T,
                   na.action = "na.impute", 
                   always.use = NULL, 
                   nrep = nrep,
                   prefit =  list(action = T, 
                                  ntree = ntree,
                                  nodesize = nodesize, 
                                  nsplit = nsplit),
                   verbose = TRUE)
  
  trc <- rcorr.cens(-vars$rfsrc.refit.obj$predicted, 
                    Surv(surv$OS.time, surv$OS))["C Index"]
  trcoob <- rcorr.cens(-vars$rfsrc.refit.obj$predicted.oob, 
                       Surv(surv$OS.time, surv$OS))["C Index"]
  
  res.trc <- rbind(res.trc, trc)
  res.trcoob <- rbind(res.trcoob, trcoob)
  if(length(vars$topvars) == 0) {
    topvars[[j]] <- "[Not Available]"
  } else {
    
    #topvars[[j]] <- vars$topvars   
    a <- as.data.frame(vars$varselect)
    topvars[[j]] <- rownames(a[a$rel.freq>0,,drop = FALSE])
    
  }
}


result <- data.frame(res.trc,res.trcoob,row.names = 1:nrow(res.trc))
colnames(result) <- c("res.trc.cindex","res.trcoob.cindex")
bestresult <- unique(result[result$res.trcoob.cindex == max(result$res.trcoob.cindex),] )
bestvars <- unique(topvars[[as.numeric(rownames(bestresult)[1])]]) 
rsf.res <- cbind(surv[,1:2],surv[,as.character(bestvars)])
bestvars

save(topvars,result,bestresult,bestvars,rsf.res,file = "result/RSF/CPTAC_RSF.Rdata")

