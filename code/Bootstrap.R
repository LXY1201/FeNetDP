#------------------------------------------------------------------#
# bootstrap approach to enhance the robustness of prognostic value #
#------------------------------------------------------------------#

library(survival)

#1.CPTAC_data
load("data/exp/CPTAC_surv_expr.Rdata")

load("result/Cox/Ferroptosis_suppressor.Rdata")

CPTAC_exp <- CPTAC_surv_tumor[,-c(1:2)]
CPTAC_clincial <- CPTAC_surv_tumor[,c(1:2)]

idx <- intersect(Coxoutput.km$gene,colnames(CPTAC_exp))

CPTAC_exp <- CPTAC_exp[,colnames(CPTAC_exp)%in%idx]
CPTAC_exp <- cbind(CPTAC_clincial,CPTAC_exp)

#2.筛基因

surv <- CPTAC_exp

set.seed(1234)

outTab <- NULL
for(i in 3:ncol(surv)){ # survival information (OS in this case)
  set.seed(i)
  # customized function
  display.progress = function (index, totalN, breakN=20) {
    if ( index %% ceiling(totalN/breakN)  ==0  ) {
      cat(paste(round(index*100/totalN), "% ", sep=""))
    }
  }    
  
  display.progress(index = i, totalN = ncol(surv)) # show running progression
  gene <- colnames(surv)[i]
  Mboot <- replicate(1000, expr = { # bootstrap for 1,000 times
    indices <- sample(rownames(surv), size = nrow(surv) * 0.8, replace = F) # extract 80% samples at each bootsratp
    data <- surv[indices,]
    fmla1 <- as.formula(Surv(data[,"OS.time"],data[,"OS"]) ~ data[,gene])
    mycox <- coxph(fmla1,data = data)
    coxResult <- summary(mycox)
    P <- coxResult$coefficients[,"Pr(>|z|)"]
  }
  )
  times <- length(Mboot[which(Mboot < 0.05)])
  outTab <- rbind(outTab,
                  cbind(gene = gene,
                        times = times))
}
outTab <- as.data.frame(outTab)
save(outTab,file = "result/Bootstrap/outTab.Rdata")

