rm(list=ls())
library(tidyverse)
library(ggplot2)
library(ggstatsplot)
library(survival)
library(stringr)
library(viridis)
library(scales)
library(survival)
library(randomForestSRC)
library(randomSurvivalForest)

load("data/exp/CPTAC_surv_expr.Rdata")
CPTAC_data_exp <- CPTAC_surv_tumor[,-c(1:2)]
#CPTAC_data_exp <- as.data.frame(t(CPTAC_data_exp))

Fe_gene <- rio::import("data/Ferroptosis/Ferroptosis_suppressor.xlsx")
Fe_gene <- Fe_gene[!duplicated(Fe_gene$Gene),]
names(Fe_gene) <- c("Gene Symbol","Group")
Fe_gene <- subset(Fe_gene,select = c("Gene Symbol","Group"))
idx <- intersect(colnames(CPTAC_data_exp),Fe_gene$Gene Symbol)

surv.cli <- CPTAC_surv_tumor[,c(1:2)]
surv.exp <- CPTAC_surv_tumor[,-c(1:2)]
surv.exp <- surv.exp[,colnames(surv.exp)%in%idx]
surv <- cbind(surv.cli,surv.exp)
surv[1:3,1:6]

realdata <- surv

Coxoutput=data.frame()

for(i in colnames(realdata[,3:ncol(realdata)])){
  cox <- coxph(Surv(OS.time, OS) ~ realdata[,i], data = realdata)
  coxSummary = summary(cox)
  Coxoutput=rbind(Coxoutput,cbind(gene=i,HR=coxSummary$coefficients[,"exp(coef)"],
                                  z=coxSummary$coefficients[,"z"],
                                  pvalue=coxSummary$coefficients[,"Pr(>|z|)"],
                                  lower=coxSummary$conf.int[,3],
                                  upper=coxSummary$conf.int[,4]))
}
for(i in c(2:6)){
  Coxoutput[,i] <- as.numeric(as.vector(Coxoutput[,i]))
}
Coxoutput <- arrange(Coxoutput,pvalue)  
Coxoutput <- na.omit(Coxoutput)

#cox plot
plotCoxoutput <-Coxoutput[Coxoutput$pvalue<0.05,]  
#write.csv(plotCoxoutput_6,file = "result/Cox/plotCoxoutput_74.csv",row.names = F)

ggplot(data=plotCoxoutput,aes(x=HR,y=gene,color=pvalue))+
  geom_errorbarh(aes(xmax=upper,xmin=lower),color='black',height=0,size=1.2)+
  geom_point(aes(x=HR,y=gene),size=3.5,shape=18)+   
  geom_vline(xintercept = 1,linetype='dashed',size=1.2)+
  scale_x_continuous(breaks = c(0.75,1,1.30))+
  coord_trans(x='log2')+ 
  ylab("Gene")+  
  xlab("HR")+ 
  labs(color="P value",title ="The univariate cox of ferroptosis suppressor genes" )+
  scale_color_gradient2(low = "#0052D5",mid ="#FFB6C1",high ="#D20103",midpoint = 0.025)+ 
  theme_bw(base_size = 12)+   
  theme(panel.grid =element_blank(), 
        axis.text.x = element_text(face="bold", color="black", size=9),   
        axis.text.y = element_text(face="bold",  color="black", size=9),
        axis.title.x = element_text(face="bold", color="black", size=11),
        axis.title.y = element_text(face="bold",color="black", size=11),
        legend.text= element_text(face="bold", color="black", size=9),
        legend.title = element_text(face="bold", color="black", size=11),
        panel.border = element_rect(colour = 'black',size=1.4),
        plot.title = element_text(face="bold", size=14, hjust=0.5))   

