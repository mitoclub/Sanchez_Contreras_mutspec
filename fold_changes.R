#!/usr/bin/env R
library(mratios)

final_df <- data.frame(Class=character(), Tissue=character(),Estimate=numeric(), confInt1=numeric(), confInt2=numeric(), pVal=numeric())
summary_data_tidy <- read.csv("data/imported_data/summary_data_tidy.csv")

for (i in unique(summary_data_tidy$Class)[3:8]){
  for (j in unique(summary_data_tidy$Tissue)[1:8]){
    df=subset(summary_data_tidy, Treatment=="NT" & Tissue==j & Class==i & Tissue!='B')
    df$Age<-factor(df$Age)
    ratio_data=ttestratio(Frequency~Age, data=df, alternative="two.sided", rho=1, var.equal=FALSE)
    
    final_df = rbind(
      final_df, 
      data.frame(Class=i, 
                 Tissue=j, 
                 Estimate=ratio_data$estimate[3],
                 confInt1=ratio_data$conf.int[1],
                 confInt2=ratio_data$conf.int[2],
                 qVal=ratio_data$p.value))
  }
}
write.csv(x=final_df, file="data/stats/Figure_2_ratio_statistics.csv")
