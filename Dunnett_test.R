#!/usr/bin/env R
library(DescTools)

final_df <- data.frame(Tissue=character(), 
                       Class=character(), 
                       Comparison=character(), 
                       diff=numeric(), 
                       lower_ci=numeric(), 
                       upper_ci=numeric(),
                       p_adj=numeric())

summary_data_tidy <- read.csv("data/imported_data/summary_data_tidy.csv")

for (j in unique(summary_data_tidy$Tissue)[1:8]){
  
  for (i in unique(summary_data_tidy$Class)[3:8]){
    df=subset(summary_data_tidy, Treatment!="perf" & Tissue==j & Class==i & Tissue!="B" & Age!="Young")
    df$Treatment<-factor(df$Treatment)
    result=DunnettTest(Frequency ~ Treatment, data=df, control="NT")
    
    attrib=unlist(attributes(result$NT), use.names=FALSE)
    
    final_df = rbind(
      final_df, 
      data.frame(Class=i, 
                 Tissue=j, 
                 Comparison=attrib[3],
                 diff=result$NT[1],
                 lower_ci=result$NT[3],
                 upper_ci=result$NT[5],
                 p_adj=result$NT[7]))
    
    final_df = rbind(
      final_df, 
      data.frame(Class=i, 
                 Tissue=j, 
                 Comparison=attrib[4],
                 diff=result$NT[2],
                 lower_ci=result$NT[4],
                 upper_ci=result$NT[6],
                 p_adj=result$NT[8]))
  }
}

write.csv(x=final_df, file="data/stats/Figure_5_Dunnett_statistics.csv")