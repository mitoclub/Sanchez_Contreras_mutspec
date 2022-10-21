rm(list=ls(all=TRUE))
wd = getwd()
setwd(gsub("/scripts", "/data/mut_files", wd))
filelist <- list.files(pattern=".*\\.dcs\\.mut") 

finalTable = data.frame()
for (f in filelist){
  fileTable = read.table(f, header = TRUE, fill = TRUE)
  fileTable$MouseID = unlist(strsplit(fileTable$sample,"_"))[1]
  fileTable$Treatment = unlist(strsplit(fileTable$sample,"_"))[2]
  fileTable$Tissue = unlist(strsplit(fileTable$sample,"_"))[3]
  finalTable = rbind(finalTable, fileTable)
}
str(finalTable)

length(unique(finalTable$sample)) #172
length(filelist) #172
file = read.csv("../Mouse_aging_mtDNA_summary.csv", header = TRUE)
file = file[,c(1,2)]
file = unique(file)

final = merge(file, finalTable, by=c("MouseID"), all.y = TRUE)
final = final[final$Tissue != "B",]
final = final[final$Treatment != "PERF",]

length(unique(final$MouseID))
table(final$Age)


write.csv(final, file = "../MergedData.csv", row.names = F)


