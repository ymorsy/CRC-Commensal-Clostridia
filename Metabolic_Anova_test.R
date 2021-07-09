#load the required libraries 
library(MetaboAnalystR)
library(stringr)
library(ggplot2)
library(tidyverse)
#create the results directories 
dirs <- c( "00_Data scaling","01_PCA","02_PLSDA","03_Correlation","04_Stats","05_Heatmap","06_Enrichment_analysis","07_Pathways","09_Biomarker")
lapply(dirs, dir.create)
file <- list.files(pattern = "*.csv")
file.copy(file , "06_Enrichment_analysis")
file.copy(file , "07_Pathways")
file.copy(file , "09_Biomarker")

mSet<-InitDataObjects("pktable", "stat", FALSE)
mSet<-Read.TextData(mSet, file , "colu", "disc")
mSet$dataSet$cls
mSet$dataSet$cls <- factor(mSet$dataSet$cls, levels = c("CTRL","CC4","NaBu","Pepto"))
mSet$dataSet$orig.cls <- factor(mSet$dataSet$orig.cls, levels = c("CTRL","CC4","NaBu","Pepto"))
mSet$dataSet$cls
mSet$dataSet$orig.cls
mSet<-SanityCheckData(mSet)
mSet<-ReplaceMin(mSet)
mSet<-ContainMissing(mSet) 
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", "NULL", "AutoNorm", ratio=FALSE, ratioNum=20)
#plot the normalized features and samples 
setwd("00_Data scaling")
mSet<-PlotNormSummary(mSet, "Features_normalized", "pdf", 100, width=NA)
mSet<-PlotSampleNormSummary(mSet, "Samples_normalized", "png", 100, width=NA)
setwd("..")

# perform PCA and PLSDA
setwd("01_PCA")

mSet<-PCA.Anal(mSet)
colVec<-c("#008000", "#FF0000","#0000FF","#800080")
#shapeVec<-c(19,1)
mSet<-UpdateGraphSettings(mSet, colVec)

mSet<-PlotPCA2DScore(mSet, "PCA_", "pdf", 100, width=NA, 1,2,0.95,0,0)
mSet<-PlotPCA2DScore(mSet, "PCA2_", "pdf", 100, width=NA, 1,2)
setwd("..")


setwd("02_PLSDA")

mSet<-PLSR.Anal(mSet, reg=TRUE)
mSet<-PlotPLS2DScore(mSet, "PCA2_", "pdf", 100, width=NA, 1,2)
mSet<-PlotPLS2DScore(mSet, "PCA_", "pdf", 100, width=NA, 1,2,0.95,0,0)
setwd("..")
#Perform correlation 
setwd("03_Correlation")

mSet<-FeatureCorrelation(mSet, "pearson", "C40 Butyric acid")
mSet<-PlotCorr(mSet, "NaBu_pattern_", "pdf", 100, width=NA)
setwd("..")

#ANOVA_test 
setwd("04_Stats")
mSet<-ANOVA.Anal(mSet, F, 0.05, "tukey", T)
setwd("..")

#heatmap 
setwd("05_Heatmap")


mSet<-PlotSubHeatMap(mSet, "Top_50_heatmap_", "pdf", 72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", "tanova", 50, "overview", T, T, T, F)
mSet<-PlotSubHeatMap(mSet, "Top_50_heatmap_2_", "pdf", 72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", "tanova", 50, "overview", F, T, T, F)

setwd("..")


#Biomarker analysis
setwd("09_Biomarker")
rm(list = ls())

file <- list.files(pattern = "*.csv")
all_data <- read.csv2(file)

all_data2 <- data.table::transpose(all_data)
colnames(all_data2) <- all_data2[1,]
#all_data2$Samples <- colnames(all_data)
rownames(all_data2) <- colnames(all_data)
all_data2 <- all_data2[-1,]
all_data2 <- rownames_to_column(all_data2,"Samples")

# all_data2$Label <- factor(all_data2$Label, levels = c("CTRL","CC4","Anaerost",  "Eubact", "Faecelb", "Roseb"))
# str(all_data2$Label)

levels = c("CTRL","CC4","NaBu","Pepto")

for (i in 2:length(levels)){
  df <- filter(all_data2, Label == levels[1] |Label == levels[i])
  write.csv(df, paste0(levels[i],".csv"),row.names = F)
}
file.remove(file)

files <- list.files(pattern = "*.csv")

for (f in files){
  file_directory <- str_remove(f , ".csv")
  dir.create(file_directory)
  file.copy(f,file_directory)
  file.remove(f)
  setwd(file_directory)
  mSet<-InitDataObjects("pktable", "roc", FALSE)
  mSet<-Read.TextData(mSet, f, "rowu", "disc")
  mSet<-SanityCheckData(mSet)
  mSet<-ReplaceMin(mSet)
  mSet<-ContainMissing(mSet) 
  mSet<-PreparePrenormData(mSet)
  mSet<-Normalization(mSet, "NULL", "NULL", "AutoNorm", ratio=FALSE, ratioNum=20)
  
  mSet<-SetAnalysisMode(mSet, "univ")
  mSet<-PrepareROCData(mSet)
  mSet<-CalculateFeatureRanking(mSet)
  
  mets <- read.table("metaboanalyst_roc_univ.csv",sep = ",",header = T)
  mets_x <- as.character(mets[,1])

  for (i in 1:20){
    tryCatch({
  print(mets_x[i])
  mSet<-Perform.UnivROC(mSet, mets_x[i], 0, "png", 72, F, T, "closest.topleft", F, "sp", 0.2)
  mSet<-PlotRocUnivBoxPlot(mSet, mets_x[i], 0, "png", 72, T, FALSE)
  
  
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
    
  }
  
  rm(mSet)
  setwd("..")
  
}
setwd("..")

# 
# #pathway analysis
rm(list = ls())
setwd("04_Stats")
anova_results <- read.delim2("anova_posthoc.csv",sep = ",", row.names = 1) 
setwd("..")
setwd("07_Pathways")
file <- list.files(pattern = "*.csv")
all_data <- read.csv2(file,row.names = 1)
row.names(all_data) <- str_remove_all(row.names(all_data),":")
row.names(all_data) <- str_remove_all(row.names(all_data),"\\(")
row.names(all_data) <- str_remove_all(row.names(all_data),"\\)")
row.names(all_data) <- str_remove_all(row.names(all_data),"\\{")
row.names(all_data) <- str_remove_all(row.names(all_data),"\\}")
row.names(all_data) <- str_remove_all(row.names(all_data),"\\[")
row.names(all_data) <- str_remove_all(row.names(all_data),"\\]")


data_anova <- all_data[c("Label",row.names(anova_results)),]
data_anova <- na.omit(rownames_to_column(data_anova, "Metabolite"))

write.csv(data_anova,"data_anova.csv",row.names = F)


mSet<-InitDataObjects("conc", "pathqea", FALSE)
mSet<-Read.TextData(mSet, "data_anova.csv", "colu", "disc")
mSet<-SanityCheckData(mSet)
mSet<-ContainMissing(mSet)
mSet<-ReplaceMin(mSet)
mSet<-CrossReferencing(mSet, "name")
mSet<-CreateMappingResultTable(mSet)
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", "NULL", "AutoNorm", ratio=FALSE, ratioNum=20)

mSet<-SetKEGG.PathLib(mSet, "hsa", "current")
mSet<-SetMetabolomeFilter(mSet, F);
#perform Relative-betweenness Centrality topology
mSet<-CalculateQeaScore(mSet, "rbc", "gt")
mSet<-PlotPathSummary(mSet, F, "path_view_rbc_", "pdf", 200, width=20)
file.rename("pathway_results.csv","Relative-betweenness Centrality topology_result.csv")

# #perform Out-degree Centrality topology
mSet<-CalculateQeaScore(mSet, "dgr", "gt")
mSet<-PlotPathSummary(mSet, F, "path_view_odc_", "pdf", 200, width=20)
file.rename("pathway_results.csv","Out-degree Centrality topology_result.csv")

#draw the barplot for the two different pathway analysis 

pathways_results <-  c("Out-degree Centrality topology_result.csv","Relative-betweenness Centrality topology_result.csv")
draw_barplot <- function(df,name){
  df$ratio <- as.character(paste(df$Hits,"/",df$Total.Cmpd))
  df2 <- df[df$Impact > 0.2,]
  p <- ggplot(df2, aes_string(x = "Impact"  , y = "reorder(X,Impact)")) + geom_bar(stat = "identity" , fill="black" )+geom_col(aes(fill=FDR))+geom_text(aes_string(label=paste0("ratio")),  hjust = 1.5 , color="white", size=3)+theme_bw()+labs(title = name,x = "Impact on pathway", y = "Pathways") 
  return(p)
}
draw_barplot2 <- function(df,name){
  df$ratio <- as.character(paste(df$Hits,"/",df$Total.Cmpd))
  p <- ggplot(df, aes_string(x = "Impact"  , y = "reorder(X,Impact)")) + geom_bar(stat = "identity" , fill="black",width = 0.5 )+geom_text(aes_string(label=paste0("ratio")),  hjust = 1.5 , color="white", size=3)+theme_bw()+labs(title = name,x = "Impact on pathway", y = "Pathways") 
  return(p)
}
for (pathway in pathways_results){
  df <- read.csv(pathway)
  name <- str_remove(pathway,"_result.csv")
  p <- draw_barplot(df,name)
  p2 <- draw_barplot2(df,name)
  
  pdf(paste0(name,".pdf"))
  print(p2)
  print(p)
  dev.off()
}


setwd("..")

#Enrichment analysis
#select only significant data 
rm(list = ls())
setwd("04_Stats")
anova_results <- read.delim2("anova_posthoc.csv",sep = ",", row.names = 1) 
setwd("..")
setwd("06_Enrichment_analysis")
file <- list.files(pattern = "*.csv")
all_data <- read.csv2(file,row.names = 1)
row.names(all_data) <- str_remove_all(row.names(all_data),":")
row.names(all_data) <- str_remove_all(row.names(all_data),"\\(")
row.names(all_data) <- str_remove_all(row.names(all_data),"\\)")
row.names(all_data) <- str_remove_all(row.names(all_data),"\\{")
row.names(all_data) <- str_remove_all(row.names(all_data),"\\}")
row.names(all_data) <- str_remove_all(row.names(all_data),"\\[")
row.names(all_data) <- str_remove_all(row.names(all_data),"\\]")


data_anova <- all_data[c("Label",row.names(anova_results)),]
data_anova <- na.omit(rownames_to_column(data_anova, "Metabolite"))
write.csv(data_anova,"data_anova.csv",row.names = F)


mSet<-InitDataObjects("conc", "msetqea", FALSE)
mSet<-Read.TextData(mSet, "data_anova.csv", "colu", "disc")
mSet<-SanityCheckData(mSet)
mSet<-ContainMissing(mSet)
mSet<-ReplaceMin(mSet)
mSet<-CrossReferencing(mSet, "name")
mSet<-CreateMappingResultTable(mSet)
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", "NULL", "AutoNorm", ratio=FALSE, ratioNum=20)
mSet<-SetMetabolomeFilter(mSet, F)
#mSet<-PlotQEA.MetSet(mSet, "Fatty Acyls", "pdf", 200, width=40)


#Supper class
mSet<-SetCurrentMsetLib(mSet, "super_class", 2)

mSet<-CalculateGlobalTestScore(mSet)

mSet<-PlotQEA.Overview(mSet, "qea_super_", "net", "pdf", 200, width=20)
mSet<-PlotEnrichDotPlot(mSet, "qea", "qea_dot_super_", "pdf", 200, width=10)
file.rename("msea_qea_result.csv","Super_msea_qea_result.csv")

#Main class
mSet<-SetCurrentMsetLib(mSet, "main_class", 2);
mSet<-CalculateGlobalTestScore(mSet)
mSet<-PlotQEA.Overview(mSet, "qea_main_", "net", "pdf", 200, width=20)
mSet<-PlotEnrichDotPlot(mSet, "qea", "qea_dot_main_", "pdf", 200, width=10)
file.rename("msea_qea_result.csv","Main_msea_qea_result.csv")

# #smpdb pathway
# mSet<-SetCurrentMsetLib(mSet, "smpdb_pathway", 2);
# mSet<-CalculateGlobalTestScore(mSet)
# mSet<-PlotQEA.Overview(mSet, "qea_smpdb_", "net", "pdf", 200, width=20)
# mSet<-PlotEnrichDotPlot(mSet, "qea", "qea_dot_smpdb_", "pdf", 200, width=10)
# file.rename("msea_qea_result.csv","Smpdb_msea_qea_result.csv")
# 

# 
# #Kegg pathway
# mSet<-SetCurrentMsetLib(mSet, "kegg_pathway", 2);
# mSet<-CalculateGlobalTestScore(mSet)
# mSet<-PlotQEA.Overview(mSet, "qea_kegg_", "net", "pdf", 200, width=20)
# mSet<-PlotEnrichDotPlot(mSet, "qea", "qea_dot_kegg_", "pdf", 200, width=10)
# file.rename("msea_qea_result.csv","Kegg_msea_qea_result.csv")



#draw the barplot for the two different enrichment analysis 

MSEA_results <-  c("Main_msea_qea_result.csv","Super_msea_qea_result.csv")
draw_barplot <- function(df,name){
  df$ratio <- df$Hits / df$Expected.Q
  df$ratio2 <- as.character(paste(df$Hits,"/",df$Total.Cmpd))
  df2 <- df[df$ratio > 0.5,]
  
  p <- ggplot(df2, aes_string(x = "ratio"  , y = "reorder(X,ratio)")) + geom_bar(stat = "identity" , fill="black"  )+geom_col(aes(fill=FDR))+geom_text(aes_string(label=paste0("ratio2")),  hjust = 1.5 , color="white", size=3)+theme_bw()+labs(title = name,x = "Enrichment score", y = "") 
  return(p)
}
draw_barplot2 <- function(df,name){
  df$ratio <- df$Hits / df$Expected.Q
  df$ratio2 <- as.character(paste(df$Hits,"/",df$Total.Cmpd))
  
  p <- ggplot(df, aes_string(x = "ratio"  , y = "reorder(X,ratio)")) + geom_bar(stat = "identity" , fill="black" ,width = 0.5)+geom_text(aes_string(label=paste0("ratio2")),  hjust = 1.5 , color="white", size=3)+theme_bw()+labs(title = name,x = "Enrichment score", y = "") 
  return(p)
}
for (msea in MSEA_results){
  df <- read.csv(msea)
  name <- str_remove(msea,"_qea__result.csv")
  p <- draw_barplot(df,name)
  p2 <- draw_barplot2(df,name)
  
  pdf(paste0(name,".pdf"),width = 8 , height = 10)
  print(p2)
  print(p)
  
  dev.off()
}

setwd("..")


