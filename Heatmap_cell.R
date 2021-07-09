library(biomformat)
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(ggtree)
library(taxize)
library(microbiome) # data analysis and visualisation
library(RColorBrewer) # nice color options
library(ggpubr) # publication quality figures, based on ggplot2
library(picante)
library(pheatmap)
cls1 <- c( "AOM-DSS" = "blue2", "H2O" = "green4")
cls2 <- c("Vendor_1" = "black", "Vendor_2" = "blue3", "Co-housed" = "red3")

otutable <- read.delim2("otu.txt",row.names=1)
colnames(otutable) <- str_remove_all(colnames(otutable),"X")
colnames(otutable) <- str_replace_all(colnames(otutable),"\\.","-")

otutable[] <- lapply(otutable,as.numeric)
OTU = otu_table(as.matrix(otutable), taxa_are_rows = TRUE)
taxonomy <- read.delim2("taxa.txt",row.names=1)
TAX = tax_table(as.matrix(taxonomy))
colnames(TAX) <- c("Kingdom", "Phylum", "Class", "Order", "Family","Genus", "Species")


sample <- read.delim2("sample_metadata.txt",row.names=1)
SAM = sample_data(sample, errorIfNULL = T)
phylob <- phyloseq(OTU,TAX, SAM)
phylob
GP <- phylob
sample_names(GP)
rank_names(GP)
sample_variables(GP)
taxa_names(GP)

GP@sam_data[["Facility"]] <- str_replace_all(GP@sam_data[["Facility"]] , "Facility","Vendor")

write.table(taxa_names(GP),"taxaid.txt",row.names = F,quote = F)

#add the phylogeny tree 
taxize_class <- classification(taxa_names(GP), db = "ncbi")
taxize_tree <- class2tree(taxize_class, check = TRUE)
taxize_phy <- taxize_tree[["phylo"]]
taxa_names(GP) <- taxa_names(taxize_phy)
GP <- merge_phyloseq(GP, taxize_phy)
#quality control 
summary(sample_sums(GP))

##Beta diversity 

Bray_cruits <- ordinate(GP, "PCoA", "bray")
Jaccard <- ordinate(GP, "PCoA", "jaccard")
Unweighted_UniFrac <- ordinate(GP, "PCoA", "unifrac", weighted=F)
Weighted_UniFrac <- ordinate(GP, "PCoA", "unifrac", weighted=T)

betas <- list(Bray_cruits, Jaccard, Unweighted_UniFrac ,  Weighted_UniFrac)
names(betas) <- c("Bray_cruits", "Jaccard", "Unweighted_UniFrac" ,  "Weighted_UniFrac")

pdf("PCoA_beta.pdf")
for (i in 1:length(betas)){
  
  p <- plot_ordination(GP,betas[[i]] , color="Facility",shape = "Treatment") 
  p <- p + ggtitle(names(betas[i])) + geom_point(size = 3)
  p <- p + theme_classic() + scale_color_manual(values = cls2)
  print(p)
  
}
for (i in 1:length(betas)){
  
  p <- plot_ordination(GP,betas[[i]] ,color="Facility",shape = "Treatment") 
  p <- p + ggtitle(names(betas[i])) + geom_point(size = 3)
  p <- p + theme_classic() + scale_color_manual(values = cls2) + stat_ellipse()
  print(p)
  
}

dev.off()
# 

####Heatmap
sample$Facility <- str_replace_all(sample$Facility , "Facility","Vendor")
genus <- read.delim2("genus.txt",row.names=1)
colnames(genus) <- str_remove_all(colnames(genus),"X")
colnames(genus) <- str_replace_all(colnames(genus),"\\.","-")
genus2 <- genus[-c(1:2),]
genus2[] <- lapply(genus2,as.numeric)



col_anno2 <- sample[,3:4]
col_anno2$Facility <- factor(col_anno2$Facility, levels = c("Vendor_1", "Vendor_2" , "Co-housed" ))
col_anno2$Treatment <- factor(col_anno2$Treatment, levels = c( "H2O" ,"AOM-DSS"  ))
my_colour = list( Facility = cls2, Treatment =  cls1 )


distt <- c('correlation', 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski')
clus_method <- c( 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', 'centroid')


pdf("all_possible_heatmaps_genus2.pdf", height = 10, width = 12)

for (i in c(1:7)){
  for (j in c(1:7)) {
    
    print(
      pheatmap(genus2,scale = "row",color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdBu")))(256), cluster_cols = T, show_colnames = T, annotation_col = col_anno2,cellwidth = 9,cellheight =  9,annotation_colors = my_colour,clustering_distance_rows = distt[i],clustering_distance_cols = distt[j] )
          )  
  }
}
dev.off()

pdf("genus2.pdf", height = 10, width = 12)

pheatmap(genus2,scale = "row",color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdBu")))(256), cluster_cols = F, show_colnames = T, annotation_col = col_anno2,cellwidth = 9,cellheight =  9,annotation_colors = my_colour)
dev.off()

