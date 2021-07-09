library(biomformat)
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(ggtree)
library(taxize)
library(microbiome) # data analysis and visualisation
# library(microbiomeutilities) # some utility tools 
library(RColorBrewer) # nice color options
library(ggpubr) # publication quality figures, based on ggplot2
# library(DT) # interactive tables in html and markdown
# library(data.table) # alternative to data.frame
# library(dplyr) # data handling  
library(picante)

mydata<-read_biom("all_samples.biom")
otutable<-as.data.frame(as.matrix(biom_data(mydata)))
taxonomy <- observation_metadata(mydata)
OTU = otu_table(as.matrix(otutable), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(taxonomy))
colnames(TAX) <- c("Kingdom", "Phylum", "Class", "Order", "Family","Genus", "Species")
sample <- read.delim2("sample_metadata.txt",row.names=1)[1:4]
SAM = sample_data(sample, errorIfNULL = T)
phylob <- phyloseq(OTU,TAX, SAM)
phylob
GP <- phylob
sample_names(GP)
rank_names(GP)
sample_variables(GP)
taxa_names(GP)
write.table(taxa_names(GP),"taxaid.txt",row.names = F,quote = F)

#add the phylogeny tree 
taxize_class <- classification(taxa_names(GP), db = "ncbi")
taxize_tree <- class2tree(taxize_class, check = TRUE)
taxize_phy <- taxize_tree[["phylo"]]
taxa_names(GP) <- taxa_names(taxize_phy)
GP <- merge_phyloseq(GP, taxize_phy)
#quality control 
summary(sample_sums(GP))

#rarefaction 
otu_tab <- t(abundances(GP))
p <- vegan::rarecurve(otu_tab, 
                      step = 50, label = FALSE, 
                      sample = min(rowSums(otu_tab), 
                                   col = "blue", cex = 0.6))

##filter sample to the lowest depth 
GP_rar <- rarefy_even_depth(GP, sample.size = min(sample_sums(GP)))

#rarefaction for filtered basesd on sample size 
otu_tab_rar <- t(abundances(GP_rar))
p_rar <- vegan::rarecurve(otu_tab_rar, 
                      step = 50, label = FALSE, 
                      sample = min(rowSums(otu_tab_rar), 
                                   col = "blue", cex = 0.6))


# quick check taxa prevalence

p_rar_t <- plot_taxa_prevalence(GP.rar, "Family")

pdf("family_prevalance.pdf", width = 40, height = 20)
print(p_rar_t)
dev.off()
###Alpha analysis 
GP_rar_a_div <- alpha(GP_rar, index = "all")
GP_rar_a_div <- rownames_to_column(GP_rar_a_div,var = "Sample_name")
# get the metadata out as seprate object
GP_rar_meta <- meta(GP_rar)
GP_rar_meta <- rownames_to_column(GP_rar_meta,var = "Sample_name")
# merge these two data frames into one
GP_rar_a_div_meta <- left_join(GP_rar_a_div,GP_rar_meta, by = "Sample_name")



# Now use this data frame to plot 

alpha_dis <- c("observed", "diversity_shannon","evenness_pielou" )
pdf("alpha_diversity.pdf")
for (alpha in alpha_dis){
p <- ggboxplot(GP_rar_a_div_meta, 
               x = "Experiment", 
               y = alpha,
               fill = "Treatment", 
               palette = "jco")

#p <- p + rotate_x_text()

print(p)
}
dev.off()

##Beta diversity 

Bray_cruits <- ordinate(GP_rar, "PCoA", "bray")
Jaccard <- ordinate(GP_rar, "PCoA", "jaccard")
Unweighted_UniFrac <- ordinate(GP_rar, "PCoA", "unifrac", weighted=F)
Weighted_UniFrac <- ordinate(GP_rar, "PCoA", "unifrac", weighted=T)

betas <- list(Bray_cruits, Jaccard, Unweighted_UniFrac ,  Weighted_UniFrac)
names(betas) <- c("Bray_cruits", "Jaccard", "Unweighted_UniFrac" ,  "Weighted_UniFrac")

pdf("PCoA_beta.pdf")
for (i in 1:length(betas)){

p <- plot_ordination(GP_rar,betas[[i]] , color="Experiment") 
p <- p + ggtitle(names(betas[i])) + geom_point(size = 2)
p <- p + theme_classic() + scale_color_brewer("Experiment", palette = "Set2")
print(p)

}
for (i in 1:length(betas)){
  
  p <- plot_ordination(GP_rar,betas[[i]] , color="Experiment") 
  p <- p + ggtitle(names(betas[i])) + geom_point(size = 2)
  p <- p + theme_classic() + scale_color_brewer("Experiment", palette = "Set2") + stat_ellipse()
  print(p)
  
}

dev.off()
# 
# alpha_meas = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")
# 
# pdf("Richness.pdf")
# par(mar = c(1, 1, 1, 1))
# plot_richness(GP,"Experiment","Treatment", measures=alpha_meas)
# plot_richness(GP,"Treatment","Experiment", measures=alpha_meas)
# dev.off()


#taxonomy bar plot 
plot_bar(GP_pd1, fill = "taxonomy3")


#for sub sampling 

for (e in c("aPD1", "Ss", "ST") ){
  dir.create(e)
  setwd(e)
  GP_sub <- subset_samples(GP, Experiment == e)
  
  pdf ("rarefaction.pdf")
  #rarefaction 
  otu_tab <- t(abundances(GP_sub))
  p <- vegan::rarecurve(otu_tab, 
                        step = 50, label = FALSE, 
                        sample = min(rowSums(otu_tab), 
                                     col = "blue", cex = 0.6))
  
  ##filter sample to the lowest depth 
  GP_rar <- rarefy_even_depth(GP_sub, sample.size = min(sample_sums(GP_sub)))
  
  #rarefaction for filtered basesd on sample size 
  otu_tab_rar <- t(abundances(GP_rar))
  p_rar <- vegan::rarecurve(otu_tab_rar, 
                            step = 50, label = FALSE, 
                            sample = min(rowSums(otu_tab_rar), 
                                         col = "blue", cex = 0.6))
  dev.off()
  
  # quick check taxa prevalence
  
  #p_rar_t <- plot_taxa_prevalence(GP.rar, "Family")
  
  pdf("family_prevalance.pdf", width = 40, height = 20)
  print(p_rar_t)
  dev.off()
  ###Alpha analysis 
  GP_rar_a_div <- alpha(GP_rar, index = "all")
  GP_rar_a_div <- rownames_to_column(GP_rar_a_div,var = "Sample_name")
  # get the metadata out as seprate object
  GP_rar_meta <- meta(GP_rar)
  GP_rar_meta <- rownames_to_column(GP_rar_meta,var = "Sample_name")
  # merge these two data frames into one
  GP_rar_a_div_meta <- left_join(GP_rar_a_div,GP_rar_meta, by = "Sample_name")
  
  
  
  # Now use this data frame to plot 
  
  alpha_dis <- c("observed", "diversity_shannon","evenness_pielou" )
  pdf("alpha_diversity.pdf")
  for (alpha in alpha_dis){
    p <- ggboxplot(GP_rar_a_div_meta, 
                   x = "Treatment", 
                   y = alpha,
                   fill = "Treatment", 
                   palette = "jco")
    
    #p <- p + rotate_x_text()
    
    print(p)
  }
  dev.off()
  
  ##Beta diversity 
  
  Bray_cruits <- ordinate(GP_rar, "PCoA", "bray")
  Jaccard <- ordinate(GP_rar, "PCoA", "jaccard")
  Unweighted_UniFrac <- ordinate(GP_rar, "PCoA", "unifrac", weighted=F)
  Weighted_UniFrac <- ordinate(GP_rar, "PCoA", "unifrac", weighted=T)
  
  betas <- list(Bray_cruits, Jaccard, Unweighted_UniFrac ,  Weighted_UniFrac)
  names(betas) <- c("Bray_cruits", "Jaccard", "Unweighted_UniFrac" ,  "Weighted_UniFrac")
  
  pdf("PCoA_beta.pdf")
  for (i in 1:length(betas)){
    
    p <- plot_ordination(GP_rar,betas[[i]] , color="Treatment") 
    p <- p + ggtitle(names(betas[i])) + geom_point(size = 2)
    p <- p + theme_classic() + scale_color_brewer("Treatment", palette = "Set2")
    print(p)
    
  }
  for (i in 1:length(betas)){
    
    p <- plot_ordination(GP_rar,betas[[i]] , color="Treatment") 
    p <- p + ggtitle(names(betas[i])) + geom_point(size = 2)
    p <- p + theme_classic() + scale_color_brewer("Treatment", palette = "Set2") + stat_ellipse()
    print(p)
    
  }
  
  dev.off()
  setwd("..")
}










pdf("Richness_GP_pd1.pdf",width = 20)
plot_richness(GP_pd1,"Treatment","Treatment", measures=alpha_meas)+ geom_boxplot() + theme_bw()
dev.off()


#build network 
ig <- make_network(GP, max.dist=0.3)
plot_network(ig, GP)





# 
# 
# 
# GP_rar_asvtab <- as.data.frame(GP_rar@otu_table)
# 
# GP_rar_tree <- as.phylo(GP_rar@phy_tree)
# 
# # We first need to check if the tree is rooted or not 
# 
# GP_rar_tree
# 
# # it is a rooted tree
# GP_pd <- pd(t(GP_rar_asvtab), GP_rar_tree, include.root = T) # t(ou_table) transposes the table for use in picante and the tre file comes from the first code chunck we used to read tree file (see making a phyloseq object section).
# datatable(df.pd)
# 
# 
# library(taxize)
# 
# taxize_tree2 <- class2tree(taxize_class, check = F)
# taxize_tree2 <- as.phylo(taxize_tree)
# 
# 
# taxize_tree2 <- as.phylo(taxize_tree)
# 
# GP <- merge_phyloseq(GP, x)
# 
# plot(taxize_oil_tree)
# #dna <- Biostrings::DNAStringSet(taxa_names(GP))
# #head(phy_tree(GP))
# refseqs <- bold_search(id = taxa_names(GP))
# dna <- Biostrings::DNAStringSet(taxa_names(GP))
# 
# 
# 
# library(metacoder)
# 
# taxon_metacoder <- extract_taxonomy(taxa_names(GP), key = "id", database = "ncbi",  allow_na = TRUE)
# ref <- metacoder::ncbi_taxon_sample(id = taxa_names(GP),target_rank = "species" ,max_counts = 1)
# ref <-lapply(taxa_names(GP) , function (x) {ncbi_taxon_sample(id = x,target_rank = "species" ,max_counts = 1)})
# x <- metacoder::parse_phyloseq(GP)
# y <- metacoder::as_phyloseq(x)
# 
# 
# library(tidyverse)
# library(ggtree)
# BiocManager::install("ggtree")
# 
# library(phylotools)
# x <- treeio::read.phylip("phyliptree.phy")
