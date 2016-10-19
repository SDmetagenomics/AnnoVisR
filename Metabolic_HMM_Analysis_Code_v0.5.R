library(plyr)
library(ggplot2)
library(pheatmap)
library(reshape)
library(vegan)

setwd("~/Desktop/testbed/AnnotateR_Out/Metabolic_HMM_output/")




#######
#### Load Metadata
#######


Genome_metadata <- read.delim(file.choose(), header = FALSE)
colnames(Genome_metadata) <- c("genome_name", "protein_count", "completeness", "contamination", "Factor")
heat_anno <- data.frame(Class = Genome_metadata$Factor)
row.names(heat_anno) <- Genome_metadata$genome_name




#######
####  Subset Metadata For Good Genomes (>= 70% Complete ; <= 10% Contamination)
#######


Genome_metadata <- subset(Genome_metadata, completeness >= 70 & contamination <= 10)
tmp_metadata_factor1 <- subset(Genome_metadata, Factor == "Increase")
tmp_metadata_factor2 <- subset(Genome_metadata, Factor == "Decrease")


ks_test_p <- ks.test(tmp_metadata_factor1$completeness, tmp_metadata_factor2$completeness, exact = FALSE)


ggplot(Genome_metadata, aes(x = completeness, fill = Factor)) + 
  geom_density( alpha = 0.5) +
  ggtitle("Genome Completeness Density") + 
  scale_fill_manual(values = c("steelblue", "firebrick3")) +
  theme(axis.text.x = element_text(colour = "black", size = 12, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 12)) +
  annotate("text", x = 80, y = 0.055, size = 5, label = paste0("KS Test p-Value: ",round(ks_test_p$p.value, 4)))

#ggsave(paste0("~/Desktop/Temp_R_Plots/Sample_Distribution_Stats/Genome_Completeness_dens.pdf"))


ks_test_p <- ks.test(tmp_metadata_factor1$contamination, tmp_metadata_factor2$contamination, exact = FALSE)


ggplot(Genome_metadata, aes(x = contamination, fill = Factor)) + 
  geom_density( alpha = 0.5) +
  ggtitle("Genome Contamination Density") + 
  scale_fill_manual(values = c("steelblue", "firebrick3")) +
  theme(axis.text.x = element_text(colour = "black", size = 12, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 12)) +
  annotate("text", x = 2.5, y = 0.20, size = 5, label = paste0("KS Test p-Value: ",round(ks_test_p$p.value, 4)))

#ggsave(paste0("~/Desktop/Temp_R_Plots/Sample_Distribution_Stats/Genome_Contam_dens.pdf"))

rm("ks_test_p", "tmp_metadata_factor1", "tmp_metadata_factor2")




#######
####  Count Hits to Karthik Metabolic HMMs from Filtered Genomes as well as Main and Sub Role Hits
#######


### Count HMM Hits per Genome and Build Master Table

Master_Metabolic_Hits <- read.delim("~/Desktop/testbed/AnnoVisR/Annotation_Templates/Metabolic_HMM_Master.txt", header = TRUE)
genome_names <- as.vector(Genome_metadata$genome_name)

for (i in 1:length(genome_names)){
  temp_table <- read.table(genome_names[i])
  temp_counts <- count(temp_table$V2)
  colnames(temp_counts) <- c("Codes", genome_names[i])
  Master_Metabolic_Hits <- merge(Master_Metabolic_Hits, temp_counts, by="Codes", all= TRUE)
}

rm("temp_counts","temp_table", "i")


### Make NA == 0s and Parse out Rows of All 0s and Report Genes that had No Hit

Master_Metabolic_Hits[is.na(Master_Metabolic_Hits)] <- 0
No_Hit <- rowSums(Master_Metabolic_Hits[,5:ncol(Master_Metabolic_Hits)]) == 0
No_Hit <- Master_Metabolic_Hits[No_Hit,1:4]
Master_Metabolic_Hits <- Master_Metabolic_Hits[!!rowSums(abs(Master_Metabolic_Hits[-c(1:4)])),]


### Calculate Numbers of Hits within Main and Sub Roles

Main_Categories <- read.delim("~/Desktop/testbed/AnnoVisR/Annotation_Templates/Metabolic_HMM_Main_Category_Template.txt", header = TRUE)
Sub_Categories <- read.delim("~/Desktop/testbed/AnnoVisR/Annotation_Templates/Metabolic_HMM_Sub_Category_Template.txt", header = TRUE)


Main_Category_Hits <- data.frame()

for (i in 1:nrow(Main_Categories)){
  Main_Subset <- Master_Metabolic_Hits[grepl(Main_Categories[i,1], Master_Metabolic_Hits$Main_Category),]
  temp_row <- colSums(Main_Subset[5:ncol(Main_Subset)])
  Main_Category_Hits <- rbind(Main_Category_Hits, temp_row)
}

rm("Main_Subset", "temp_row")


Sub_Category_Hits <- data.frame()

for (i in 1:nrow(Sub_Categories)){
  Sub_Subset <- Master_Metabolic_Hits[grepl(Sub_Categories[i,1], Master_Metabolic_Hits$Sub_Category),]
  temp_row <- colSums(Sub_Subset[5:ncol(Sub_Subset)])
  Sub_Category_Hits <- rbind(Sub_Category_Hits, temp_row)
}

rm("Sub_Subset", "temp_row")


rownames(Sub_Category_Hits) <- Sub_Categories$Sub_Category
colnames(Sub_Category_Hits) <- genome_names


rownames(Main_Category_Hits) <- Main_Categories$Main_Category
colnames(Main_Category_Hits) <- genome_names

rm("Main_Categories", "Sub_Categories")




#######
#### Calculate Summary Statistics for Genomes with Respect to Metabolic Enzyme Number and Diversity
#######


Master_Metabolic_Matrix <- Master_Metabolic_Hits[,5:ncol(Master_Metabolic_Hits)]

Metabolic_Variety <- data.frame(Genome_name = Genome_metadata$genome_name, Factor = Genome_metadata$Factor,
                           Protein_count = Genome_metadata$protein_count,
                           Unique_Enzymes = NA, All_Enzymes = NA)


### Populate Metabolic_Variety DF with Absolute and Unique Domain Hit Numbers 

for (i in 1:ncol(Master_Metabolic_Matrix)){
  tmp_count <- sum(ifelse(Master_Metabolic_Matrix[,i] > 0, 1, 0))
  Metabolic_Variety[i,4] <- tmp_count
  tmp_count <- sum(Master_Metabolic_Matrix[,i])
  Metabolic_Variety[i,5] <- tmp_count
}


Metabolic_Variety <- data.frame(Metabolic_Variety, 
                           Unique_Enzymes_protnorm = (Metabolic_Variety$All_Enzymes/Metabolic_Variety$Protein_count)*1000, 
                           All_Enzymes_protnorm = (Metabolic_Variety$Unique_Enzymes/Metabolic_Variety$Protein_count)*1000)


### Produce Diagnostic Scatter Plots Comparing Various Count Statistics

ggplot(Metabolic_Variety, aes(x = Protein_count, y = All_Enzymes, color = Factor)) + 
  geom_point(size = 3, alpha = 0.8) +
  stat_smooth(method = "lm") +
  scale_color_manual(values = c("steelblue", "firebrick3")) +
  ggtitle("Total CAZy Enzymes vs Genome Protein Number")


ggplot(Metabolic_Variety, aes(x = Protein_count, y = Unique_Enzymes, color = Factor)) + 
  geom_point(size = 3, alpha = 0.8) +
  stat_smooth(method = "lm") +
  scale_color_manual(values = c("steelblue", "firebrick3")) +
  ggtitle("Unique CAZy Enzymes vs Genome Protein Number")


ggplot(Metabolic_Variety, aes(x = Unique_Enzymes , y = All_Enzymes, color = Factor)) + 
  geom_point(size = 3, alpha = 0.8) +
  stat_smooth(method = "lm") +
  scale_color_manual(values = c("steelblue", "firebrick3")) +
  ggtitle("Unique CAZy vs Total CAZy")


ggplot(Metabolic_Variety, aes(x = All_Enzymes_protnorm, y = Unique_Enzymes_protnorm, color = Factor)) + 
  geom_point(size = 3, alpha = 0.8) +
  stat_smooth(method = "lm") +
  scale_color_manual(values = c("steelblue", "firebrick3")) +
  ggtitle("Total vs Unique CAZy Enzymes Prot_Normalized")


ggplot(melt(Metabolic_Variety[,2:7]), aes(x = Factor, y = value, fill = Factor)) +
  geom_boxplot() +
  facet_wrap(~ variable, scales = "free_y") +
  ggtitle("Group CAZy Count Statistics") +
  scale_fill_manual(values = c("steelblue", "firebrick3"))






############





Enzyme_Diversity <- data.frame(Factor = Genome_metadata$Factor, 
                               Shannon = diversity(t(Master_Metabolic_Matrix)),
                               Simpson = diversity(t(Master_Metabolic_Matrix), index = "simpson"),
                               InvSimpson = diversity(t(Master_Metabolic_Matrix), index = "invsimpson"))


ggplot(melt(Enzyme_Diversity), aes(x = Factor, y = value, fill = Factor)) +
  geom_boxplot() + 
  facet_wrap(~variable, scales = "free_y") +
  ggtitle("Alpha Diversity Indicies for CAZy Enzymes in Sample Groups") +
  scale_fill_manual(values = c("steelblue", "firebrick3"))


#######
####  Count Hits to Karthik Metabolic HMMs from Filtered Genomes as well as Main and Sub Role Hits
#######



pheatmap(log2(Main_Category_Hits+1), 
         annotation_col = heat_anno,
         clustering_method = "ward.D"
         )
pheatmap(Sub_Category_Hits)

Log_Main_Category <- Main_Category_Hits + 1
Log_Main_Category <- apply(Log_Main_Category, 2, log2)
pheatmap(Log_Main_Category, clustering_method = "ward.D", annotation_col = `20cm_Effect`)


Log_Sub_Category <- Sub_Category_Hits + 1
Log_Sub_Category <- apply(Log_Sub_Category, 2 , log2)
pheatmap(Log_Sub_Category, clustering_method = "ward.D2", annotation_col = `20cm_Effect`)

pheatmap(t(ordered_matrix), cluster_cols = FALSE, cluster_rows = FALSE, annotation_col = `20cm_Effect`)
