library(plyr)
library(ggplot2)
library(pheatmap)
library(reshape)
library(vegan)

# gets wd for local execution and sets it to the working directory

#PWD <- getwd()
setwd("~/Desktop/testbed/AnnotateR_Out/CAZy_HMM_output/") # This should be changed to PWD for later




#######
#### Load Metadata
#######


Genome_metadata <- read.delim(file.choose(), header = FALSE)
colnames(Genome_metadata) <- c("genome_name", "protein_count", "completeness", "contamination", "exp_plot")
heat_anno <- data.frame(Class = Genome_metadata$exp_plot)
row.names(heat_anno) <- Genome_metadata$genome_name




#######
####  Subset Metadata For Good Genomes (>= 70% Complete ; <= 10% Contamination)
#######


Genome_metadata <- subset(Genome_metadata, completeness >= 70 & contamination <= 10)

ggplot(Genome_metadata, aes(x = completeness, fill = exp_plot)) + 
  geom_density( alpha = 0.5) +
  ggtitle("Genome Completeness Density") + 
  scale_fill_manual(values = c("steelblue", "firebrick3")) +
  theme(axis.text.x = element_text(colour = "black", size = 12, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 12))

ggsave(paste0("~/Desktop/Temp_R_Plots/Sample_Distribution_Stats/Genome_Completeness_dens.pdf"))


ggplot(Genome_metadata, aes(x = contamination, fill = exp_plot)) + 
  geom_density( alpha = 0.5) +
  ggtitle("Genome Contamination Density") + 
  scale_fill_manual(values = c("steelblue", "firebrick3")) +
  theme(axis.text.x = element_text(colour = "black", size = 12, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 12))

ggsave(paste0("~/Desktop/Temp_R_Plots/Sample_Distribution_Stats/Genome_Contam_dens.pdf"))




#######
####  Build Single Domain Hit Table
#######


CAZy_SingleDom <- read.delim("~/Desktop/testbed/AnnoVisR/Annotation_Templates/CAZy_Annotation_Template.txt", header = TRUE)

genome_names <- as.vector(Genome_metadata$genome_name)

for (i in 1:length(genome_names)){
  temp_table1 <- read.table(genome_names[i])
  temp_table1 <- temp_table1[!duplicated(temp_table1[c("V1", "V3")]),]
  temp_counts1 <- count(temp_table1$V1)
  colnames(temp_counts1) <- c("CAZy_Family", genome_names[i])
  CAZy_SingleDom <- merge(CAZy_SingleDom, temp_counts1, by="CAZy_Family", all=TRUE)
}

rm("temp_counts1","temp_table1")


## Make NA == 0 & Parse out CAZy enzymes with no hits in any organism (All Zero Rows)

CAZy_SingleDom[is.na(CAZy_SingleDom)] <- 0
CAZy_SingleDom <- CAZy_SingleDom[!!rowSums(abs(CAZy_SingleDom[-c(1:4)])),]

write.table(CAZy_SingleDom, "~/Desktop/Temp_R_Plots/Output_Tables/Domain_Hits_per_Genome.txt", sep = "\t", row.names = FALSE)




#######
#### Calculate Summary Statistics for Genomes and CAZy Enzymes
#######


Single_Dom_Matrix <- CAZy_SingleDom[,5:ncol(CAZy_SingleDom)]

CAZy_Variety <- data.frame(Genome_name = Genome_metadata$genome_name, Factor = Genome_metadata$exp_plot,
                           Protein_count = Genome_metadata$protein_count,
                           Unique_cazy = NA, All_cazy = NA)


### Populate CAZy_Variety DF with Absolute and Unique Domain Hit Numbers 

for (i in 1:ncol(Single_Dom_Matrix)){
  tmp_count <- sum(ifelse(Single_Dom_Matrix[,i] > 0, 1, 0))
  CAZy_Variety[i,4] <- tmp_count
  tmp_count <- sum(Single_Dom_Matrix[,i])
  CAZy_Variety[i,5] <- tmp_count
}

### Add cazy normalized by genome/protein counts

CAZy_Variety <- data.frame(CAZy_Variety, 
                           Unique_cazy_protnorm = (CAZy_Variety$All_cazy/CAZy_Variety$Protein_count)*1000, 
                           All_cazy_protnorm = (CAZy_Variety$Unique_cazy/CAZy_Variety$Protein_count)*1000)


### Produce Box-Plots of CAZy Count Statistics by Group Factor

ggplot(melt(CAZy_Variety[,2:7]), aes(x = Factor, y = value, fill = Factor)) +
         geom_boxplot() +
         facet_wrap(~ variable, scales = "free_y") +
         ggtitle("Group CAZy Count Statistics") +
         scale_fill_manual(values = c("steelblue", "firebrick3"))

ggsave(paste0("~/Desktop/Temp_R_Plots/CAZy_Summary_Stats/Group_CAZy_Count_Stats.pdf"))


### Calculate F and t Statistics Comparing CAZy Counts by Group Factor

CAZy_Count_Fandt_Stats <- data.frame(Test_Param = (1:5),F_test = (1:5), t_test = (1:5), Fold_Diff = (1:5))

tmp_decrease <- subset(CAZy_Variety, Factor == "Decrease")[,3:7]
tmp_increase <- subset(CAZy_Variety, Factor == "Increase")[,3:7]

for (i in 1:5){
 test_param_tmp <- colnames(tmp_decrease[i])
 F_test_tmp <- var.test(tmp_increase[,i], tmp_decrease[,i])
 var_equal <- ifelse(F_test_tmp$p.value < 0.05, FALSE, TRUE)
 T_test_tmp <- t.test(tmp_increase[,i], tmp_decrease[,i], var.equal = var_equal)
 CAZy_Count_Fandt_Stats[i,1] <- test_param_tmp
 CAZy_Count_Fandt_Stats[i,2] <- F_test_tmp$p.value
 CAZy_Count_Fandt_Stats[i,3] <- T_test_tmp$p.value
 CAZy_Count_Fandt_Stats[i,4] <- unname(T_test_tmp$estimate[1]/T_test_tmp$estimate[2])
}


rm("tmp_decrease", "tmp_increase", "F_test_tmp", "T_test_tmp", "var_equal", "test_param_tmp")

write.table(CAZy_Count_Fandt_Stats, "~/Desktop/Temp_R_Plots/CAZy_Summary_Stats/Group_Count_T_Tests.txt", 
            sep = "\t", row.names = FALSE)


#### Produce Diagnostic Scatter Plots Comparing Various Count Statistics

ggplot(CAZy_Variety, aes(x = Protein_count, y = Unique_cazy, color = Factor)) + 
  geom_point(size = 3, alpha = 0.8) +
  stat_smooth(method = "lm") +
  scale_color_manual(values = c("steelblue", "firebrick3")) +
  ggtitle("Unique CAZy Enzymes vs Genome Protein Number")

ggsave(paste0("~/Desktop/Temp_R_Plots/CAZy_Summary_Stats/Unique_CAZy_vs_Protein.pdf"))

ggplot(CAZy_Variety, aes(x = Protein_count, y = All_cazy, color = Factor)) + 
  geom_point(size = 3, alpha = 0.8) +
  stat_smooth(method = "lm") +
  scale_color_manual(values = c("steelblue", "firebrick3")) +
  ggtitle("Total CAZy Enzymes vs Genome Protein Number")

ggsave(paste0("~/Desktop/Temp_R_Plots/CAZy_Summary_Stats/Total_CAZy_vs_Protein.pdf"))

ggplot(CAZy_Variety, aes(x = All_cazy, y = Unique_cazy, color = Factor)) + 
  geom_point(size = 3, alpha = 0.8) +
  stat_smooth(method = "lm") +
  scale_color_manual(values = c("steelblue", "firebrick3")) +
  ggtitle("Unique CAZy vs Total CAZy")

ggsave(paste0("~/Desktop/Temp_R_Plots/CAZy_Summary_Stats/Unique_CAZy_vs_Total_CAZy.pdf"))

ggplot(CAZy_Variety, aes(x = All_cazy/Protein_count, y = Unique_cazy/Protein_count, color = Factor)) + 
  geom_point(size = 3, alpha = 0.8) +
  stat_smooth(method = "lm") +
  scale_color_manual(values = c("steelblue", "firebrick3")) +
  ggtitle("Total vs Unique CAZy Enzymes Prot_Normalized")

ggsave(paste0("~/Desktop/Temp_R_Plots/CAZy_Summary_Stats/Unique_CAZy_vs_Total_CAZy_ProtNorm.pdf"))


#### Compare alpha-diversity Between Factors Based on CAZy Enzyme Abundance

rownames(Single_Dom_Matrix) <- CAZy_SingleDom$CAZy_Family

Enzyme_Diversity <- data.frame(Factor = Genome_metadata$exp_plot, 
                               Shannon = diversity(t(Single_Dom_Matrix)),
                               Simpson = diversity(t(Single_Dom_Matrix), index = "simpson"),
                               InvSimpson = diversity(t(Single_Dom_Matrix), index = "invsimpson"))

ggplot(melt(Enzyme_Diversity), aes(x = Factor, y = value, fill = Factor)) +
  geom_boxplot() + 
  facet_wrap(~variable, scales = "free_y") +
  ggtitle("Alpha Diversity Indicies for CAZy Enzymes in Sample Groups") +
  scale_fill_manual(values = c("steelblue", "firebrick3"))

ggsave(paste0("~/Desktop/Temp_R_Plots/CAZy_Summary_Stats/CAZy_Enzyme_Diversity_by_Group.pdf"))



#test_matrix <- t(Single_Dom_Matrix)
#test_matrix <- data.frame(Factor = Genome_metadata$exp_plot, test_matrix)

#test_matrix2 <- subset(test_matrix, Factor == "Decrease")

#All_CAZy = 20677
#All_GH109 = 1264
#GH_109_In_Set = 881
#GH_109_Not_in_Set = 383   # Col 35
#All_CAZy_In_UP = 12642
#All_CAZy_In_Down = 
#All_GH109_In_Down = 383

# example matrix
#
#                         PART OF GH109      NOT PART OF GH109
#         
# IN ENRICHED SET         881                 11761
# NOT ENRICHED SET        383                 8035
#

#fisher.test(matrix(c(881,383,11761,8035), nrow=2, ncol=2))


####
#### Calculate the number of general substrate class hits per sample
####

Substrate_Classes <- read.delim("~/Desktop/testbed/AnnoVisR/Annotation_Templates/CAZy_Substrate_Class.txt", header = TRUE) # Change to biotite dir
Substrate_Class_Hits <- data.frame()

for (i in 1:nrow(Substrate_Classes)){
  Class_Subset <- CAZy_SingleDom[grepl(Substrate_Classes[i,1], CAZy_SingleDom$Substrate_Class),]
  temp_row <- colSums(Class_Subset[5:ncol(Class_Subset)])
  Substrate_Class_Hits <- rbind(Substrate_Class_Hits, temp_row)
}

Substrate_Class_Hits <- data.frame(Substrate_Class = Substrate_Classes$Substrate_Class, Substrate_Class_Hits)
colnames(Substrate_Class_Hits)[2:ncol(Substrate_Class_Hits)] <- genome_names

rm("Substrate_Classes", "Class_Subset")

write.table(Substrate_Class_Hits, "~/Desktop/Temp_R_Plots/Output_Tables/Substrate_Class_per_Genome.txt", sep = "\t", row.names = FALSE)






####
#### Plot total number of proteins in genome set acting on a general substrate class
####

Total_Class_Counts <- data.frame(Counts = rowSums(Substrate_Class_Hits[,2:ncol(Substrate_Class_Hits)]))
Total_Class_Counts <- data.frame(Substrate_Class = Substrate_Class_Hits$Substrate_Class, Counts = Total_Class_Counts$Counts)

ggplot(Total_Class_Counts, aes(x = reorder(Substrate_Class, -Counts), y = Counts)) + 
                     geom_bar(stat = "identity", fill = "steelblue4", color = "black") +
                     ggtitle("Enzymes by Substrate Class") +
                     xlab(NULL) + ylab("Protein Counts") +
                     theme(axis.text.x = element_text(colour = "black", size = 12, angle = 45, hjust = 1),
                           axis.text.y = element_text(colour = "black", size = 12))

ggsave(paste0("~/Desktop/Temp_R_Plots/Emzymes_by_Substrate_Class.pdf"))


ggplot(Total_Class_Counts, aes(x = reorder(Substrate_Class, -Counts), y = Counts/sum(Total_Class_Counts$Counts))) + 
  geom_bar(stat = "identity", fill = "steelblue4", color = "black") +
  ggtitle("Fraction of Enzymes per Substrate Class") +
  xlab(NULL) + ylab("Fraction of Proteins") +
  theme(axis.text.x = element_text(colour = "black", size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 12))

ggsave(paste0("~/Desktop/Temp_R_Plots/Enzyme_Fraction_for_Substrate_Class.pdf"))

Total_Class_Counts <- data.frame(Total_Class_Counts, Fraction = Total_Class_Counts$Counts/sum(Total_Class_Counts$Counts))

write.table(Total_Class_Counts, "~/Desktop/Temp_R_Plots/Output_Tables/Substrate_Class_Counts.txt", sep = "\t", row.names = FALSE)






####
#### Build Heatmap of Substrate_Class_Hits per genome
####  

rownames(Substrate_Class_Hits) <- Substrate_Class_Hits$Substrate_Class
Substrate_Class_Hits <- as.matrix(Substrate_Class_Hits[,2:ncol(Substrate_Class_Hits)])

ann_cols <- list(Class = c(Decrease = "steelblue", Increase = "firebrick3"))

drows <- vegdist(Substrate_Class_Hits, method = "binomial")
dcols <- vegdist(t(Substrate_Class_Hits), method = "binom")

pheatmap(Substrate_Class_Hits, scale = "column", clustering_distance_rows = drows,
         clustering_distance_cols = dcols,
         clustering_method = "ward.D",
         margins = (c(10,10)),
         annotation_col = heat_anno, 
         annotation_colors = ann_cols)






####  
#### Calculation of Specific Substrates Hit per Genome
####

Substrates <- read.delim("~/Desktop/testbed/AnnoVisR/Annotation_Templates/CAZy_Substrate.txt", header = TRUE)
Substrate_Hits <- data.frame()

for (i in 1:nrow(Substrates)){
  Class_Subset <- CAZy_SingleDom[grepl(Substrates[i,1], CAZy_SingleDom$Substrate),]
  temp_row <- colSums(Class_Subset[5:ncol(Class_Subset)])
  Substrate_Hits <- rbind(Substrate_Hits, temp_row)
}

rm("Class_Subset", "temp_row", "i")


Total_Substrate_Counts <- data.frame(Counts = rowSums(Substrate_Hits))
Total_Substrate_Counts <- data.frame(Substrate = Substrates$Substrate, Counts = Total_Substrate_Counts$Counts)

ggplot(Total_Substrate_Counts, aes(x = reorder(Substrate, -Counts), y = Counts)) + 
  geom_bar(stat = "identity", fill = "steelblue4", color = "black") + xlab(NULL) +
  theme(axis.text.x = element_text(colour = "black", size = 12, angle = 90, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 12))

ggsave(paste0("~/Desktop/Temp_R_Plots/Emzymes_by_Substrate_Class.pdf"))

Total_Substrate_Counts <- data.frame(Total_Substrate_Counts, Fraction = Total_Substrate_Counts$Counts/sum(Total_Substrate_Counts$Counts))

write.table(Total_Substrate_Counts, "~/Desktop/Temp_R_Plots/Output_Tables/Specific_Substrate_Counts.txt", sep = "\t", row.names = FALSE)






####
#### Build table and heatmap of Specific Substrate Hits in each Genome
####

export_table <- data.frame(Substrate_Hits)
colnames(export_table) <- genome_names
export_table <- cbind(Substrates$Substrate, export_table)

write.table(Total_Substrate_Counts, "~/Desktop/Temp_R_Plots/Output_Tables/Specific_Substrate_by_Genome.txt", sep = "\t", row.names = FALSE)
rm(export_table)


Substrate_Hits <- as.matrix(Substrate_Hits)
rownames(Substrate_Hits) <- Substrates$Substrate
colnames(Substrate_Hits) <- genome_names


pheatmap(Substrate_Hits, margins = (c(10,10)),clustering_method = "ward.D",
         annotation_col = heat_anno, annotation_colors = ann_cols)

Substrate_Hits_Log <- Substrate_Hits + 1
Substrate_Hits_Log <- apply(Substrate_Hits_Log, 2, log2)

pheatmap(Substrate_Hits_Log, margins = (c(10,10)),clustering_method = "ward.D",
         annotation_col = heat_anno, annotation_colors = ann_cols)





####
#### Calculate and plot the total number of each enzyme type in this dataset
####

Enzyme_Classes <- c("GH", "GT", "CE", "PL", "AA")

for (i in 1:length(Enzyme_Classes)){
  keep_rows <- grep(Enzyme_Classes[i], CAZy_SingleDom$CAZy_Family)
  functional_subset <- CAZy_SingleDom[keep_rows,]
  function_sums <- data.frame(CAZy_Family = functional_subset$CAZy_Family, 
                              Counts = rowSums(functional_subset[,5:ncol(functional_subset)]))

  ggplot(function_sums, aes(x = reorder(CAZy_Family, -Counts), y = Counts)) +
         geom_bar(stat = "identity", fill = "steelblue4", color = "black") +
         xlab(NULL) +
         theme(axis.text.x = element_text(colour = "black", size = 10, angle = 90, hjust = 1),
         axis.text.y = element_text(colour = "black", size = 12))
  
  ggsave(paste0("~/Desktop/Temp_R_Plots/",Enzyme_Classes[i],"_Protein_Counts.pdf"), width = 14)
}


# ## Build Multiple Heatmaps for number of each CAZy Enzyme Class
# 
# for (i in 1:length(Enzyme_Classes)){
#   keep_rows <- grep(Enzyme_Classes[i], CAZy_SingleDom$CAZy_Family)
#   functional_subset <- CAZy_SingleDom[keep_rows,]
#   subset_matrix <- as.matrix(functional_subset[,5:ncol(functional_subset)])
#   rownames(subset_matrix) <- functional_subset$CAZy_Family
# 
#   pheatmap(subset_matrix,scale = "none",clustering_method = "ward.D2", margins = c(10,10), 
#           annotation_col = `20cm_Effect`)
# 
# }
# 
# # One big Heatmap for all CAZy
# subset_matrix <- as.matrix(CAZy_SingleDom[,5:ncol(CAZy_SingleDom)])
# rownames(subset_matrix) <- CAZy_SingleDom$CAZy_Family
# pheatmap(subset_matrix,clustering_method = "ward.D", margins = c(10,10), annotation_col = `20cm_Effect`)

### Development Notes ###
## Need a way to label and specify significantly altered organisms from some background set 
## Add plots that asses class based on specific substrate in addition to broad substrate group
## 






####
#### Perform Enrichment Statistics between Metadata Labled Groups
####

wilx_stats_output <- data.frame()

for (i in 1: nrow(CAZy_SingleDom)){
  test <- subset(CAZy_SingleDom, CAZy_Family == CAZy_SingleDom[i,1])
  test <- melt(test)
  test <- data.frame(test, exp_plot = Genome_metadata$exp_plot)
  GH_up <- subset(test, exp_plot == "Increase")
  GH_down <- subset(test, exp_plot == "Decrease")
  t_test_out <- t.test(GH_down$value, GH_up$value)
  wilcox_test_out <- wilcox.test(GH_down$value, GH_up$value)
  CAZy_Family <- CAZy_SingleDom[i,1]
  p_value_wilx <- wilcox_test_out$p.value
  genome_mean_up <- mean(GH_up$value)
  genome_mean_down <- mean(GH_down$value)
  stats_tmp <- data.frame(CAZy_Family = CAZy_Family, genome_mean_up = genome_mean_up, 
                          genome_mean_down = genome_mean_down, 
                          p_value = p_value_wilx)
  wilx_stats_output <- rbind(wilx_stats_output, stats_tmp)
}

wilx_stats_output <- data.frame(wilx_stats_output, FDR = p.adjust(wilx_stats_output$p_value))

write.table(wilx_stats_output, "~/Desktop/Temp_R_Plots/Output_Tables/Wilx_Stats_Output.txt", sep = "\t", row.names = FALSE)

wilx_significant <- subset(wilx_stats_output, FDR <= 0.1)

wilx_significant <- merge(wilx_significant, CAZy_SingleDom, by = "CAZy_Family")


####
#### Build Heatmap of Significantly Different CAZy_Families for all Genomes
####

sig_heat <- as.matrix(wilx_significant[,10:ncol(wilx_significant)])
rownames(sig_heat) <- wilx_significant$CAZy_Family



pheatmap(sig_heat, margins = (c(10,10)),clustering_method = "ward.D",
         annotation_col = heat_anno, annotation_colors = ann_cols)


pheatmap(sig_heat, margins = (c(10,10)),clustering_method = "ward.D", scale = "row",
         annotation_col = heat_anno, annotation_colors = ann_cols)

sig_heat_Log <- sig_heat + 1
sig_heat_Log <- apply(sig_heat_Log, 2, log2)

pheatmap(sig_heat_Log, margins = (c(10,10)), clustering_method = "ward.D",
         annotation_col = heat_anno, annotation_colors = ann_cols)


####
#### Calculate number of Enzymes by General Substrate Class and Treatment in Significant Set 
####

Substrate_Classes <- read.delim("~/Desktop/testbed/AnnoVisR/Annotation_Templates/CAZy_Substrate_Class.txt", header = TRUE) # Change to biotite dir
Significant_Gen_Class <- data.frame()

for (i in 1:nrow(Substrate_Classes)){
  Class_Subset <- wilx_significant[grepl(Substrate_Classes[1,1], wilx_significant$Substrate_Class),]
  temp_row <- colSums(Class_Subset[10:ncol(Class_Subset)])
  temp_row <- as.data.frame(t(t(temp_row)))
  Significant_Gen_Class <- cbind(Significant_Gen_Class, temp_row$V1)
}
colnames(Significant_Gen_Class) <- Substrate_Classes$Substrate_Class


Substrate_Class_Hits <- data.frame(Substrate_Class = Substrate_Classes$Substrate_Class, Substrate_Class_Hits)
colnames(Substrate_Class_Hits)[2:ncol(Substrate_Class_Hits)] <- genome_names

rm("Substrate_Classes", "Class_Subset")

write.table(Substrate_Class_Hits, "~/Desktop/Temp_R_Plots/Output_Tables/Substrate_Class_per_Genome.txt", sep = "\t", row.names = FALSE)






