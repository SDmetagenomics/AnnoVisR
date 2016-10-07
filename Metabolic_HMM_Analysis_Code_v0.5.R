library(plyr)
library(ggplot2)
library(pheatmap)

setwd("~/Desktop/Metabolic_HMM_summary/")

## Template Tables ##
Metabolic_Hits <- read.delim("~/Desktop/testbed/Metabolic_HMM_Template1.txt", header = TRUE)
genome_names <- list.files()
for (i in 1:length(genome_names)){
  temp_table <- read.table(genome_names[i])
  temp_counts <- count(temp_table$V2)
  colnames(temp_counts) <- c("Codes", genome_names[i])
  Metabolic_Hits <- merge(Metabolic_Hits, temp_counts, by="Codes", all= TRUE)
}

rm("temp_counts","temp_table")

# Make NA == 0s and parse out rows of all 0s
Metabolic_Hits[is.na(Metabolic_Hits)] <- 0

# Load template sheets where numbers of genes within main and sub roles will be calculated
Main_Categories <- read.delim("~/Desktop/testbed/Metabolic_HMM_Main_Category_Template.txt", header = TRUE)
Sub_Categories <- read.delim("~/Desktop/testbed/Metabolic_HMM_Sub_Category_Template.txt", header = TRUE)

Main_Category_Hits <- data.frame()

for (i in 1:nrow(Main_Categories)){
  Main_Subset <- Metabolic_Hits[grepl(Main_Categories[i,1], Metabolic_Hits$Main_Category),]
  temp_row <- colSums(Main_Subset[5:ncol(Main_Subset)])
  Main_Category_Hits <- rbind(Main_Category_Hits, temp_row)
}

rm("Main_Subset", "temp_row")

###

Sub_Category_Hits <- data.frame()

for (i in 1:nrow(Sub_Categories)){
  Sub_Subset <- Metabolic_Hits[grepl(Sub_Categories[i,1], Metabolic_Hits$Sub_Category),]
  temp_row <- colSums(Sub_Subset[5:ncol(Sub_Subset)])
  Sub_Category_Hits <- rbind(Sub_Category_Hits, temp_row)
}

rm("Sub_Subset", "temp_row")

rownames(Sub_Category_Hits) <- Sub_Categories$Sub_Category
colnames(Sub_Category_Hits) <- genome_names


rownames(Main_Category_Hits) <- Main_Categories$Main_Category
colnames(Main_Category_Hits) <- genome_names

pheatmap(Main_Category_Hits, annotation_col = `20cm_Effect`)
pheatmap(Sub_Category_Hits)

Log_Main_Category <- Main_Category_Hits + 1
Log_Main_Category <- apply(Log_Main_Category, 2, log2)
pheatmap(Log_Main_Category, clustering_method = "ward.D", annotation_col = `20cm_Effect`)


Log_Sub_Category <- Sub_Category_Hits + 1
Log_Sub_Category <- apply(Log_Sub_Category, 2 , log2)
pheatmap(Log_Sub_Category, clustering_method = "ward.D2", annotation_col = `20cm_Effect`)

pheatmap(t(ordered_matrix), cluster_cols = FALSE, cluster_rows = FALSE, annotation_col = `20cm_Effect`)
