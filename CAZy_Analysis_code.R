library(plyr)
library(ggplot2)
library(pheatmap)


# gets wd for local execution and sets it to the working directory
#PWD <- getwd()
setwd("~/Desktop/CAZy_summary/") # This should be changed to PWD for later

## Template Tables##
CAZy_ALL <- read.delim("~/Desktop/testbed/CAZy_hit_template.txt", header = TRUE) # Change to biotite dir
CAZy_SingleDom <- read.delim("~/Desktop/testbed/CAZy_hit_template.txt", header = TRUE) # Change to biotite dir
  
## Build Tables##
genome_names <- list.files() # will pull from working dir by default
for (i in 1:length(genome_names)){
  temp_table1 <- read.table(genome_names[i])
  temp_table2 <- temp_table1
  temp_table2 <- temp_table2[!duplicated(temp_table2[c("V1", "V3")]),]
  temp_counts1 <- count(temp_table1$V1)
  temp_counts2 <- count(temp_table2$V1)
  colnames(temp_counts1) <- c("CAZy_Family", genome_names[i])
  colnames(temp_counts2) <- c("CAZy_Family", genome_names[i])
  CAZy_ALL <- merge(CAZy_ALL, temp_counts1, by="CAZy_Family", all=TRUE)
  CAZy_SingleDom <- merge(CAZy_SingleDom, temp_counts2, by="CAZy_Family", all=TRUE)
}
rm("temp_counts1","temp_counts2", "temp_table1", "temp_table2")


## Make NA == 0 & Parse out All Zero Rows
CAZy_ALL[is.na(CAZy_ALL)] <- 0
CAZy_ALL <- CAZy_ALL[!!rowSums(abs(CAZy_ALL[-c(1:4)])),]
CAZy_SingleDom[is.na(CAZy_SingleDom)] <- 0
CAZy_SingleDom <- CAZy_SingleDom[!!rowSums(abs(CAZy_SingleDom[-c(1:4)])),]

## Make Presence/Absence Dataframe
#CAZy_PresAbs <- CAZy_ALL
#CAZy_PresAbs[,3:ncol(CAZy_PresAbs)] <- ifelse(CAZy_PresAbs[,3:ncol(CAZy_PresAbs)] > 0 ,1,0)

## Should also consider 
## making a normalized to 
## gene number dataframe...need sample info file


## Now we will calculate the number of Substrate Class hits per sample

Substrate_Classes <- read.delim("~/Desktop/testbed/Substrate_Class_Template.txt", header = TRUE) # Change to biotite dir

Substrate_Class_Hits <- data.frame()

for (i in 1:nrow(Substrate_Classes)){
  Class_Subset <- CAZy_SingleDom[grepl(Substrate_Classes[i,1], CAZy_SingleDom$Substrate_Class),]
  temp_row <- colSums(Class_Subset[5:ncol(Class_Subset)])
  Substrate_Class_Hits <- rbind(Substrate_Class_Hits, temp_row)
}


## Plot total number of proteins in genome set acting on a substrate class
Total_Class_Counts <- data.frame(Counts = rowSums(Substrate_Class_Hits))
Total_Class_Counts <- data.frame(Substrate_Class = Substrate_Classes$Substrate_Class, Counts = Total_Class_Counts$Counts)

ggplot(Total_Class_Counts, aes(x = reorder(Substrate_Class, -Counts), y = Counts)) + 
                     geom_bar(stat = "identity", fill = "steelblue4", color = "black") +
                     theme(axis.text.x = element_text(colour = "black", size = 12, angle = 45, hjust = 1),
                           axis.text.y = element_text(colour = "black", size = 12))

ggsave(paste0("~/Desktop/Temp_R_Plots/Emzymes_by_Substrate_Class.pdf"))



## Plot number of proteins normalized to number of genomes acting on a substrate class
Class_Counts_Normalized <- Total_Class_Counts
Class_Counts_Normalized$Counts <- Class_Counts_Normalized$Counts / length(genome_names)

ggplot(Class_Counts_Normalized, aes(x = reorder(Substrate_Class, -Counts), y = Counts)) + 
  geom_bar(stat = "identity", fill = "steelblue4", color = "black") +
  theme(axis.text.x = element_text(colour = "black", size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 12))

ggsave(paste0("~/Desktop/Temp_R_Plots/Emzymes_by_Substrate_Class_GenomeNorm.pdf"))



## Build Heatmap of Substrate_Class_Hits per genome 
Substrate_Class_Hits <- as.matrix(Substrate_Class_Hits)
rownames(Substrate_Class_Hits) <- Substrate_Classes$Substrate_Class
colnames(Substrate_Class_Hits) <- colnames(CAZy_SingleDom[,5:ncol(CAZy_SingleDom)])

pheatmap(Substrate_Class_Hits, scale = "none", clustering_method = "ward.D", margins = (c(10,10)),
         annotation_col = `20cm_Effect`)


  
## Now we will calculate the number of Substrate hits per sample
Substrates <- read.delim("~/Desktop/testbed/Substrate_Template.txt", header = TRUE)

Substrate_Hits <- data.frame()

for (i in 1:nrow(Substrates)){
  Class_Subset <- CAZy_SingleDom[grepl(Substrates[i,1], CAZy_SingleDom$Substrate),]
  temp_row <- colSums(Class_Subset[5:ncol(Class_Subset)])
  Substrate_Hits <- rbind(Substrate_Hits, temp_row)
}

Total_Substrate_Counts <- data.frame(Counts = rowSums(Substrate_Hits))
Total_Substrate_Counts <- data.frame(Substrate = Substrates$Substrate, Counts = Total_Substrate_Counts$Counts)

ggplot(Total_Substrate_Counts, aes(x = reorder(Substrate, -Counts), y = Counts)) + 
  geom_bar(stat = "identity", fill = "steelblue4", color = "black") +
  theme(axis.text.x = element_text(colour = "black", size = 12, angle = 90, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 12))

#ggsave(paste0("~/Desktop/Temp_R_Plots/Emzymes_by_Substrate_Class.pdf"))


## Build table of class counts per genome and build heatmap off table

Substrate_Hits <- as.matrix(Substrate_Hits)
rownames(Substrate_Hits) <- Substrates$Substrate
colnames(Substrate_Hits) <- genome_names


pheatmap(Substrate_Hits, margins = (c(10,10)),clustering_method = "ward.D",
         annotation_col = `20cm_Effect`)

Substrate_Hits_Log <- Substrate_Hits + 1
Substrate_Hits_Log <- apply(Substrate_Hits_Log, 2, log2)

pheatmap(Substrate_Hits_Log, margins = (c(10,10)),clustering_method = "ward.D", cluster_cols = FALSE,
         annotation_col = `20cm_Effect`)


Substrate_Classes <- as.factor(Substrate_Class_Hits$Substrate_Class)
Substrate_Class_Hits <- data.frame(Substrate_Class = Substrate_Classes, Substrate_Class_Hits)
colnames(Substrate_Class_Hits)[3:ncol(Substrate_Class_Hits)] <- genome_names

Class_Counts_Matrix <-  as.matrix(Substrate_Class_Hits[,3:ncol(Substrate_Class_Hits)])
row.names(Class_Counts_Matrix) <- Substrate_Class_Hits$Substrate_Class

pheatmap(Class_Counts_Matrix, scale = "row", margins = c(10,10))




## Calculate and plot the total number of each enzyme type in this dataset
Enzyme_Classes <- c("GH", "GT", "CE", "PL", "AA")
for (i in 1:length(Enzyme_Classes)){
  keep_rows <- grep(Enzyme_Classes[i], CAZy_SingleDom$CAZy_Family)
  functional_subset <- CAZy_SingleDom[keep_rows,]
  function_sums <- data.frame(CAZy_Family = functional_subset$CAZy_Family, 
                              Counts = rowSums(functional_subset[,5:ncol(functional_subset)]))

  ggplot(function_sums, aes(x = reorder(CAZy_Family, -Counts), y = Counts)) +
         geom_bar(stat = "identity", fill = "steelblue4", color = "black") +
         theme(axis.text.x = element_text(colour = "black", size = 10, angle = 90, hjust = 1),
         axis.text.y = element_text(colour = "black", size = 12))
  
  ggsave(paste0("~/Desktop/Temp_R_Plots/",Enzyme_Classes[i],"_Protein_Counts.pdf"))
}


## Build Multiple Heatmaps for number of each CAZy Enzyme Class

for (i in 1:length(Enzyme_Classes)){
  keep_rows <- grep(Enzyme_Classes[i], CAZy_SingleDom$CAZy_Family)
  functional_subset <- CAZy_SingleDom[keep_rows,]
  subset_matrix <- as.matrix(functional_subset[,5:ncol(functional_subset)])
  rownames(subset_matrix) <- functional_subset$CAZy_Family

  pheatmap(subset_matrix,scale = "none",clustering_method = "ward.D2", margins = (c(10,10)))

}

# One big Heatmap for all CAZy
subset_matrix <- as.matrix(CAZy_SingleDom[,5:ncol(CAZy_SingleDom)])
rownames(subset_matrix) <- CAZy_SingleDom$CAZy_Family
pheatmap(subset_matrix, scale = "row",clustering_method = "ward.D", margins = c(10,10), annotation_col = Test_poop)

### Development Notes ###
## Need a way to label and specify significantly altered organisms from some background set 
## Add plots that asses class based on specific substrate in addition to broad substrate group
## 


