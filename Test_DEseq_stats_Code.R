library(DESeq2)



CAZy_Hits <- CAZy_SingleDom[,5:ncol(CAZy_SingleDom)]
rownames(CAZy_Hits) <- CAZy_SingleDom$CAZy_Family
Sample_Metadata <- data.frame(Sample_name = Genome_metadata$genome_name, Factor = Genome_metadata$exp_plot)


countData <- as.matrix(CAZy_Hits) # Use col1 as rownames
colData <- as.data.frame(Sample_Metadata)


## Build DEseq Object
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ Factor)
print(dds)

dds_out <- DESeq(dds, fitType = "local")
sizeFactors(dds_out)
plotDispEsts(dds_out)
rld <- rlog(dds)
vsd <- varianceStabilizingTransformation(dds)

library("vsn")
notAllZero <- (rowSums(counts(dds))>0)
meanSdPlot(log2(counts(dds_out,normalized=TRUE)[notAllZero,] + 1))
meanSdPlot(assay(rld[notAllZero,]))
meanSdPlot(assay(vsd[notAllZero,]))


## Generate PCA and other Diagnostic Plots
rld <- rlog(dds_out)
plotPCA(rld, intgroup=c("Factor"))


resultsNames(dds_out)
DEseq_results <- results(dds_out)


summary(DEseq_results)

DEseq_significant <- subset(DEseq_results, padj < 0.1)
DEseq_significant_out <- as.data.frame(DEseq_significant)


DEseq_heat_obj <- assay(dds_out)[DEseq_significant@rownames,]
df <- as.data.frame(colData(dds_out)[,c("Factor")])

dcols <- vegdist(t(log2(DEseq_heat_obj+1)), method = "bray")

plot(hclust(dcols, method = "ward.D"), hang = -1)

test <- data.frame(cmdscale(dcols), Factor = Genome_metadata$exp_plot)

ggplot(test, aes(x = X1, y = X2, color = Factor)) + 
  geom_point(size = 3) +
  scale_color_manual(values = c("steelblue", "firebrick3"))
  
pheatmap(log2(DEseq_heat_obj+1), clustering_distance_rows = drows,
         clustering_distance_cols = dcols,
         clustering_method = "ward.D",
         treeheight_col = 150,
         margins = (c(10,10)),
         annotation_col = heat_anno)



#gene_to_plot <- "PL9"

#d <- plotCounts(dds_out, gene=gene_to_plot, intgroup=c("Factor"),
                returnData=TRUE)
#ggplot(d, aes(x=Factor, y=count, fill=Factor)) + 
#  ggtitle(gene_to_plot) +
#  stat_summary(fun.y="mean", geom="point", size = 4, shape = 21) + 
#  geom_jitter() +
#  ylim(c(0,5))
#  #scale_y_log10() +
#  #scale_x_discrete(limits = c("control20", "treatment20", "control30", "treatment30", "control40", "treatment40")) +
#  #geom_boxplot()  
#  #scale_fill_manual(values = colourList)



