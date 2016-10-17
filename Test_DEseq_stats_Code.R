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
dds

dds_out <- DESeq(dds)
sizeFactors(dds_out)
plotDispEsts(dds_out)

## Generate PCA and other Diagnostic Plots
rld <- rlog(dds_out, blind=TRUE)
plotPCA(rld, intgroup=c("Factor"))


resultsNames(dds_out)
DEseq_results <- results(dds_out)


summary(DEseq_results)

DEseq_significant <- subset(DEseq_results, padj < 0.1)
DEseq_significant_out <- as.data.frame(DEseq_significant)


DEseq_heat_obj <- assay(rld)[DEseq_significant@rownames,]
df <- as.data.frame(colData(dds_out)[,c("Factor")])


pheatmap(DEseq_heat_obj,
         #clustering_method = "ward.D2",
         margins = (c(10,10)),
         annotation_col = heat_anno, 
         annotation_colors = ann_cols )

DEseq_heat_obj <- assay(dds_out)[DEseq_significant@rownames,]



library(ropls)
pls_DA_model <- opls(t(DEseq_heat_obj), Genome_metadata$exp_plot)


Scores <- data.frame(class = Genome_metadata$exp_plot, pls_DA_model@scoreMN)

ggplot(Scores, aes(x = p1, y = p2, color = class, fill = class)) + geom_point(size = 3) + 
  stat_ellipse(geom="polygon", alpha = 0.5)


Loadings <- data.frame(Variable = rownames(pls_DA_model@loadingMN), pls_DA_model@loadingMN[,1:2])

ggplot(Loadings, aes(x = p1, y = p2, label = Variable)) + geom_point(alpha = 0.5, size = 3, color = "purple") +
  geom_text(nudge_x = 0.05) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  expand_limits(x = c(-.5,.5), y = c(-.75,.75)) + ggtitle("Loadings")


ggplot() + geom_point(data = Scores, aes(x = p1, y = p2, color = class, fill = class), size = 3) + 
  #stat_ellipse(data = Scores, aes(x = p1, y = p2, color = class, fill = class), geom="polygon", alpha = 0.5) + 
  geom_segment(data = Loadings, aes(x = 0, y = 0, xend = 8*p1, yend = 8*p2), linetype = 5) +
  geom_text(data = Loadings, aes(x = 8*p1, y = 8*p2, label = Variable), nudge_x = -.5)


VIP_scores <- data.frame(Names = rownames(as.data.frame(pls_DA_model@vipVn)), VIP = pls_DA_model@vipVn)

ggplot(VIP_scores, aes(x = reorder(Names, -VIP), y = VIP)) + geom_bar(stat="identity", fill = "steelblue")


for_ggplot <- data.frame(pls_DA_model@scoreMN, Factor = Genome_metadata$exp_plot)

poop <- ggplot(for_ggplot, aes(x = p1, y = p2, color = Factor)) + geom_point()
poop
