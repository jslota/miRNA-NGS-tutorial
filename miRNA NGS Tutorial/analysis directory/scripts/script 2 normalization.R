#miRNA NGS data analysis
#
#Script 2
#
#Normalizing the data and assessing variation


library(DESeq2)
library(ggplot2)

#load raw data
read_counts <- read.csv("raw data/raw_read_counts.csv", row.names = 1)
sample_info <- read.csv("raw data/sample_info.csv", row.names = 1)
rownames(sample_info) <- sample_info$Sample

#make sure samples are in order
summary(colnames(read_counts)==rownames(sample_info))

#make DEseq data object
dds <- DESeqDataSetFromMatrix(countData = read_counts, colData = sample_info, design = ~CWD_status)
dds <- DESeq(dds)

#plot dispersion estimates to examine normalization
plotDispEsts(dds)

#extract normalized read counts
norm_counts <- varianceStabilizingTransformation(dds)

#make a basic PCA plot
plotPCA(norm_counts, intgroup="CWD_status")

#make a custom PCA plot with ggplot
pcaData <- plotPCA(norm_counts, intgroup=c("Sample", "CWD_status"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=CWD_status, label=Sample)) +
  geom_point(size=3) +
  geom_text(nudge_y = 1) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_color_manual(values = c("navy", "firebrick")) +
  coord_fixed() +
  theme_classic()

#Save normalized read counts for visualization later
write.csv(assay(norm_counts), "raw data/normalized_read_counts.csv")
