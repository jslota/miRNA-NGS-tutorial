#miRNA NGS data analysis
#
#Script 3
#
#Differential expression analysis

library(DESeq2)

#load raw data
read_counts <- read.csv("raw data/raw_read_counts.csv", row.names = 1)
sample_info <- read.csv("raw data/sample_info.csv", row.names = 1)
rownames(sample_info) <- sample_info$Sample

#make DEseq data object
dds <- DESeqDataSetFromMatrix(countData = read_counts, colData = sample_info, design = ~CWD_status)
dds <- DESeq(dds)

#get differential expression results
resultsNames(dds)
res <- results(object = dds, contrast = c("CWD_status", "Pos", "Neg"))

#clean up results file
res <- res[order(res$padj),]
res <- na.omit(res)
res <-as.data.frame(res)

summary(res$padj < 0.05)

#Save differential expression results
if (dir.exists("DE results")==FALSE) { dir.create("DE results") }
write.csv(res, "DE results/CWD_DE_miRNAs.csv")
