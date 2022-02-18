#miRNA NGS data analysis
#
#Script 4
#
#Common data visualizations

library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(dplyr)

#The volcano plot
#Load differential expression results from 150 dpi
res <- read.csv("DE results/CWD_DE_miRNAs.csv")

#a basic volcano plot
ggplot(res, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point()

#a nicer volcano plot
ggplot(res, aes(x=log2FoldChange, y=-log10(padj), color=stat)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color="grey50") +
  geom_vline(xintercept = 0.85, linetype="dashed", color="grey50") +
  geom_vline(xintercept = -0.85, linetype="dashed", color="grey50") +
  scale_color_gradient2(low = "navy", high = "firebrick", mid="grey95", midpoint = 0) +
  theme_classic()


#The heatmap
#Identify DE miRNAs
miRNAs <- read.csv("DE results/CWD_DE_miRNAs.csv") %>%  filter(padj < 0.05, abs(log2FoldChange) > 0.85, baseMean > 15) %>% pull(X)

#We will use normalized read-counts to calculate z-scores
zscores <- read.csv("raw data/normalized_read_counts.csv", row.names = 1)
zscores <- as.matrix(zscores[miRNAs,])
zscores <- (zscores-rowMeans(zscores))/matrixStats::rowSds(zscores)

#basic heatmap with hierachichal clustering
pheatmap(zscores)

#nicer heatmap

#specify additional variables required by pheatmap
plot_colors <- rev(colorRampPalette(brewer.pal(11,"PuOr"))(100))
column_annotation <- read.csv("raw data/sample_info.csv", row.names = 1)
column_annotation <- data.frame(row.names = column_annotation$Sample,
                                CWD_status = column_annotation$CWD_status)
annotation_colors <- list(`CWD_status`=c(`Pos`="firebrick", `Neg`="navy"))

pheatmap(zscores, color = plot_colors, annotation_col = column_annotation, annotation_colors = annotation_colors,
         treeheight_row = 25, treeheight_col = 25, border_color = FALSE)
