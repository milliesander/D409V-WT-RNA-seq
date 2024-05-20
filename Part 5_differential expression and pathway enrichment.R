# D409V/WT medial septum RNA seq analysis 
# Author: Millie Sander
# Part 5: differential expression analysis 

# package installation: 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("edgeR")
BiocManager::install("sva")
BiocManager::install("DEGreport")
BiocManager::install("apeglm")
BiocManager::install("rrvgo")
BiocManager::install("pheatmap")
BiocManager::install("rafalib")
BiocManager::install("WGCNA")

# loading required packages: 
library(DESeq2)
library(edgeR)
library(pheatmap)
library(RColorBrewer)
library(WGCNA)
library(sva)
library(rafalib)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(EnhancedVolcano)
library(writexl)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(DEGreport)
library(rrvgo)

dds_results <- DESeq(combat_dds2)
results <- results(dds_results)
res_filt <- na.omit(results)

res_tibble <- res_filt %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

# filtering results to retain only significant DEGs in table 

GBA_DEGs <- res_filt[res_filt$pvalue < 0.05,]

# keeping the first 18 digits per gene symbol in order to convert to ensembl ID 
rownames(GBA_DEGs) <- substr(rownames(GBA_DEGs), 1, 18)

DEGs <- rownames(GBA_DEGs)
DEGs <- substr(DEGs, 1, 18)
DEGs <- mapIds(org.Mm.eg.db, keys = DEGs, keytype = "ENSEMBL", column = "SYMBOL")
DEGs <- as.data.frame(DEGs)

GBA_DEGs <- as.data.frame(GBA_DEGs)
DEGs <- merge(GBA_DEGs, DEGs, by = "row.names", all = TRUE)

library(xlsx)
write.xlsx(DEGs, "/lustre/home/ms739/GBA_RNA-seq/DEGs.xlsx")

# visualising results 
# volcano plot with FDR significant genes labelled
res_col <- res_tibble
res_col <- res_col %>% 
  mutate(
    Expression = case_when(log2FoldChange >= 0 & pvalue <= 0.05 ~ "Up-regulated",
                           log2FoldChange <= 0 & pvalue <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged"))
head(res_col)

p1 <- ggplot(res_col, aes(log2FoldChange, -log(pvalue,10))) +
  geom_point(aes(color = Expression), size = 2/5) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*" p value")) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) 

p2 <- p1 + geom_label_repel(data = FDR_genes,
                   mapping = aes(log2FoldChange, -log(pvalue,10), label = symbol),
                   size = 2)

pdf(file = "/lustre/home/ms739/GBA_RNA-seq/ggplot_volcano_labelled.pdf", height = 4, width = 5)
p2
dev.off()

# extracting FDR significant genes to plot as violin plots (using external software, OriginPro)

FDR_genes <- (res_tb[res_tb$padj < 0.05,])
FDR_genes <- as.data.frame(FDR_genes)
rownames(FDR_genes) <- FDR_genes$gene
FDR_genes$gene <- substr(rownames(FDR_genes), 1, 18)
FDR_genes <- FDR_genes[order(FDR_genes$padj),]

# adding symbol name 
FDR_genes$symbol <- mapIds(org.Mm.eg.db, keys = FDR_genes$gene, keytype = "ENSEMBL", column = "SYMBOL")

counts <- norm.cts
counts$gene <- substr(rownames(norm.cts), 1, 18)
FDR_counts <- merge(FDR_genes, counts, by = "gene")

write_xlsx(FDR_counts, "/lustre/home/ms739/GBA_RNA-seq/FDR_counts.xlsx")


# pathway enrichment 
# create dataframes from the DE results table 
GO_data <- as.data.frame(res_tibble)
GO_data_shrunken <- as.data.frame(res_tb)

# set genes as the rownames
rownames(GO_data) <- GO_data$gene
rownames(GO_data_shrunken) <- GO_data_shrunken$gene

# remove the previous gene column 
GO_data <- subset(GO_data, select = -gene)
GO_data_shrunken <- subset(GO_data_shrunken, select = -gene)

## for enrichGO, we need to convert to entrez ID
## to do this, we need to remove the decimal point from ensembl IDs
GO_list <- row.names(GO_data)
GO_list <- substr(GO_list, 1, 18)
GO_data_entrez <- mapIds(org.Mm.eg.db, keys = GO_list, keytype = "ENSEMBL", column = "ENTREZID")

GO_results <- enrichGO(gene = GO_data_entrez,
                       OrgDb = "org.Mm.eg.db",
                       ont = "BP")

# binding GO terms using rrvgo 

simMatrix <- calculateSimMatrix(GO_results$ID,
                                orgdb="org.Mm.eg.db",
                                ont="BP",
                                method="Rel")

# In calculateSimMatrix(GO_results$ID, orgdb = "org.Mm.eg.db", ont = "BP",  :
#   Removed 14 terms that were not found in orgdb for BP

scores <- setNames(-log10(GO_results$qvalue), GO_results$ID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Mm.eg.db")

reduced_GO <- ggplot(reducedTerms, aes(x = as.factor(Description), y = -log10(qvalue)))+
  geom_point(fill = "#c0ff80",colour = "#2d5900", shape = 21, size = 5)+
  coord_flip()+
  theme_bw()

pdf(file = "/lustre/home/ms739/GBA_RNA-seq/WGCNA/GO_reduced_top10.pdf", height = 6, width = 6)
reduced_GO
dev.off()
