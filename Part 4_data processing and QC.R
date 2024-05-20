# D409V/WT medial septum RNA seq analysis 
# Author: Millie Sander
# Part 4: data processing and QC

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

# new function: pca x trait correlation plot 
pca.correlation <- function(pca_score, pheno_df, n_PC, file_name){
  sub_pheno <- pheno_df
  princcomps <- as.data.frame(pca_score[,1:n_PC])
  princcomps_PCA <- merge(sub_pheno, princcomps,  by ="row.names")
  princcomps_PCA <- princcomps_PCA[,-1]
  #### create correlation matrix for p_value
  cor_p <- cor.mtest(princcomps_PCA, conf.level=0.95)
  p.val_df <- cor_p$p
  p.val_df <- p.val_df[1:ncol(sub_pheno),(ncol(sub_pheno)+1):ncol(p.val_df)]
  #### create raw correlation matrix
  cor_PCA <- cor(princcomps_PCA[,(ncol(sub_pheno)+1):ncol(princcomps_PCA)],
                 princcomps_PCA[,1:ncol(sub_pheno)],
                 use = "complete.obs")
  colnames(cor_PCA) <- colnames(sub_pheno)
  FullPCA_Corr <- corrplot(cor_PCA,p.mat = t(p.val_df), addrect = 2,method = 'color',
                           insig = 'label_sig', sig.level = c(0.001, 0.01, 0.05),
                           tl.cex = 0.6, cl.cex = 0.4, cl.offset = 0.2, cl.align.text = "l",
                           pch.cex = 0.6, pch.col = 'grey20')

# PART 1: quality control 
# create a path to state where the data is stored
directory <- "/lustre/home/ms739/GBA_RNA-seq"

# create a path to state where htseq files are stored
htseq_directory <- "/lustre/home/ms739/GBA_RNA-seq/htseq_files"

# load sample table text file 
sampleTable <- read.table(file.path(directory, "sample_table.txt"), header = TRUE)

# creating the DESeq data set 
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = htseq_directory, design= ~ condition)
count_df <- dds@assays@data@listData$counts
meta_df <- dds@colData 

# looking at average counts per sample to compare read depth
avg_9531 <- mean(count_df[,1])
avg_13136 <- mean(count_df[,2])
avg_13137 <- mean(count_df[,3])
avg_13138 <- mean(count_df[,4])
avg_13165 <- mean(count_df[,5])
avg_13166 <- mean(count_df[,6])
avg_13167 <- mean(count_df[,7])
avg_13168 <- mean(count_df[,8])
avg_13169 <- mean(count_df[,9])
avg_13170 <- mean(count_df[,10])
avg_13171 <- mean(count_df[,11])
avg_13174 <- mean(count_df[,12])

avg <- c(avg_9531, avg_13136, avg_13137, avg_13138, avg_13165, avg_13166, avg_13167, avg_13168, avg_13169, avg_13170, avg_13171, avg_13174)
avg <- as.data.frame(avg)
rownames(avg) <- colnames(count_df)

# removing htseq files from the sample table and re-labelling sex/genotype as numeric 
# F = 1; M = 2
# HET = 1; WT = 2

sampleTableNew <- sampleTable
sampleTableNew <- sampleTableNew[,-2]
levels(sampleTableNew$condition) <- c("1", "2")
levels(sampleTableNew$sex) <- c("1", "2")

sampleTable$sampleName <- as.factor(sampleTable$sampleName)
sampleTable$condition <- as.factor(sampleTable$condition)
sampleTable$sex <- as.factor(sampleTable$sex)
sampleTable$TIN <- as.numeric(sampleTable$TIN)
sampleTable$batch <- as.factor(sampleTable$batch)
sampleTable$age <- as.factor(sampleTable$age)

# setting quantiles for the TIN scores (cannot use as co-variate as a continuous variable)

TIN_int <- quantile(sampleTableNew$TIN)
sampleTableNew[which(sampleTableNew$TIN >= TIN_int[[1]] & sampleTableNew$TIN < TIN_int[[2]]), "TIN_quant"] <- 1
sampleTableNew[which(sampleTableNew$TIN >= TIN_int[[2]] & sampleTableNew$TIN < TIN_int[[3]]), "TIN_quant"] <- 2
sampleTableNew[which(sampleTableNew$TIN >= TIN_int[[3]] & sampleTableNew$TIN < TIN_int[[4]]), "TIN_quant"] <- 3
sampleTableNew[which(sampleTableNew$TIN >= TIN_int[[4]] & sampleTableNew$TIN <= TIN_int[[5]]), "TIN_quant"] <- 4
sampleTableNew$TIN_quant <- as.factor(sampleTableNew$TIN_quant)

expression_matrix <- dds@assays@data@listData[["counts"]]

# PCA x trait correlation plot before batch correction 

vsd <- vst(dds, blind = FALSE)
vsd_matrix <- vsd@assays@data@listData[[1]]
pca <- prcomp(t(vsd_matrix))

sampleTableNumeric <- as.data.frame(apply(sampleTableNew, 2, as.numeric))
rownames(sampleTableNumeric) <- sampleTableNumeric[,1]
# setting data as numeric, female = 1, male = 1 and wt = 1, D409V/WT (het) = 2
sampleTableNumeric$sex <- c(1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 1, 2)
sampleTableNumeric$condition <- c(1, 2, 2, 2, 1, 2, 2, 1, 1, 1, 2, 1)

pca_corr <- pca.correlation(pca$x, sampleTableNumeric , 10, "PCA_correlations.png")

pdf(file = "/lustre/home/ms739/GBA_RNA-seq/pca_corr.pdf", height = 3, width = 3)
pca.correlation(pca$x, sampleTableNumeric , 10, "PCA_correlations.png")
dev.off()

# correcting both batch and TIN quantile effects 
combat <- ComBat_seq(counts = expression_matrix, 
	batch = sampleTableNew$batch,
	full_mod = FALSE)

combat2 <- ComBat_seq(counts = combat, 
	batch = sampleTableNew$TIN_quant,
	full_mod = FALSE)

# creating new dataframe using corrected data
rownames(sampleTableNew) <- sampleTableNew$sampleName
all(rownames(sampleTableNew) %in% colnames(combat2))
all(rownames(sampleTableNew) == colnames(combat2))

combat_dds <- DESeqDataSetFromMatrix(countData = combat2,
	colData = sampleTableNew, 
	design = ~ sex + condition)

combat_dds$condition <- relevel(combat_dds$condition, ref = "wt")

# filtering using the group size equivilant to the smallest group in the data (HET and WT both have 6)
# keeping rows (genes) which have at least 10 counts 

smallestGroupSize <- 6
keep <- rowSums(counts(combat_dds) >= 10) >= smallestGroupSize
combat_dds2 <- combat_dds[keep,]

# visiualising results
# estimating library size correction and saving normalised counts
dds2 <- estimateSizeFactors(combat_dds2)
norm.cts <- counts(dds2, normalized=TRUE)

# export normalized counts
norm.cts <- as.data.frame(norm.cts)
write_xlsx(norm.cts, "/lustre/home/ms739/GBA_RNA-seq/norm_counts.xlsx")

sizeFactors <- sizeFactors(dds2)
# normalising using variance stabilising transformation
vsd <- vst(dds2, blind = FALSE)
datExpr0 <- as.data.frame(assay(vsd))
datExpr0 <- t(datExpr0)
names(datExpr0) <- colnames(datExpr0)

# check genes with too many missing samples: 
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

# eucladian distance heatmap
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
colnames(sampleDistMatrix) <- sampleTable$condition
# or without genotype labelled as columns: 
# colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
sampleDist_heatmap <- pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists, col = colors)

pdf(file = "/lustre/home/ms739/GBA_RNA-seq/sampleDist_heatmap.pdf", height = 10, width = 15)
sampleDist_heatmap
dev.off()

# PCA x trait correlation plot 
vsd_matrix <- vsd@assays@data@listData[[1]]
pca <- prcomp(t(vsd_matrix))

sampleTableNumeric <- as.data.frame(apply(sampleTableNew, 2, as.numeric))
rownames(sampleTableNumeric) <- sampleTableNumeric[,1]
# setting data as numeric, female = 1, male = 1 and wt = 1, D409V/WT (het) = 2
sampleTableNumeric$sex <- c(1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 1, 2)
sampleTableNumeric$condition <- c(1, 2, 2, 2, 1, 2, 2, 1, 1, 1, 2, 1)

pca_corr_batchcorrected <- pca.correlation(pca$x, sampleTableNumeric , 10, "PCA_correlations.png")

pdf(file = "/lustre/home/ms739/GBA_RNA-seq/pca_corr_batchcorrected.pdf", height = 3, width = 3)
pca.correlation(pca$x, sampleTableNumeric , 10, "PCA_correlations.png")
dev.off()

# dispersion estimates
pdf(file = "/lustre/home/ms739/GBA_RNA-seq/DispEst.pdf", height = 3, width = 3)
plotDispEsts(dds2)
dev.off()

# sample dendrogram with trait heatmap
# loading trait data 
levels(sampleTableNew$condition) <- c("het", "wt")
levels(sampleTableNew$sex) <- c("f", "m")
traitData <- sampleTableNew

allTraits = traitData[, -c(4)]

# Form a data frame analogous to expression data that will hold the clinical traits.
GBA_samples = rownames(datExpr0)
traitRows = match(GBA_samples, allTraits$sampleName)
datTraits = allTraits[traitRows, -1]
rownames(datTraits) = allTraits[traitRows, 1]

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr0), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry

levels(datTraits$condition) <- c("1", "2")
levels(datTraits$sex) <- c("1", "2")
datTraits$condition <- as.numeric(datTraits$condition)
datTraits$sex <- as.numeric(datTraits$sex)
datTraits$batch <- as.numeric(datTraits$batch)
datTraits$TIN_quant <- as.numeric(datTraits$TIN_quant)
traitColors = numbers2colors(datTraits, signed = FALSE);

# Plot the sample dendrogram and the colors underneath.
traitClustTree <- plotDendroAndColors(sampleTree2, traitColors,
groupLabels = names(datTraits),
main = "Sample dendrogram and trait heatmap")

pdf(file = "/lustre/home/ms739/GBA_RNA-seq/WGCNA/TraitClustTree.pdf", height = 9, width = 12)
plotDendroAndColors(sampleTree2, traitColors,
groupLabels = names(datTraits),
main = "Sample dendrogram and trait heatmap")
dev.off()
