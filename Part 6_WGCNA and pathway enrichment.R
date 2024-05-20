# D409V/WT medial septum RNA seq analysis 
# Author: Millie Sander
# Part 6: WGCNA

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

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)

sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# network construction
net = blockwiseModules(datExpr0, power = 18,
TOMType = "signed", minModuleSize = 30, maxBlockSize = 20000,
reassignThreshold = 0, mergeCutHeight = 0.25,
numericLabels = TRUE, pamRespectsDendro = FALSE,
saveTOMs = TRUE,
saveTOMFileBase = "GBA_TOM",
verbose = 3)

# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath

pdf(file = "/lustre/home/ms739/GBA_RNA-seq/WGCNA/moduleTree.pdf", height = 9, width = 12)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
"Module colors",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
dev.off()

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]

# relating modules to external traits 
# Define numbers of genes and samples
nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "(",
signif(moduleTraitPvalue, 1), ")", sep = " ");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 10, 3, 3))

# Display the correlation values within a heatmap plot

pdf(file = "/lustre/home/ms739/GBA_RNA-seq/WGCNA/traitHeatMap.pdf", height = 16, width = 6)
par(mar = c(6, 10, 3, 3))

labeledHeatmap(Matrix = moduleTraitCor,
xLabels = names(datTraits),
yLabels = names(MEs),
ySymbols = names(MEs),
colorLabels = FALSE,
colors = greenWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.8,
cex.lab.y = 1,
zlim = c(-1,1),
main = paste("Module-trait relationships"))
dev.off()

# Define variable weight containing the weight column of datTrait
genotype = as.data.frame(datTraits$condition)
names(genotype) = "genotype"

# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
geneTraitSignificance = as.data.frame(cor(datExpr0, genotype, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(genotype), sep="")
names(GSPvalue) = paste("p.GS.", names(genotype), sep="")

# MM vs GS scatter plots per significant module
module = "plum2"
column = match(module, modNames)
moduleGenes = moduleColors==module

pdf(file = "/lustre/home/ms739/GBA_RNA-seq/WGCNA/plum2_GS_MM.pdf", height = 5, width = 5)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
abs(geneTraitSignificance[moduleGenes, 1]),
xlab = paste("Module Membership in", module, "module"),
ylab = "Gene significance for genotype",
main = paste("Module membership vs. gene significance\n"),
cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "deeppink3")
# abline(h = 0.6)
# abline(v = 0.95)
dev.off()

module = "red"
column = match(module, modNames)
moduleGenes = moduleColors==module

pdf(file = "/lustre/home/ms739/GBA_RNA-seq/WGCNA/red_GS_MM.pdf", height = 5, width = 5)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
abs(geneTraitSignificance[moduleGenes, 1]),
xlab = paste("Module Membership in", module, "module"),
ylab = "Gene significance for genotype",
main = paste("Module membership vs. gene significance\n"),
cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
# abline(h = 0.7)
# abline(v = 0.9)
dev.off()

module = "yellowgreen"
column = match(module, modNames)
moduleGenes = moduleColors==module

pdf(file = "/lustre/home/ms739/GBA_RNA-seq/WGCNA/yellowgreen_GS_MM.pdf", height = 5, width = 5)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
abs(geneTraitSignificance[moduleGenes, 1]),
xlab = paste("Module Membership in", module, "module"),
ylab = "Gene significance for genotype",
main = paste("Module membership vs. gene significance\n"),
cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
# abline(h = 0.7)
# abline(v = 0.95)
dev.off()

# summary output of module analysis 
yellowgreen_genes <- names(datExpr0)[moduleColors=="yellowgreen"]
yellowgreen_genes <- as.data.frame(yellowgreen_genes)
yellowgreen_genes <- na.omit(yellowgreen_genes)
yellowgreen_genes <- as.character(yellowgreen_genes$yellowgreen_genes)

red_genes <- names(datExpr0)[moduleColors=="red"]
red_genes <- as.data.frame(red_genes)
red_genes <- na.omit(red_genes)
red_genes <- as.character(red_genes$red_genes)

plum2_genes <- names(datExpr0)[moduleColors=="plum2"]
plum2_genes <- as.data.frame(plum2_genes)
plum2_genes <- na.omit(plum2_genes)
plum2_genes <- as.character(plum2_genes$plum2_genes)

# identifying hub genes
geneTraitSignificance <- as.data.frame(geneTraitSignificance)
GS_MM <- merge(geneModuleMembership, geneTraitSignificance, by = rownames)

# red module
GS_MM_red <- GS_MM[,c("Row.names", "GS.genotype", "MMred")]

#yellowgreen module
GS_MM_yellowgreen <- GS_MM[,c("Row.names", "GS.genotype", "MMyellowgreen")]

#plum2 module
GS_MM_plum2 <- GS_MM[,c("Row.names", "GS.genotype", "MMplum2")]

# preparing data for network analysis in cytoscape
load("/lustre/home/ms739/GBA_RNA-seq/GBA_TOM-block.1.RData")
TOM <- as.matrix(TOM)

# Select modules
modules = c("plum2") # edit this per module 
# Select module probes
probes = names(datExpr0)
inModule = is.finite(match(moduleColors, modules))
modProbes = probes[inModule]
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]

modProbes <- as.data.frame(modProbes)
modProbes2 <- na.omit(modProbes)

modProbes2$modProbes2 <- substr(modProbes2$modProbes, 1, 18)
modProbes2$gene_id <- mapIds(org.Mm.eg.db, keys = modProbes2$modProbes2, keytype = "ENSEMBL", column = "SYMBOL")

dimnames(modTOM) = list(modProbes2$modProbes, modProbes2$modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                              edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                              nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                              weighted = TRUE,
                              threshold = 0.02,
                              nodeNames = modProbes2$modProbes,
                              altNodeNames = modProbes2$gene_id,
                              nodeAttr = moduleColors[inModule])

# pathway enrichment analysis per significant module

GO_yellowgreen <- substr(yellowgreen_genes, 1, 18)
GO_yellowgreen_entrez <- mapIds(org.Mm.eg.db, keys = GO_yellowgreen, keytype = "ENSEMBL", column = "ENTREZID")
GO_yellowgreen_results <- enrichGO(gene = GO_yellowgreen_entrez,
                       OrgDb = "org.Mm.eg.db",
                       ont = "BP")
# zero enriched terms found

GO_red <- substr(red_genes, 1, 18)
GO_red_entrez <- mapIds(org.Mm.eg.db, keys = GO_red, keytype = "ENSEMBL", column = "ENTREZID")
GO_red_results <- enrichGO(gene = GO_red_entrez,
                       OrgDb = "org.Mm.eg.db",
                       ont = "BP")
# 99 enriched terms found

GO_plum2 <- substr(plum2_genes, 1, 18)
GO_plum2_entrez <- mapIds(org.Mm.eg.db, keys = GO_plum2, keytype = "ENSEMBL", column = "ENTREZID")
GO_plum2_results <- enrichGO(gene = GO_plum2_entrez,
                       OrgDb = "org.Mm.eg.db",
                       ont = "BP")
# zero enriched terms found

# plotting top 20 terms for the red module 
top20_red <- head(GO_red_results, 20)
red_GO <- ggplot(top20_red, aes(x = as.factor(Description), y = -log10(qvalue)))+
  geom_point(fill = "#90ee90",colour = "#228b22", shape = 21, size = 5)+
  coord_flip()+
  theme_bw()

pdf(file = "/lustre/home/ms739/GBA_RNA-seq/WGCNA/GO_red_top20.pdf", height = 6, width = 6)
red_GO
dev.off()

# exporting WGCNA modules for processing in cytoscape 
# Select modules
modules = c("red")
# Select module probes
probes = names(datExpr)
inModule = is.finite(match(moduleColors, modules))
modProbes = probes[inModule]
modGenes = annot$gene_name[match(modProbes, annot$gene_id)]
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                              edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                              nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                              weighted = TRUE,
                              threshold = 0.02,
                              nodeNames = modProbes,
                              altNodeNames = modGenes,
                              nodeAttr = moduleColors[inModule])

 
