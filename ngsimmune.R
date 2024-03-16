#set the directory
setwd("C:/Users/aruna/Documents/Projects/NGSimmunecells")


## load the libraries
library(DESeq2) 
library(GEOquery)
library(canvasXpress)
library(ggplot2)
library(clinfun)
library(GGally) 
library(factoextra)
library(tidyverse)


#load the data
rawCounts <- read.csv("GSE60424_counts.csv")
sampleData <- read.csv("design.csv")
sampleData_v2 <- sampleData

#count data to a matrix of appropriate for DEseq2 
geneID <- rawCounts$Gene.ID
sampleIndex <- grepl("SRR\\d+", colnames(rawCounts))
rawCounts <- as.matrix(rawCounts[,sampleIndex])
rownames(rawCounts) <- geneID
head(rawCounts)
rawCounts <- as.data.frame(rawCounts)
#sample variable mappings to an appropriate for DESeq2 
head(sampleData)
rownames(sampleData) <- sampleData$Run
keep <- c("Sample.Characteristic.disease.", "Sample.Characteristic.phenotype.")
sampleData <- sampleData[,keep]
colnames(sampleData) <- c("tissueType", "individualID")
sampleData$individualID <- factor(sampleData$individualID)
head(sampleData)
sampleData <- as.data.frame(sampleData)
#Checking the columns of the count data in the same order as rows names of the sample mapping, then makeing sure it worked
rawCounts <- rawCounts[,unique(rownames(sampleData))]
all(colnames(rawCounts) %in% rownames(sampleData))
all(colnames(rawCounts) == rownames(sampleData))
# rename the Individual types
rename_tissues <- function(x){
  x <- switch(as.character(x), "normal"="healthy_subject", "multiple sclerosis"= "diseased_state", "juvenile diabetes"= "diseased_state", "Sepsis"= "diseased_state", "lateral sclerosis"= "diseased_state" )
  return(x)
}
sampleData$individualID <- unlist(lapply(sampleData$individualID, rename_tissues))
#Order the Individual Ids types so that it is sensible and make sure the control sample is first: normal sample -> disease sample 
sampleData$individualID <- factor(sampleData$individualID, levels=c("healthy_subject", "diseased_state"))


#Createing the DEseq2DataSet object
deseq2Data <- DESeqDataSetFromMatrix(countData=rawCounts, colData=sampleData, design= ~individualID)

dim(deseq2Data)
dim(deseq2Data[rowSums(counts(deseq2Data)) > 5, ])
# Performing pre-filtering of the data
deseq2Data <- deseq2Data[rowSums(counts(deseq2Data)) > 5, ]

library(BiocParallel)
register(MulticoreParam(4))
#differential expression
deseq2Data <- DESeq(deseq2Data, modelMatrixType = T, full = design(deseq2Data))
# Extract differential expression results
deseq2Results <- results(deseq2Data)

#to Store results from dds

write.csv(deseq2Results,"deseq2Results.csv")
summary(deseq2Results)# View summary of results
plotMA(deseq2Results) 


# Load libraries
library(ggplot2)
library(scales) # needed for oob parameter
library(viridis)

# Coerce to a data frame
deseq2ResDF <- as.data.frame(deseq2Results)

# Examine this data frame
head(deseq2ResDF)

# Set a boolean column for significance
deseq2ResDF$significant <- ifelse(deseq2ResDF$padj < .1, "Significant", NA)

# Plot the results similar to DEseq2
ggplot(deseq2ResDF, aes(baseMean, log2FoldChange, colour=significant)) + geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=squish) + scale_x_log10() + geom_hline(yintercept = 0, colour="tomato1", size=2) + labs(x="mean of normalized counts", y="log fold change") + scale_colour_manual(name="q-value", values=("Significant"="red"), na.value="grey50") + theme_bw()

# Let's add some more detail
ggplot(deseq2ResDF, aes(baseMean, log2FoldChange, colour=padj)) + geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=squish) + scale_x_log10() + geom_hline(yintercept = 0, colour="darkorchid4", size=1, linetype="longdash") + labs(x="mean of normalized counts", y="log fold change") + scale_colour_viridis(direction=-1, trans='sqrt') + theme_bw() + geom_density_2d(colour="black", size=2)

# Transform count data using the variance stablilizing transform
deseq2VST <- vst(deseq2Data)

# Convert the DESeq transformed object to a data frame
deseq2VST <- assay(deseq2VST)
deseq2VST <- as.data.frame(deseq2VST)
deseq2VST$Gene <- rownames(deseq2VST)
head(deseq2VST)

# Keep only the significantly differentiated genes where the fold-change was at least 3
sigGenes <- rownames(deseq2ResDF[deseq2ResDF$padj <= .05 & abs(deseq2ResDF$log2FoldChange) > 3,])
deseq2VST <- deseq2VST[deseq2VST$Gene %in% sigGenes,]

# Convert the VST counts to long format for ggplot2
library(reshape2)

# First compare wide vs long version
deseq2VST_wide <- deseq2VST
deseq2VST_long <- melt(deseq2VST, id.vars=c("Gene"))

head(deseq2VST_wide)
head(deseq2VST_long)

# Now overwrite our original data frame with the long format
deseq2VST <- melt(deseq2VST, id.vars=c("Gene"))

# Make a heatmap
heatmap <- ggplot(deseq2VST, aes(x=variable, y=Gene, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
heatmap


# Convert the significant genes back to a matrix for clustering
deseq2VSTMatrix <- dcast(deseq2VST, Gene ~ variable)
rownames(deseq2VSTMatrix) <- deseq2VSTMatrix$Gene
deseq2VSTMatrix$Gene <- NULL

# Compute a distance calculation on both dimensions of the matrix
distanceGene <- dist(deseq2VSTMatrix)
distanceSample <- dist(t(deseq2VSTMatrix))

# Cluster based on the distance calculations
clusterGene <- hclust(distanceGene, method="average")
clusterSample <- hclust(distanceSample, method="average")

# Construct a dendogram for samples
#install.packages("ggdendro")
library(ggdendro)
sampleModel <- as.dendrogram(clusterSample)
sampleDendrogramData <- segment(dendro_data(sampleModel, type = "rectangle"))
sampleDendrogram <- ggplot(sampleDendrogramData) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + theme_dendro()

# Re-factor samples for ggplot2
deseq2VST$variable <- factor(deseq2VST$variable, levels=clusterSample$labels[clusterSample$order])

# Construct the heatmap. note that at this point we have only clustered the samples NOT the genes
heatmap <- ggplot(deseq2VST, aes(x=variable, y=Gene, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
heatmap

# Combine the dendrogram and the heatmap
#install.packages("gridExtra")
library(gridExtra)
grid.arrange(sampleDendrogram, heatmap, ncol=1, heights=c(1,5))


# Load in libraries necessary for modifying plots
#install.packages("gtable")
library(gtable)
library(grid)

# Modify the ggplot objects
sampleDendrogram_1 <- sampleDendrogram + scale_x_continuous(expand=c(.0085, .0085)) + scale_y_continuous(expand=c(0, 0))
heatmap_1 <- heatmap + scale_x_discrete(expand=c(0, 0)) + scale_y_discrete(expand=c(0, 0))

# Convert both grid based objects to grobs
sampleDendrogramGrob <- ggplotGrob(sampleDendrogram_1)
heatmapGrob <- ggplotGrob(heatmap_1)

# Check the widths of each grob
sampleDendrogramGrob$widths
heatmapGrob$widths

# Add in the missing columns
sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob, heatmapGrob$widths[7], 6)
sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob, heatmapGrob$widths[8], 7)

# Make sure every width between the two grobs is the same
maxWidth <- unit.pmax(sampleDendrogramGrob$widths, heatmapGrob$widths)
sampleDendrogramGrob$widths <- as.list(maxWidth)
heatmapGrob$widths <- as.list(maxWidth)

# Arrange the grobs into a plot
finalGrob <- arrangeGrob(sampleDendrogramGrob, heatmapGrob, ncol=1, heights=c(2,5))

# Draw the plot
grid.draw(finalGrob)


# VOLCANO PLOT
# Create volcano plot
with(deseq2ResDF, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano Plot", xlim=c(-5,5)))
# Add significance threshold line
abline(h=-log10(0.05/length(res)), col="red", lwd=2)

# deseq2ResDF$significant <- ifelse(deseq2ResDF$padj < 0.1 & abs(deseq2ResDF$log2FoldChange) > 1, "Significant", "Not significant")
deseq2ResDF$significant <- ifelse(deseq2ResDF$padj < .1, "Significant", "Not significant")

# Create volcano plot
ggplot(deseq2ResDF, aes(x=log2FoldChange, y=-log10(padj), color=significant)) +
  geom_point(size=1) +
  scale_color_manual(values=c("Not significant"="grey", "Significant"="red")) +
  labs(x="Log2 Fold Change", y="-Log10 adjusted p-value") +
  theme_bw()

# install gage
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("gage","GO.db","AnnotationDbi","org.Hs.eg.db"))
library(gage)
# set up kegg database
kg.hsa <- kegg.gsets(species="hsa")
kegg.sigmet.gs <- kg.hsa$kg.sets[kg.hsa$sigmet.idx]
kegg.dise.gs <- kg.hsa$kg.sets[kg.hsa$dise.idx]

# set up go database
go.hs <- go.gsets(species="human")
go.bp.gs <- go.hs$go.sets[go.hs$go.subs$BP]
go.mf.gs <- go.hs$go.sets[go.hs$go.subs$MF]
go.cc.gs <- go.hs$go.sets[go.hs$go.subs$CC]

# load in libraries to annotate data
library(AnnotationDbi)
library(org.Hs.eg.db)
columns(org.Hs.eg.db)

deseq2Results$symbol = mapIds(org.Hs.eg.db,
                              keys=row.names(deseq2Results), 
                              column="SYMBOL",
                              keytype="ENSEMBL",
                              multiVals="first")
deseq2Results$entrez = mapIds(org.Hs.eg.db,
                              keys=row.names(deseq2Results), 
                              column="ENTREZID",
                              keytype="ENSEMBL",
                              multiVals="first")
deseq2Results$name =   mapIds(org.Hs.eg.db,
                              keys=row.names(deseq2Results), 
                              column="GENENAME",
                              keytype="ENSEMBL",
                              multiVals="first")

head(deseq2Results, 10)
library(dplyr)
library(pathview)
library(gage)
library(gageData)
data(kegg.sets.hs)
data(sigmet.idx.hs)
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
head(kegg.sets.hs, 3)

foldchanges = deseq2Results$log2FoldChange
names(foldchanges) = deseq2Results$entrez
head(foldchanges)

keggres = gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)
lapply(keggres, head)

install.packages("magrittr")
library(magrittr)
library(dplyr)
# Get the pathways
keggrespathways = data.frame(id=rownames(keggres$greater), keggres$greater) %>% 
  tbl_df() %>% 
  filter(row_number()<=5) %>% 
  .$id %>% 
  as.character()
keggrespathways


# Get the IDs.
keggresids = substr(keggrespathways, start=1, stop=8)
keggresids

keggres$stats


# Define plotting function for applying later
plot_pathway = function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa", new.signature=FALSE)

# plot multiple pathways (plots saved to disk and returns a throwaway list object)
tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa"))


data(go.sets.hs)
data(go.subs.hs)
gobpsets = go.sets.hs[go.subs.hs$BP]
gobpres = gage(foldchanges, gsets=gobpsets, same.dir=TRUE)
lapply(gobpres, head)
