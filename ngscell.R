#set the directory
#setwd("/home/shadow/Projects/NGSCell")
setwd("C:/Users/aruna/Documents/Projects/NGSCell")

## load the libraries
library(DESeq2)
library(readxl)
library(DESeq2)
library(GEOquery) 
library(canvasXpress)
library(ggplot2)
library(clinfun)
library(GGally)
library(factoextra)
library(tidyverse)
library(pheatmap)
library(readxl)
library(openxlsx)
library(scales) 
library(viridis)
library(dplyr)
library(reshape2)
library(tidyr)
library(ggdendro)
library(gridExtra)
library(gtable)
library(grid)
library(Glimma)
library(EnhancedVolcano)
library(gage)

#load the data
rawCounts <- read.csv("GSE60424_counts.csv")
sampleData <- read.csv("design.csv")
sampleData_v2 <- sampleData

#count data to a matrix of appropriate for DEseq2 
geneID <-rawCounts$geneID
#geneID <-rownames(rawCounts) 
head(rawCounts)
head(sampleData)
rownames(sampleData) <- sampleData$Run
sampleData$group<- factor(sampleData$Condition, levels = c("Normal", "MultipleSclerosis", "JuvenileDiabetes", "Sepsis", "LateralSclerosis"))
sampleData <- as.data.frame(sampleData)
head(sampleData)
#Checking the columns of the count data in the same order as rows names of the sample mapping, then makeing sure it worked
rawCounts <- rawCounts[,unique(rownames(sampleData))]
all(colnames(rawCounts) %in% rownames(sampleData))
all(colnames(rawCounts) == rownames(sampleData))


#Createing the DEseq2DataSet object
dds<- DESeqDataSetFromMatrix(countData=rawCounts, colData=sampleData, design= ~group)

# Performing pre-filtering of the data
dim(dds)
dds<-dds[rowSums(counts(dds)) > 5, ]
dim(dds)
library(BiocParallel)
#register(MulticoreParam(20))
register(SnowParam(20))
#differential expression
dds <- DESeq(dds)

# Extract differential expression results
res <- results(dds)


#to Store results from dds
write.csv(res,"Newdeseq2Results.csv")
###plot type 1
# Find the top 5 genes based on padj values
top_genes <- rownames(res)[order(res$padj)][1:5]

# Initialize an empty data frame to store the counts
count_data <- data.frame()

# Loop through the top genes and generate count data
for (gene in top_genes) {
  d <- plotCounts(dds, gene = gene, intgroup = "group", returnData = TRUE)
  count_data <- rbind(count_data, d)
}

library("ggplot2")

# Create a plot with all top 5 genes
ggplot(count_data, aes(x = group, y = count)) +
  geom_point(position = position_jitter(w = 0.1, h = 0)) +
  scale_y_log10(breaks = c(25, 100, 400)) +
  facet_wrap(~gene, scales = "free_y")  # This will create separate panels for each gene

###plolt type 2
# Find the top 5 genes based on padj values
top_genes <- rownames(res)[order(res$padj)][1:5]

# Loop through the top genes and generate plotCounts
for (gene in top_genes) {
  plotCounts(dds, gene = gene, intgroup = "group")
  d <- plotCounts(dds, gene = gene, intgroup = "group", returnData = TRUE)
  
  # Create a simple plot for each gene using base R plotting
  plot(d$count, type = 'b', pch = 16, col = 'blue', xlab = 'Group', ylab = 'Count', main = paste("Gene:", gene))
  
  # Add a legend for groups
  legend("topright", legend = levels(d$group), col = 'blue', pch = 16)
}




#plotcounts
plotCounts(dds, gene=which.min(res$padj), intgroup="group")
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="group", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=group, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))





#to Store results from dds
write.csv(res,"deseq2Results.csv")
summary(res)# View summary of results
plotMA(res) 
 
 
# Plot Dispersions:
png("qc-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds, main = "Dispersion plot")
dev.off()

# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)
head(assay(rld))
hist(assay(rld))

# Principal Components Analysis
plotPCA(rld, intgroup = "Condition")


#res <- results(dds, contrast = c("group","control","affected"))
#res <- results(dds, name="group_affected_vs_control")
res <- results(dds)
resultsNames(res)
table(res$padj<0.05)
#FALSE  TRUE 
#6128    23 

## Order by adjusted p-value
res <- res[order(res$padj), ]
## Merge with normalized count data
resdata <- merge(as.data.frame(res), 
                 as.data.frame(counts(dds, normalized = TRUE)), 
                 by = "row.names", 
                 sort = FALSE)
names(resdata)[1] <- "gene"
head(resdata)

## Write results
write.csv(resdata, file = "diffexpr-results.csv", quote = FALSE, row.names = F)

save(res, file = "res.RData")

## Examine plot of p-values
hist(res$pvalue, breaks=50, col="grey")

## MA plot
## Could do with built-in DESeq2 function:
## DESeq2::plotMA(dds, ylim=c(-1,1), cex=1)
## This is Stephen Turner's code:
maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), points(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}
png("diffexpr-maplot.png", 1500, 1000, pointsize=20)
maplot(resdata, main="MA Plot")
dev.off()

## Plots to Examine Results:


## Volcano plot with "significant" genes labeled
library(MASS)
library(calibrate)
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), labs=gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
png("diffexpr-volcanoplot.png", 1200, 1000, pointsize=20)
volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-5, 5))
dev.off()




library("BiocParallel")
register(MulticoreParam(4))

resOrdered <- res[order(res$pvalue),]
summary(res)
sum(res$padj < 0.1, na.rm=TRUE)
res05 <- results(dds, alpha=0.05)
sum(res05$padj < 0.05, na.rm=TRUE)

#MA Plot
plotMA(res, ylim=c(-2,2))

#After calling plotMA, one can use the function identify to interactively detect the row number of individual genes by clicking on the plot. One can then recover the gene identifiers by saving the resulting indices:
#idx <- identify(res$baseMean, res$log2FoldChange)
#rownames(res)[idx]

resNorm <- lfcShrink(dds, coef=2, type="normal")
resAsh <- lfcShrink(dds, coef=2, type="ashr")
resLFC <- lfcShrink(dds, coef=2, type="apeglm")#diid  not worked
summary(res)
plotMA(resLFC, ylim=c(-2,2))
par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm")#did not work
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")

#plot counts
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="group", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=group, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

mcols(res)$description

#Rich visualization and reporting of results
library(Glimma)
glimmaMDS(dds)
glimmaMA(dds)
glimmaVolcano(dds)
glimmaVolcano(res)####
glimmaMD(dds)

#order our results table by the smallest p value:
resOrdered <- res[order(res$pvalue),]
summary(res)
sum(res$padj < 0.05, na.rm=TRUE)

write.csv(as.data.frame(resOrdered), 
          file="group_treated_results.csv")
resSig <- subset(resOrdered, padj < 0.05)
resSig
write.csv(as.data.frame(resSig), 
          file="significant_results.csv")
# Coerce to a data frame
res <- as.data.frame(res)
View(res)



# Examine this data frame
head(res)
res$significant <- ifelse(res$padj < 0.05, "Significant", NA)

ggplot(res, aes(baseMean, log2FoldChange, colour=significant)) + geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=squish) + scale_x_log10() + geom_hline(yintercept = 0, colour="tomato1", size=2) + labs(x="mean of normalized counts", y="log fold change") + scale_colour_manual(name="q-value", values=("Significant"="red"), na.value="grey50") + theme_bw()

# Let's add some more detail
ggplot(res, aes(baseMean, log2FoldChange, colour=padj)) + geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=squish) + scale_x_log10() + geom_hline(yintercept = 0, colour="darkorchid4", size=1, linetype="longdash") + labs(x="mean of normalized counts", y="log fold change") + scale_colour_viridis(direction=-1, trans='sqrt') + theme_bw() + geom_density_2d(colour="black", size=2)


# Transform count data using the variance stablilizing transform
deseq2VST <- vst(dds)

# Convert the DESeq transformed object to a data frame
deseq2VST <- assay(deseq2VST)
deseq2VST <- as.data.frame(deseq2VST)
deseq2VST$Gene <- rownames(deseq2VST)
head(deseq2VST)

# Keep only the significantly differentiated genes where the fold-change was at least 3
# (since my file contains the genes which have the fold-change less than 3)
sigGenes <- rownames(res[res$padj <= .05 & abs(res$log2FoldChange) > 3,])
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
clusterGene <- hclust(distanceGene, method="complete")
clusterSample <- hclust(distanceSample, method="complete")

# Construct a dendogram for samples
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
library(gridExtra)
grid.arrange(sampleDendrogram, heatmap, ncol=1, heights=c(1,5))

# Load in libraries necessary for modifying plots
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
# deseq2ResDF$significant <- ifelse(deseq2ResDF$padj < 0.1 & abs(deseq2ResDF$log2FoldChange) > 1, "Significant", "Not significant")
res$significant <- ifelse(res$padj < .1, "Significant", "Not significant")

# Create volcano plot
ggplot(res, aes(x=log2FoldChange, y=-log10(padj), color=significant)) +
  geom_point(size=1) +
  scale_color_manual(values=c("Not significant"="grey", "Significant"="red")) +
  labs(x="Log2 Fold Change", y="-Log10 adjusted p-value") +
  theme_bw()

rld <- rlog( dds, fitType='mean', blind=TRUE)
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 25)
top25Counts<-assay(rld)[topVarGenes,]
write.csv(top25Counts, file="top25counts1.rld.csv", quote=FALSE)
top25Counts1<-read.table("top25counts1.rld.csv", sep=",",header=TRUE, row.names=1)
pheatmap(top25Counts1,scale = "row")

vsd <- vst(dds, blind = FALSE,fitType = "mean")
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 25)
top25Counts<-assay(vsd)[topVarGenes,]
write.csv(top25Counts, file="top25counts2.rld.csv", quote=FALSE)
top25Counts2<-read.table("top25counts2.rld.csv", sep=",",header=TRUE, row.names=1)

#PLOT HEATMAP USING PHEATMAP
#INCLUDE NEXT LINE IF YOU WANT TO SAVE THE FIGURE IN A FILE
pdf(file="gene.heatmap.pdf")
#heatmap using top 25 varying genes
pheatmap(top25Counts2,scale = "row")

# this gives log2(n + 1)
ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))
ye1vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))

plotDispEsts(dds, main="Dispersion plot")

res <- res[order(res$padj), ]
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)
write.csv(resdata, file="diffexpr-results.csv")
hist(res$pvalue, breaks=50, col="grey", main = "P-value Distribution")
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, topGene, "group")

library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("Sample.Characteristic.disease.","Condition")])
ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))

pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)
meanSdPlot(assay(vsd))
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
sampleDists <- dist(t(assay(vsd)))

library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$group, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
plotPCA(vsd, intgroup=c("Sample.Characteristic.disease.","Condition"))
pcaData <- plotPCA(vsd, intgroup=c("Sample.Characteristic.disease.","Condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=group)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

#could not figure it  outt
EnhancedVolcano(res,
                lab = res,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.01,
                FCcutoff = 2,
                xlim = c(-2, 4),
                col=c('#f6d746', 'black', '#e55c30', '#84206b'), labSize=3.0, colAlpha=.4, 
                title="Mock vs. Infected Immune Response", subtitle = NULL,
                gridlines.major = FALSE, gridlines.minor = FALSE, boxedLabels = TRUE, 
                drawConnectors = TRUE, widthConnectors = 0.25, colConnectors = 'black',
                endsConnectors= 'last')




# VOLCANO PLOT
# Create volcano plot
#with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano Plot", xlim=c(-5,5)))
# Add significance threshold line
#abline(h=-log10(0.05/length(res)), col="red", lwd=2)

# deseq2ResDF$significant <- ifelse(deseq2ResDF$padj < 0.1 & abs(deseq2ResDF$log2FoldChange) > 1, "Significant", "Not significant")
#res$significant <- ifelse(res$padj < .1, "Significant", "Not significant")

# Create volcano plot
#ggplot(res, aes(x=log2FoldChange, y=-log10(padj), color=significant)) +
  #geom_point(size=1) +
  #scale_color_manual(values=c("Not significant"="grey", "Significant"="red")) +
  #labs(x="Log2 Fold Change", y="-Log10 adjusted p-value") +
  #theme_bw()

# install gage
#if (!requireNamespace("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")
#BiocManager::install(c("gage","GO.db","AnnotationDbi","org.Hs.eg.db"))
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

res$symbol = mapIds(org.Hs.eg.db,
                              keys=row.names(res), 
                              column="SYMBOL",
                              keytype="ENSEMBL",
                              multiVals="first")
res$entrez = mapIds(org.Hs.eg.db,
                              keys=row.names(res), 
                              column="ENTREZID",
                              keytype="ENSEMBL",
                              multiVals="first")
res$name =   mapIds(org.Hs.eg.db,
                              keys=row.names(res), 
                              column="GENENAME",
                              keytype="ENSEMBL",
                              multiVals="first")

head(res, 10)
library(dplyr)
library(pathview)
library(gage)
library(gageData)
data(kegg.sets.hs)
data(sigmet.idx.hs)
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
head(kegg.sets.hs, 3)

foldchanges = res$log2FoldChange
names(foldchanges) = res$entrez
head(foldchanges)

keggres = gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)
lapply(keggres, head)

#install.packages("magrittr")
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
library(pathview)
# plot multiple pathways (plots saved to disk and returns a throwaway list object)
tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa"))


data(go.sets.hs)
data(go.subs.hs)
gobpsets = go.sets.hs[go.subs.hs$BP]
gobpres = gage(foldchanges, gsets=gobpsets, same.dir=TRUE)
lapply(gobpres, head)
