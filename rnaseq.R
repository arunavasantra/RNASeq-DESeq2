####_____Set the Directory____#####
setwd("C:/Users/aruna/Desktop/Bulk RNA seq")

#####____Load The Libraries____#####
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

#########______Load The Data____#####

#load the rawcountsdata
rawCounts <- read.delim("GSE185372_RawCounts.csv",header=TRUE, sep=",")
#load the design file
sampleData <- read.delim("design.csv",header=TRUE, sep=",")

#checking rows and columns Length
all( sampleData$Run == colnames(rawCounts$counts) )

#rawcount data to a matrix of appropriate for DEseq2 
geneID <- rawCounts$Gene_ID
sampleIndex <- grepl("GSM\\d+", colnames(rawCounts)) #grabs the gene name of GSM_IDS
rawCounts <- as.matrix(rawCounts[,sampleIndex])
rownames(rawCounts) <-geneID 
head(rawCounts)
rawCounts <- as.data.frame(rawCounts)
head(rawCounts)

head(sampleData)
rownames(sampleData) <- sampleData$Run
keep <- c("Disease_group", "Phenotype")
sampleData <- sampleData[,keep]
#Renaming The column names -not necessary
colnames(sampleData) <- c("Disease", "Condition")
sampleData$Condition <- factor(sampleData$Condition)
head(sampleData)
sampleData <- as.data.frame(sampleData)
##Checking the columns of the count data in the same order as rows names of the sample mapping, then makeing sure it worked
rawCounts <- rawCounts[,unique(rownames(sampleData))]
all(colnames(rawCounts) %in% rownames(sampleData))
all(colnames(rawCounts) == rownames(sampleData))
## rename the Condition
rename_tissues <- function(x){
  x <- switch(as.character(x), "Healthy"="healthy_subject","Active_TB"= "diseased_state" ,"Latent_TB"="diseased_state")
  return(x)
}
sampleData$Condition <- unlist(lapply(sampleData$Condition, rename_tissues))
#Order the Condition so that it is sensible and make sure the control sample is first: normal sample -> disease sample 
sampleData$Condition <- factor(sampleData$Condition, levels=c("healthy_subject", "diseased_state"))

#################__________Creating DESEQ2 Object__________##############

#creating matrix
dds <-DESeqDataSetFromMatrix(countData = rawCounts,colData = sampleData,design = ~ Condition)
dim(dds)
dim(dds[rowSums(counts(dds)) > 5, ])
# Performing pre-filtering of the data
dds <- dds[rowSums(counts(dds)) > 5, ]

library(BiocParallel)
register(MulticoreParam(5))
#differential expression
dds <- DESeq(dds, modelMatrixType = T, full = design(dds))
# Extract differential expression results
deseq2Results <- results(dds)


# Check if pvalue column exists
if (!"pvalue" %in% colnames(deseq2Results)) {
  stop("pvalue column not found in the data. Please ensure it exists.")
}

# Split the row names (assuming they're in the format "ensembleID|geneName")
split_names <- strsplit(rownames(deseq2Results), "\\|")

# Extract ensemble IDs and gene names into separate columns
deseq2Results$ENSEMBLE_ID <- sapply(split_names, `[[`, 1)
deseq2Results$Gene_ID <- sapply(split_names, `[[`, 2)

# Add log10pvalue column (assuming pvalue is numeric)
deseq2Results$log10pvalue <- -log10(deseq2Results$pvalue)

# Check the modified table
head(deseq2Results)


#to Store results from dds
write.csv(deseq2Results,"deseq2Results.csv")

summary(deseq2Results)# View summary of results

#################___Visualization__#################

#MA Plot
plotMA(deseq2Results) 

# Filter by significance threshold (adjust p-value cutoff)
significant <- subset(deseq2Results, pvalue < 0.05)

# Explore significant genes and expression patterns
head(significant)

#Dispersion Estimates
png("qc-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds, main = "Dispersion plot")
dev.off()





