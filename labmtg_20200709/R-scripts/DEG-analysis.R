# Differential gene expression analysis using edgeR package
# Goal - perform pairwise comparisons of expression data for Col-0 and hy5 genotypes

# Install and load packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")

library(edgeR)

# Specify paths to data files
justCounts <- "..."
geneLengths <- "..."

# Load gene expression and gene length data
expressionData <- readDGE(paste(justCounts, dir(justCounts), sep=""),
  labels=sapply(strsplit(dir(justCounts), "-"), `[`, 1),
  group=c(...))
expressionData$genes <- read.delim(geneLengths, row.names=1)

# Only keep genes with cpm > 1 in at least two samples
keep <- rowSums(cpm(expressionData) > 1) >= 2
expressionData <- expressionData[keep, ]
expressionData$samples$lib.size <- colSums(expressionData$counts)

# Normalize for library size
expressionData <- calcNormFactors(expressionData)

# Create MDS plot to verify integrity of RNAseq data
plotMDS(expressionData, method="bcv")

# Calculate RPKM values
rpkm <- rpkm(expressionData$counts, gene.length=expressionData$genes$Length, normalized.lib.sizes=TRUE, log=FALSE)