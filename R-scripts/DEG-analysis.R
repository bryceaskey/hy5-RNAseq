# Differential gene expression analysis using edgeR package
# Goal - perform pairwise comparisons of expression data for Col-0 and hy5 genotypes

# Install and load packages, and set working directory ----
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("edgeR")

library(edgeR)
library(dplyr)

setwd("C:/Users/Bryce/Documents/hy5-RNAseq/")


# Load sample metadata ----
metadata <- read.delim("data/sample-metadata.txt", sep=",")
metadata$genotype <- factor(metadata$genotype)
metadata$treatment <- factor(metadata$treatment)

# For each sample, sum expression data for different transcripts of the same gene together ----
# Use SRA-number in metadata to loop through each transcript-counts file
for(SRA.number in metadata$SRA.number){
  transcript.expression <- read.delim(paste("data/", SRA.number, "_transcript-counts.txt", sep=""), header=FALSE)
  colnames(transcript.expression) <- c("transcript.id", "count")
  
  # Create new column with just gene id
  transcript.expression$gene.id <- sapply(strsplit(transcript.expression$transcript.id, split="[.]"), "[", 1)
  transcript.expression$gene.id <- factor(transcript.expression$gene.id)
  
  # Sum together all transcript counts for each unique gene id to get gene counts
  # ex: AT1G01080.1 and AT1G01080.2 should be summed together as AT1G01080
  gene.expression <- transcript.expression %>%
    group_by(gene.id) %>%
    summarize(count=sum(count))
  
  
}

#read.delim()