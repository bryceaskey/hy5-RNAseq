library(dplyr)

annotations <- "D:/bioinformatics/hy5/bedAnnotated/"
peaks <- "D:/bioinformatics/hy5/peaks/"
peaksAnnotated <- "D:/bioinformatics/hy5/peaksAnnotated/"
samples <- c("hy5_1", "hy5_2", "hy5_3", "hy5-VP16_1", "hy5-VP16_2", "hy5-VP16_3", "hy5-SRDX_1", "hy5-SRDX_2")

for(sample in samples){
  print(sample)
  
  geneInfo <- read.delim(paste(annotations, sample, ".txt", sep=""), header=TRUE, sep="", row.names=NULL)
  geneInfo <- geneInfo[, c(4, 15)]
  colnames(geneInfo) <- c("PeakID", "Tags")
  geneInfo$TagID <- sapply(strsplit(geneInfo$Tags, ";"), `[`, 1)
  geneInfo$GeneID <- sapply(strsplit(geneInfo$TagID, "="), `[`, 2)
  geneInfo <- geneInfo[, c(1, 4)]
  geneInfo <- aggregate(GeneID ~ PeakID, data=geneInfo, toString)
  
  peakInfo <- read.delim(paste(peaks, sample, ".txt", sep=""), header=TRUE, sep="\t", row.names=NULL)
  
  merged <- merge(geneInfo, peakInfo, by="PeakID", all.y=TRUE)
  merged <- merged[order(merged$GeneID),]
  finalData <- merged[, c(2:5,11:13)]
  colnames(finalData) <- c("gene_id", "peak_chr", "peak_start", "peak_end", "normalized_tag_count", "fold_change", "p_value")
  write.csv(finalData, file=paste(peaksAnnotated, sample, ".csv", sep=""), row.names=FALSE)
}

