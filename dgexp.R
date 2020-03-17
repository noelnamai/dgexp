
setwd("C:/Users/noel.namai/personal/asimov/results/")

########################################################################################################

library("DESeq2")

########################################################################################################

sample.files = grep("state", list.files("."), value = TRUE)
file.names   = sub("\\.tsv$", "", sample.files)
sample.table = as.data.frame(matrix(unlist(strsplit(file.names, "_")), ncol = 2, byrow = TRUE))
sample.table = cbind(file.names, sample.files, sample.table)
colnames(sample.table) = c("samplename", "filename", "state", "replicate")

########################################################################################################

alpha  = 0.05
log2fc = 1.50

ddshtseq = DESeqDataSetFromHTSeqCount(sampleTable = sample.table,
                                      directory = ".",
                                      design = ~ state)
dds = DESeq(ddshtseq)
res = results(dds, lfcThreshold = log2(log2fc), alpha = alpha)
dexp.genes = as.data.frame(res[!is.na(res$padj) &
                                 res$padj < alpha &
                                 res$log2FoldChange > log2(log2fc),])

########################################################################################################

write.table(
  x = as.data.frame(res),
  file = "genes-results.tsv",
  sep = "\t",
  quote = FALSE
)

write.table(
  x = dexp.genes,
  file = "dexp-genes.tsv",
  sep = "\t",
  quote = FALSE
)

########################################################################################################

#ma-plot from base means and log fold changes
png(filename = "ma-plot.png")
plotMA(res, ylim = c(-2, 2))
dev.off()

#plot of normalized counts for a single gene
png(filename = "max-counts-plot.png")
plotCounts(dds, gene = which.max(res$padj), intgroup = "state")
dev.off()

png(filename = "min-counts-plot.png")
plotCounts(dds, gene = which.min(res$padj), intgroup = "state")
dev.off()

########################################################################################################