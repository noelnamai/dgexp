
setwd("C:/Users/noel.namai/personal/asimov")

library("DESeq2")

sample.files = grep("state", list.files("./results"), value = TRUE)
file.names   = sub("\\.tsv$", "", sample.files)
sample.table = as.data.frame(matrix(unlist(strsplit(file.names, "_")), ncol = 2, byrow = TRUE))
sample.table = cbind(file.names, sample.files, sample.table)
colnames(sample.table) = c("samplename", "filename", "state", "replicate")

ddshtseq = DESeqDataSetFromHTSeqCount(sampleTable = sample.table,
                                      directory = "./results",
                                      design = ~ state)
dds = DESeq(ddshtseq)
res = results(dds, lfcThreshold = log2(1.5), alpha = 0.1)
dexp.genes = as.data.frame(res[!is.na(res$padj) &
                                 res$padj < 0.1 &
                                 res$log2FoldChange > log2(1.5),])

write.table(x = as.data.frame(res),
            file = "./results/genes_results.tsv", 
            sep = "\t",
            na = "0")

write.table(x = dexp.genes,
            file = "./results/dexp_genes.tsv", 
            sep = "\t",
            na = "0")
