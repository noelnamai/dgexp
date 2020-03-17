"

Usage:
  dgexp.R [--threads=<threads>] [--results-dir=<dir>] [--alpha=<alpha>] [--log2fc-threshold=<log2fc>]
  dgexp.R (-h | --help)

Options:
  -h --help                   Show help and exit
  --threads=<threads>         Number of workers to used to run DESeq2 [default: 1]
  --alpha=<alpha>             The adjusted p value cutoff [default: 0.01]
  --log2fc-threshold=<log2fc> Threshold for constructing Wald tests of significance [default: 1.50]
  --results-dir=<dir>         Directory where htseq-count filenames are located [default: '.']

" -> doc

########################################################################################################

suppressMessages(suppressWarnings(library("DESeq2")))
suppressMessages(suppressWarnings(library("docopt")))
suppressMessages(suppressWarnings(library("BiocParallel")))

########################################################################################################

args = docopt(doc)

results.dir = args[["--results-dir"]]
alpha       = as.double(args[["--alpha"]])
log2fc      = as.double(args[["--log2fc-threshold"]])
threads     = as.double(args[["--threads"]])

########################################################################################################

sample.files = grep("state", list.files(results.dir), value = TRUE)
file.names   = sub("\\.tsv$", "", sample.files)
sample.table = as.data.frame(matrix(unlist(strsplit(file.names, "_")), ncol = 2, byrow = TRUE))
sample.table = cbind(file.names, sample.files, sample.table)
colnames(sample.table) = c("samplename", "filename", "state", "replicate")

########################################################################################################

ddshtseq = DESeqDataSetFromHTSeqCount(sampleTable = sample.table,
                                      directory = results.dir,
                                      design = ~ state)

dds = DESeq(object = ddshtseq, parallel = TRUE)
res = results(dds, lfcThreshold = log2(log2fc), alpha = alpha)
dexp.genes = as.data.frame(res[!is.na(res$padj) &
                                 res$padj < alpha &
                                 res$log2FoldChange > log2(log2fc), ])

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