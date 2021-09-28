# Load salmon output from nf-core/rnaseq and compute differential expression
# statistics via DESeq2.
# cmdoret, 20210708
library(DESeq2)
library(apeglm)

# Load RData file with dds experiment from nf-core/rnaseq
args <- commandArgs(trailingOnly=T)
load(args[1])
outfile <- args[2]

### LOAD DATA ###

# Ensure we expression relative to uninfected
print(resultsNames(dds))
dds$condition <- relevel(dds$condition, ref="uninfected")
dds <- DESeq(dds)
print(resultsNames(dds))
# Pre-filtering genes with extremely low counts (<10 reads in total)
res <- results(dds)
# Computing log fold change with shrinkage estimator from apeglm
resLFCshrink <- lfcShrink(dds, coef="condition_infected_vs_uninfected", res=res, type="apeglm")

### SAVE OUTPUT ###
# Writing diff expression stats table
resLFCshrink_df <- as.data.frame(resLFCshrink)
write.table(resLFCshrink_df, file=outfile, sep='\t', quote=F, row.names=T)
