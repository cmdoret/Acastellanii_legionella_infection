#!/bin/env Rscript
# Prepares a hdf5 pairs file compatble with diffHiC from a BAM file.
#if (!requireNamespace("diffHic", quietly = TRUE))
#    install.packages('BiocManager', repos = "http://cran.us.r-project.org")
#    BiocManager::install(ask=F, version="devel")
#    BiocManager::install(ask=F, "diffHic")
library("yaml")
library("optparse")
library("readr")
library("diffHic")

# Parse CL arguments
option_list = list(
  make_option(c("-s", "--samples"), type="character", default="samples.tsv", help="samples metatdata file path [default= %default]", metavar="character"),
  make_option(c("-c", "--config"), type="character", default="config.yaml", help="analysis configuration file path [default= %default]", metavar="character"),
  make_option(c("-i", "--indir"), type="character", help="path to the directory containing input h5 files", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", help="path to the output directory", metavar="character"))
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Load configuration parameters
conf <- yaml.load_file(opt$config)
samples <- read_tsv(opt$samples)

# Build list of input and output files
bam <- paste(file.path(opt$indir, samples$library), ".bam", sep="")
pairs <- paste(file.path(opt$outdir, samples$library), ".h5", sep="")

# Digest genome to define genomic fragments
mm.frag <- cutGenome(conf$reference, c("GATC", "GANTC"), c(4, 3))
mm.param <- pairParam(mm.frag)


# Get bam files into pairs files
for (i in seq_along(bam)){
  prepped <- preparePairs(bam[i], mm.param, file=pairs[i],dedup=TRUE, minq=30)
}

