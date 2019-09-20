#!/bin/env Rscript
# Differential domain boundaries analysis with diffHic.

library(optparse)
library(diffHic)
library(edgeR)
library(yaml)
library(readr)
option_list = list(
  make_option(c("-s", "--samples"), type="character", default="samples.tsv", help="samples metatdata file path [default= %default]", metavar="character"),
  make_option(c("-c", "--config"), type="character", default="config.yaml", help="analysis configuration file path [default= %default]", metavar="character"),
  make_option(c("-i", "--indir"), type="character", help="path to the directory containing input h5 files", metavar="character"),
  make_option(c("-o", "--outfile"), type="character", help="path to the output file", metavar="character")) 
 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser);

# Load configuration parameters
conf <- yaml.load_file(opt$config)
samples <- read_tsv(opt$samples)

# Build list of input and output files
input <- paste(file.path(opt$indir, samples$library), ".h5", sep="")

# Load genome info
mm.frag <- cutGenome(conf$reference, c("GATC", "GANTC"), c(4, 3))
mm.param <- pairParam(mm.frag)
bin.size <- conf$contact_maps$comp_res

# Define analysis parameters
finder <- domainDirections(input, mm.param, width=bin.size, span=10)
all.counts <- cbind(assay(finder, "up"), assay(finder, "down"))
totals <- unname(sapply(input, function(x){totalCounts(x, mm.param)}))
ydom <- DGEList(all.counts, lib.size=rep(totals, 2))

# Explain experimental design
Condition <- factor(samples$condition)
Sample <- factor(seq_along(input))
Condition <- rep(Condition, 2)
Sample <- rep(Sample, 2)
Direction <- rep(c("Up", "Down"), each=length(input))
design <- model.matrix(~0 + Sample + Direction:Condition)
design <- design[,!grepl("DirectionDown", colnames(design))]
colnames(design) <- sub("DirectionUp:", "", colnames(design))

# Fit model
ab <- aveLogCPM(ydom)
keep <- ab > 0
ydom <- ydom[keep,]
#summary(keep)
cur.regions <- rowRanges(finder)[keep,]
ydom <- estimateDisp(ydom, design)
fitdom <- glmQLFit(ydom, design, robust=TRUE)

# Exploratory plots
base_out <- tools::file_path_sans_ext(opt$outfile)
pdf(paste0(base_out, "_BCV.pdf"))
plotBCV(ydom)
dev.off()
pdf(paste0(base_out, "_QLDisp.pdf"))
plotQLDisp(fitdom)
dev.off()

# Statistical significance
con <- makeContrasts(Conditioninfected - Conditionuninfected, levels=design)
resdom <- glmQLFTest(fitdom, contrast=con)
#topTags(resdom)
output <- data.frame(as.data.frame(cur.regions)[,1:3],
                     infected=fitdom$coefficients[,"Conditioninfected"]/log(2),
                     uninfected=fitdom$coefficients[,"Conditionuninfected"]/log(2),
                     resdom$table)

# Correct for multiple testing
output$FDR <- p.adjust(resdom$table$PValue, method="BH")

# Transform data to output format
o <- order(output$PValue)
output <- output[o,]
write_tsv(output, opt$outfile)
