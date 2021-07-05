# Testing for enrichment of GO Terms among HGT candidates
# 20190604, cmdoret

library(tidyverse)
library(topGO)
library(viridis)

# Parse CL args
args <- commandArgs(trailingOnly=T)

# Extract gene IDs and GO terms.
annot_tbl <- read_tsv(args[1], col_types='cccciicccccccccccccc')
select_genes <- read_tsv(
    args[2],
    col_names=c("chrom", "start", "end", "strand", "score", "name")
  ) %>%
  pull(name)


out_fig <- args[3]
out_tbl <- args[4]
tmp_mapfile <- paste(dirname(out_fig), 'id2go.tsv', sep='/')

# Write geneID - GO terms mapping to file
annot_tbl %>%
  dplyr::select(GeneID, `GO Terms`) %>%
  mutate(`GO Terms` = gsub('GO_[a-z]*: (GO:[0-9]+) - [^;]*', '\\1', `GO Terms`)) %>%
  mutate(`GO Terms` = str_replace_all(`GO Terms`, ";", ", ")) %>%
  write_tsv(tmp_mapfile)

geneID2GO <- readMappings(tmp_mapfile)
geneNames <- names(geneID2GO)
geneList <- factor(as.integer(geneNames %in% select_genes))
names(geneList) <- geneNames
GOdata <- new(
  "topGOdata",
  description = "Gene at infection-dependent regions in A. castellanii",
  ontology = "BP",
  allGenes = geneList, 
  nodeSize = 10,
  gene2GO = geneID2GO, 
  annot = annFUN.gene2GO
)


fisher.stat <- new(
  "classicCount",
  testStatistic = GOFisherTest, 
  name = "Fisher test"
)
resultFisher <- getSigGroups(GOdata, fisher.stat)
weight.stat <- new(
  "weightCount",
  testStatistic = GOFisherTest, 
  name = "Fisher test with weight algorithm",
  sigRatio = 'ratio'
)
resultWeight <- getSigGroups(GOdata, weight.stat)


# Get 'significant' terms
pvals <- score(resultWeight)
sig <- pvals[pvals<0.05]

termStat(GOdata, names(sig))

# Report top enriched GO terms
allRes <- GenTable(
  GOdata,
  classic = resultFisher,
  weight = resultWeight,
  orderBy = "weight",
  ranksOf = "classic",
  topNodes = 200
)

write_tsv(allRes, out_tbl)

allRes$GO.ID <- paste(allRes$GO.ID, allRes$Term, sep=' - ')
tidy_res = as.tibble(allRes) %>%
  mutate(
    weight = as.numeric(weight),
    classic = as.numeric(classic),
    GO.ID = factor(GO.ID, ordered=T, levels=unique(GO.ID[ordered(weight)]))
  )

pattern = strsplit(basename(out_fig), '_')[[1]][1]
# Visualise top enriched GO terms
svg(out_fig, width=13, height=9)
ggplot(data=tidy_res %>% top_n(30, -weight), 
  aes(x=GO.ID, y=-log10(weight))) + 
  geom_segment(aes(xend=GO.ID, yend=min(-log10(weight))), size=1.1) +
  geom_point(aes(size=Annotated, color=Significant / Expected)) + 
  geom_hline(aes(yintercept=2), lty=2, col='red') +
  theme_minimal() + 
  xlab("") +
  ylab("-log10 pvalue") +
  ggtitle(sprintf("GO enrichment test at infection-dependent %s\n(Fisher exact test, weight algorithm)", pattern)) +
  scale_color_viridis() +
  coord_flip() +
  theme(text=element_text(family="Liberation", size=12))

dev.off()
