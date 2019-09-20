# Download list of ENSEMBL gene ids with positions, gene symbols and associated GO terms using biomaRt
# cmdoret, 20190905
library(biomaRt)
library(tidyverse)

args <- commandArgs(trailingOnly=T)
outfile <- args[1]
# select mart and data set
bm <- useMart("ensembl")
bm <- useDataset("mmusculus_gene_ensembl", mart=bm)

# Get ensembl gene ids and GO terms
EG2GO <- getBM(mart=bm, attributes=c('ensembl_gene_id', 'chromosome_name', 'start_position', 'end_position', 'external_gene_name', 'entrezgene_id', 'go_id'))
EG2GO <- EG2GO %>%
  as_tibble %>%
  filter(str_detect(chromosome_name, '^[0-9]*$')) %>% 
  mutate(chromosome_name=paste0('chr', as.numeric(chromosome_name))) %>%
  distinct(ensembl_gene_id, .keep_all=T) %>%
  arrange(chromosome_name, start_position)

# Remove blank entries
#EG2GO <- EG2GO[EG2GO$go_id != '',]

write.table(EG2GO, outfile, row.names=F, sep='\t', quote=F)
