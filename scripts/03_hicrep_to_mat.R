# Make a correlation matrix from correlation coefficients ofa all HiCrep runs
# cmdoret, 20191015

library(readr)
library(dplyr)
library(ggplot2)
library(tibble)
library(pheatmap)

args <- commandArgs(trailingOnly=T)

# 0: Load data
corr_tbl <- read_tsv(args[1], col_names=F)
samples <- read_tsv(args[2], comment="#") %>%
  select(-protocol)
heatmap <- args[3]

# 1: Format data
colnames(corr_tbl) <- c("lib1", "lib2", "scc")

# Make matrix symmetric
rev_tbl <- corr_tbl %>%
  rename(lib3=lib1, lib1=lib2) %>%
  rename(lib2=lib3)
  
corr_tbl <- corr_tbl %>%
  bind_rows(rev_tbl)

# Add self comparisons
corr_tbl <- corr_tbl %>% bind_rows(
    tibble(lib1 = unique(corr_tbl$lib1), lib2 = unique(corr_tbl$lib1), scc = rep(1))
  )

corr_samples <- corr_tbl %>%
  inner_join(samples, by=c("lib1" = "library")) %>%
  inner_join(samples, by=c("lib2" = "library"), suffix=c('_1', '_2')) %>%
  mutate(cond_combo = paste(condition_1, "-", condition_2)) %>%
  arrange(cond_combo)


#t.test(data=corr_samples, "scc ~ cond_combo")

#ggplot(data=corr_samples, aes(x=cond_combo, y=scc)) + geom_point()
#ggplot(data = corr_samples, aes(x=lib1, y=lib2, fill=scc)) + 
#  geom_tile()
sample_annot <- column_to_rownames(samples, 'library') %>%
  select("condition", "infection_time")

pdf(heatmap, width=15, height=12)
pheatmap(xtabs(data=corr_samples, "scc ~ lib1 + lib2"), annotation_col=sample_annot)
dev.off()
