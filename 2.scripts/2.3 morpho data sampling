#===============================================================#
# Head ----------------------------------------------------------------------
setwd("~/Documents/Labis/Dados/") #defining work directory
if (!require(librarian)) install.packages("librarian"); library("librarian")
librarian::shelf(phytools, dplyr, tidyr, purrr, factoextra, vegan, tidyverse,ape, stringr, readr,
                 tidytree, ggtree, ggnewscale, treesliceR, paletteer, viridis) #installing and/or loading packages

# read phylogeny ---- 
tree_nitfix <- read.tree("3.trees/mimosoid_calibrated_clean_updated.tre")
# check species names
sp_nitfix <- tree_nitfix$tip.label# 4. reading morpho data =====
morpho <- read.csv("1.datasets/moniq/continuous_data_cleaned-20260410.csv")
#checking NA in each species 
morpho$na_percent <- rowMeans(is.na(morpho[,3:32])) * 100
#nao incluir acacia
morpho_selected <- morpho %>%
  filter(!str_detect(taxon, "Acacia"))
#2019 sp
#cortar a arvore pra 
tree_pruned <- drop.tip(tree_nitfix, setdiff(tree_nitfix$tip.label,
                                             morpho_selected$taxon))
#now we have 745


### time slice
SLICE <- 15
#para os principais grupos de mimoseae 
#get clades from time slice
time_slice <- treeSlice(tree_pruned, 
                        slice = max(nodeHeights(tree_pruned))-SLICE, 
                        trivial = F)
taxa_time_slice <- bind_rows(lapply(time_slice, function(x) data.frame(taxon = x$tip.label)), .id = 'clade') %>%
  mutate(clade = as.numeric(clade))

#filttering species from phylogeny 
morpho_phylo <- morpho_selected %>%
  filter(taxon %in% tree_nitfix$tip.label)
#now we have 850 species 
#81 subclados
## selecionar 10% de cada clado do time slice, considerando o que tiver menos NA 
set.seed(123) #permite replicabilidade

amostra_10pct <- taxa_time_slice %>%
  left_join(morpho_phylo %>% select(taxon, na_percent), by = "taxon") %>%
  group_by(clade) %>%
  arrange(na_percent, .by_group = TRUE) %>%
  mutate(n_clado = n(),
         n_selecionar = pmax(1, ceiling(n_clado * 0.10))) %>%
  filter(row_number() <= n_selecionar) %>%
  ungroup()

#conferir numero de especies por cada subclado 
amostra_10pct %>%
  group_by(clade) %>%
  summarise(n_total = unique(n_clado), n_amostrado = n()) %>%
  arrange(n_total)
#166 obs 
morpho10pct <- morpho_phylo %>%
  filter(taxon %in% amostra_10pct$taxon)

