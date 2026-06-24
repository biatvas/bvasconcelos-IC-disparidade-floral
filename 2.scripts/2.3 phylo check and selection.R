#===============================================================#
# Head ----------------------------------------------------------------------
setwd("~/Documents/Labis/Dados/") #defining work directory
if (!require(librarian)) install.packages("librarian"); library("librarian")
librarian::shelf(phytools, dplyr, tidyr, purrr, factoextra, vegan, tidyverse,ape, stringr, readr,
                 tidytree, ggtree) #installing and/or loading packages

# read phylogeny ---- 
tree_nitfix <- read.tree("3.trees/mimosoid_calibrated_clean_updated.tre")

# check species names
sp_nitfix <- tree_nitfix$tip.label

#===============================================================#
# 3. correction of names 
# Função para padronizar nomes
#clean_names <- function(x){
#  x %>%
#   str_replace_all("_", " ") %>%
#   str_trim()}

# sp_nitfix <- clean_names(sp_nitfix)

#===============================================================#
# 4. reading morpho data =====

morpho <- read.csv("1.datasets/continuous_data_cleaned-20260410.csv")

#checking NA in each species 
morpho$na_percent <- rowMeans(is.na(morpho[,3:32])) * 100

#Time Slice for taxon sampling 
# Setup =======================================================================

## libraries ------------------------------------------------------------------
library(phytools)
library(tidyverse)
library(ggtree)
library(viridis)
library(paletteer)
library(ggnewscale)
library(treesliceR)
#we have 1219 species
#morphodata contains 3172 species
#prune phylogeny with dataset 
tree_pruned <- drop.tip(tree_nitfix, setdiff(tree_nitfix$tip.label,
                                             morpho$taxon))
#now we have 1166
### time slice
SLICE <- 15
#para os principais grupos de mimoseae 
#get clades from time slice
time_slice <- treeSlice(tree_pruned, 
                        slice = max(nodeHeights(tree_pruned))-SLICE, 
                        trivial = F)
taxa_time_slice <- bind_rows(lapply(time_slice, function(x) data.frame(taxon = x$tip.label)), .id = 'clade') %>%
  mutate(clade = as.numeric(clade))

#you can use the argument criterion to choose slice for 
# million years (my) or phylogenetic diversity (pd)

#===============================================================#
# 5. checking species from phylogeny =======
#nitfix_genus_count <- sp_nitfix %>% 
#  tibble(species = .) %>% 
#  separate(species, into = c("genus", "epithet"),
#           sep = " ") %>% 
#  distinct() %>% 
#  count(genus, name = "n_nitfix_species")

#morpho$phylo <- morpho$species %in% sp_nitfix #check which species are in phylo 
#sum(morpho$phylo, na.rm = TRUE) #1283 sp, sendo que temos 1219 sp no nitfix
#conferir espécies que estão com nomes repetidos
#morpho_with_phylo <- morpho %>% filter(phylo == TRUE)
#write.csv(morpho_with_phylo, "4.outputs/species_features_phylo_clean.csv", row.names = FALSE)
#vou excluir os tips duplicados dessa tabela com apenas sp da phylo 

#especies_duplicadas <- morpho_with_phylo$species[duplicated(morpho_with_phylo$species)]
#118 especies duplicadas 
#write.table(especies_duplicadas, file = "4.outputs/duplicate_data_tips.txt", row.names = FALSE)
#clean duplicate tips from dataset 
# open cleaned dataset after manual check 
#morpho <- read.csv("4.outputs/species_features_phylo_clean_check.csv")
#agora o dataset contêm 1168 espécies 

#filttering species from phylogeny 
morpho_phylo <- morpho %>%
  filter(taxon %in% tree_nitfix$tip.label)

#===============================================================#
# 6. SELECIONAR 10% DAS ESPÉCIES POR GÊNERO ======
mimoseae_generos <- read.csv("4.outputs/mimoseae_generos.csv")#tabela com numero de especies em cada genero de mimoseae
mimoseae_generos <- mimoseae_generos %>% separate(Genero, into = c("genero", "author"), 
                                    sep = " ")
#calculating the percentage of each genus 
mimoseae_generos <- mimoseae_generos %>%
  mutate(
    spp = as.numeric(spp),
    n_target = if_else( #generos com menos de 5 sp não entram na porcentagem 
      spp > 5,
      ceiling(spp * 0.10),
      spp
    )
  )

#ao inves de fazer isso posso trabalhar com o dataset ja filtrado 

#creating a column with genus from each species 
morpho_filtered <- morpho_filtered %>%
  separate(species, into = c("genus", "epithet"),
           sep = " ", remove = FALSE)

morpho_filtered <- morpho_filtered %>% select(-"epithet")

#creating a column with genus from each species 
morpho_phylo <- morpho_phylo %>%
  separate(species, into = c("genus", "epithet"),
           sep = " ", remove = FALSE)

#corrigindo o nome da coluna de generos pra extrair no dataset
mimoseae_generos <- mimoseae_generos %>%
  rename(genus = genero)
#inserindo a porcentagem na tabela 
morpho_with_target <- morpho_phylo %>%
  left_join(
    mimoseae_generos %>% select(genus, n_target),
    by = "genus"
  )

#excluindo generos que não estão na filogenia 
morpho_with_target <- morpho_with_target %>%
  filter(!is.na(n_target))

# sorteio das espécies
set.seed(123) #permite replicabilidade

morpho_10perc <- morpho_with_target %>%
  group_by(genus) %>%
  group_modify(~ {
    
    n_needed <- unique(.x$n_target)   # 10% calculado pelo total do gênero
    n_available <- nrow(.x)           # quantas espécies desse gênero estão na filogenia
    
    n_select <- min(n_needed, n_available)
    
    slice_sample(.x, n = n_select)
  }) %>%
  ungroup()

# salvar nova tabela
write.csv(morpho_10perc,
          "4.outputs/morpho_10perc_all.csv",
          row.names = FALSE)
#===============================================================#
# 7. SELECIONAR ESPÉCIES CORTE TEMPORAL =====
data_slice <- left_join(morpho_phylo, phylo_taxa, by = "taxon")

morpho_filtered <- data_slice %>%
  filter(na_percent <= 50) #usando apenas espécies com menos de 50% de NA 
#669 especies sobraram aqui 

set.seed(123)
sampled_species <- morpho_filtered %>% 
  filter(!is.na(clade)) %>% #considerando apenas especies na filogenia 
  group_by(clade) %>%
  reframe(slice_sample(cur_data(),
                       n = max(1, round(0.2*n()))))

##===============================================================#
# visualizar como está a distribuição filogenetica das sp ====
sample <- sampled_species$taxon

tip_data <- data.frame(
  label = tree_nitfix$tip.label
)

tip_data$sampled <- ifelse(tip_data$label %in% sample,
                           "selected",
                           "not_selected")

p <- ggtree(tree_nitfix, size = 0.1, layout = "circular") %<+% tip_data +
  geom_tippoint(aes(color = sampled), size = 0.3) +
  scale_color_manual(values = c(
    selected = "red",
    not_selected = "grey80"
  )) +
  theme(
    legend.position = "right"
  )

