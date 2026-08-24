#===============================================================#
# Head ----------------------------------------------------------------------
setwd("~/Documents/GitHub/bvasconcelos-IC-disparidade-floral/") #defining work directory
if (!require(librarian)) install.packages("librarian"); library("librarian")
librarian::shelf(phytools, dplyr, tidyr, purrr, factoextra, vegan, tidyverse,ape, stringr, readr,
                 tidytree, ggtree) #installing and/or loading packages

# read phylogeny ---- 
tree_nitfix <- read.tree("4.trees/mimosoid_calibrated_clean_updated.tre")
# check species names
sp_nitfix <- tree_nitfix$tip.label

#===============================================================#
# reading morpho data =====
morpho <- read.csv("3.outputs/continuous_data_cleaned-20260410.csv")
#checking NA in each species 
morpho$na_percent <- rowMeans(is.na(morpho[,3:31])) * 100

#excluir acacia pois o sampling esta sendo validado separadamente
morpho <- morpho %>%
  filter(!str_detect(taxon, "Acacia"))

#excluindo especies duplicadas porque apos o tnrs as subspecies foram sumarizadas em especie
morpho <- morpho %>%
  group_by(taxon) %>%
  slice_min(order_by = na_percent, n = 1, with_ties = FALSE) %>%
  ungroup()

sum(duplicated(morpho$taxon))

#Time Slice for taxon sampling 
# Setup =======================================================================
## libraries ------------------------------------------------------------------
library(viridis)
library(paletteer)
library(ggnewscale)
library(treesliceR)

#we have 1219 species on the phylogeny 
#morphodata contains 1836 species
#prune phylogeny with dataset 
tree_pruned <- drop.tip(tree_nitfix, setdiff(tree_nitfix$tip.label,
                                             morpho$taxon))  #745 species 
                        
### time slice

SLICE <- 15

#para os principais grupos de mimoseae 
#get clades from time slice
time_slice <- treeSlice(tree_pruned, 
                        slice = max(nodeHeights(tree_pruned))-SLICE, 
                        trivial = T)
taxa_time_slice <- bind_rows(lapply(time_slice, function(x) data.frame(taxon = x$tip.label)), .id = 'clade') %>%
  mutate(clade = as.numeric(clade))

##165 clados para time slice de 15 

#you can use the argument criterion to choose slice for 
# million years (my) or phylogenetic diversity (pd)

#===============================================================#
#filttering species from phylogeny 
morpho_phylo <- morpho %>%
  filter(taxon %in% tree_nitfix$tip.label)
#now we have 745 species 

all(tree_pruned$tip.label %in% morpho_phylo$taxon)

#===============================================================#
#SELECIONAR 10% DAS ESPÉCIES POR GÊNERO ======
mimoseae_generos <- read.csv("1.datasets/mimoseae_generos.csv") #tabela com numero de especies em cada genero de mimoseae
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
    ))

#creating a column with genus from each species 
morpho_filtered <- morpho_phylo %>%
  separate(taxon, into = c("genero", "epithet"),
           sep = "_", remove = FALSE)

morpho_filtered <- morpho_filtered %>% select(-"epithet")

#inserindo a porcentagem na tabela 
morpho_with_target <- morpho_filtered %>%
  left_join(
    mimoseae_generos %>% select(genero, n_target),
    by = "genero"
  )

# sorteio das espécies
set.seed(123) #permite replicabilidade

##testando ajustes no script
morpho_10perc <- morpho_with_target %>%
  group_by(genero) %>%
  group_modify(~ {
    
    n_needed <- unique(.x$n_target)   #n_target (10%) desse gênero específico
    n_available <- nrow(.x)           #quantas espécies desse gênero estão disponíveis
    
    n_select <- min(n_needed, n_available, na.rm = TRUE)
    slice_sample(.x, n = n_select)
  }) %>%
  ungroup()
#284 species selecionadas
length(unique(morpho_10perc$genero) %in% unique(morpho_filtered$genero))

# salvar nova tabela
write.csv(morpho_10perc,
          "3.outputs/morpho_10perc_all.csv",
          row.names = FALSE)

#===============================================================#
# selecionar espécies por corte temporal =====
data_slice <- left_join(morpho_phylo, taxa_time_slice, by = "taxon")

# morpho_filtered <- data_slice %>%
#  filter(na_percent <= 50) #usando apenas espécies com menos de 50% de NA 
#609 especies sobraram aqui 

sampled_species <- data_slice %>% 
  filter(!is.na(clade)) %>% #considerando apenas especies na filogenia 
  group_by(clade) %>%
  reframe(slice_sample(cur_data(),
                       n = max(1, round(0.2*n())))) #20% selecionado
##clados com menos de 5 sp nao entram na porcentagem 

#10% selecionou 191 de espécies
#20% 239 especies 

write.csv(sampled_species,
          "3.outputs/morpho_10perc_timeslice.csv",
          row.names = FALSE)

##===============================================================#
# visualizar como está a distribuição filogenetica das sp ====
librarian::shelf(patchwork)  # pra juntar os dois plots

# --- Plot 1: seleção por corte temporal (sampled_species, 20% por clado) ---
selected_taxa_temporal <- sampled_species$taxon

tip_data_temporal <- data.frame(label = tree_pruned$tip.label)
tip_data_temporal$sampled <- ifelse(tip_data_temporal$label %in% selected_taxa_temporal,
                                    "selected", "not_selected")

p_temporal <- ggtree(tree_pruned, size = 0.1, layout = "circular") %<+% tip_data_temporal +
  geom_tippoint(aes(color = sampled), size = 0.3) +
  scale_color_manual(values = c(selected = "red", not_selected = "grey80")) +
  labs(title = "Corte temporal (20% por clado)") +
  theme(legend.position = "right")

# --- Plot 2: seleção por gênero (morpho_10perc, 10% por gênero) ---
selected_taxa_genus <- morpho_10perc$taxon

tip_data_genus <- data.frame(label = tree_pruned$tip.label) #pruned é a arvore sem acacia
tip_data_genus$sampled <- ifelse(tip_data_genus$label %in% selected_taxa_genus,
                                 "selected", "not_selected")

p_genus <- ggtree(tree_pruned, size = 0.1, layout = "circular") %<+% tip_data_genus +
  geom_tippoint(aes(color = sampled), size = 0.3) +
  scale_color_manual(values = c(selected = "blue", not_selected = "grey80")) +
  labs(title = "Por gênero (10%)") +
  theme(legend.position = "right")

# --- Combinar lado a lado ---
p_temporal + p_genus + plot_layout(ncol = 2, guides = "collect")

# opcional: salvar em alta resolução
ggsave("5.figuras/comparacao_amostragem.png",
       width = 14, height = 7, dpi = 300)
