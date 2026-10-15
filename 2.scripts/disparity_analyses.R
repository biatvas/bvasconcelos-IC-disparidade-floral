### REDUÇÃO DE DIMENSIONALIDADE (GOWER) E
# ANALISES DE DISPARIDADE 
library(dplyr)
library(tibble)
library(purrr)
library(ggrepel)

#read phylogenetic tree and ecological data
traits <- read.csv("3.outputs/mimoseae_flower_data.csv")
biomes <- read.csv("3.outputs/ecological_data_mimoseae.csv")
tree <- read.tree("4.trees/mimosoid_calibrated_clean_updated.tre")

#prune morphodata to match biomes dataset
traits$species <- gsub(pattern = "_"," ", traits$species)

traits <- traits[traits$species %in% biomes$species, ]
#all(traits$species %in% biomes$species)

traits$species <- gsub(pattern = " ","_", traits$species)
## prune phylogeny
tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, traits$species))

#clade_info <- traits %>%
#  select(species, clade)
traits <- traits %>%
  as.data.frame() %>%
  remove_rownames() %>%
  column_to_rownames("species")


continuous_cols <- c( grep("_mean$", colnames(traits), value = TRUE), "stamen_count") 
# Definindo quais são as colunas categóricas 
categorical_cols <- setdiff(colnames(traits), continuous_cols)

# log nos traços contínuos
traits_log <- traits %>%
mutate(across(all_of(continuous_cols), log))

#transformar character em factor
traits_log[] <- lapply(traits_log, function(x) {
  if (is.character(x)) as.factor(x) else x
})

#binarias nominais vao virar factor tambem
traits_log[, c(2, 5)] <- lapply(
  traits_log[, c(2, 5), drop = FALSE],
  factor
)

##Gower distance x PCoA =====------ 
gower <- daisy(traits_log, metric = "gower")
matrix_phylo <- as.matrix(gower)

## checagem: quantos pares tem NA na distancia (caso alguma linha nao compartilhe
## nenhuma variavel observada com outra - daria distancia NA)
#como fiz a imputação vai dar 0 
sum(is.na(as.matrix(gower)))

## PCoA com a matriz de gower =====================================
#como a matriz nao é euclidiana, 
# talvez a gente possa considerar uma correção para autovalores negativos
pcoa <- pcoa(gower, correction = "cailliez")  #ape package

## calculating variance explained explained by each principal component
eigen <- pcoa$values$Corr_eig

# PC1
round(eigen[1] / sum(eigen[eigen > 0]) * 100, 2)
# PC2
round(eigen[2] / sum(eigen[eigen > 0]) * 100, 2)
# PC3
round(eigen[3] / sum(eigen[eigen > 0]) * 100, 2)


## garantir ordem identica a arvore 
scores_pcoa <- pcoa$vectors
scores_pcoa <- scores_pcoa[match(tree_pruned$tip.label, rownames(scores_pcoa)), ]
stopifnot(identical(rownames(scores_pcoa), tree_pruned$tip.label))

#plot pcoa
pcoa_scores <- as.data.frame(pcoa$vectors)
pcoa_scores$species <- rownames(pcoa_scores)

#plotar o morfoespaço
library(ggplot2)
ggplot(pcoa_scores, aes(x = Axis.1, y = Axis.2)) +
  geom_point(size = 3) +
  geom_text(aes(label = species), vjust = -0.5) +
  theme_classic()


##incluindo a informaçao de clados
# traits_log <- traits %>% rownames_to_column("species")
# 
# pcoa_df <- pcoa_res$vectors %>%
#   as.data.frame() %>% 
#   rownames_to_column("species") %>%
#   left_join(
#     traits %>%
#       select(species, clade),
#     by = "species"
#   )

##colorindo por clado
# ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, color = clade)) +
#   geom_point(size = 3) +
#   theme_classic() +
#   labs(
#     x = "PCoA axis 1",
#     y = "PCoA axis 2",
#     color = "Clade"
#   )

#definindo por generos
genus <- data.frame("genus" = sub("_.*" ,"", 
                                  row.names(pcoa$vectors)),
                    "species" = row.names(pcoa$vectors))

pcoa_data <- as.data.frame(pcoa$vectors.cor[, c(1, 2)])
pcoa_data$species <- row.names(pcoa_data)
pcoa_data <- merge(pcoa_data, genus, by = "species")


#plotar morfoespaço colorido por genero
ggplot(pcoa_data, aes(x = Axis.1, y = Axis.2, color = genus)) +
  geom_point(size = 3) +
  geom_text(
    aes(label = species),
    size = 2.5,
    show.legend = FALSE
  ) +
  theme_classic() + 
  ggplot2::coord_equal()

ggplot(pcoa$vectors.cor[,c(1,2)], aes(x = Axis.1, y = Axis.2)) +
  geom_point(size = 3) +
  geom_text(aes(label = row.names(pcoa_res$vectors.cor)), 
           vjust = -0.5, size = 2) +
  theme_classic()

##plot a phylomorphospace 
space_axis <- pcoa$vectors.cor[, c(1,2)]

space_axis <- space_axis[rownames(space_axis) %in% tree$tip.label, ]

phylomorphospace <- phylomorphospace(
  tree_pruned,
  space_axis,
  xlab = "PC1",
  ylab = "PC2",
  label = "off"
) + ggplot2::coord_equal()

#eco-morphospace
data_vectors_eco <- pcoa$vectors.cor[, c (1, 2)]

data_vectors_eco <- as.data.frame(data_vectors_eco) %>%
  add_column(as.factor(biomes$biome), .before = T) %>%
colnames(data_vectors_eco) <- c("biome", "PC1", "PC2")

# checking species of each vegetation type #

groups_veg_type <- split(dat_eco$species, dat_eco$biome)




#####Disparity metrics
# objeto dispRity a partir dos eixos da PCoA
library(dispRity)
gower_mat <- as.matrix(gower)
disp_obj_dist <- custom.subsets(gower_mat, group = list(all = rownames(gower_mat)))

## Bootstrapping the data
bootstrapped_data <- boot.matrix(disp_obj_dist, bootstraps = 100)
## Calculating the sum of variances
sum_of_variances <- dispRity(bootstrapped_data, metric = c(sum, variances))
summary(sum_of_variances)

mpd <- dispRity(disp_obj_dist, metric = c(mean, pairwise.dist))
sov <- dispRity(disp_obj_dist, metric = c(sum, variances))
sor <- dispRity(disp_obj_dist, metric = c(sum, ranges))


##calcular metricas de disparidade =======
## SV, SR, MPD
## --- SV (sum of variances) e SR (sum of ranges) via dispRity,
##     calculadas sobre TODOS os eixos retidos do PCoA (nao so os 2 primeiros,
##     para nao perder disparidade que esta em eixos de ordem maior) ---
# objeto dispRity a partir dos eixos da PCoA



##plot clades ======


#plot ecomorphospace ======


### sinal filogenetico?
# δ (delta) statistic Borges

##phylomorphospace ======
# phylomorphospace exige matrix numerica pura, na mesma ordem da arvore
library(phytools)
scores_mat <- as.matrix(pcoa_scores)

phylomorphospace(
  tree_pruned,
  scores_mat[,1:2],
  xlab = "PC1",
  ylab = "PC2",
  label = "off"
)


## Size - SV ------------------------------------------------------------------
# sum of variances - post-ordination metric

### This part of the script was based on the code of Casali et al. (2023)
### available at https://doi.org/10.5281/zenodo.7240006

bootstraps <- 1000
rarefaction <- "min" # "min" or FALSE


#PCA Hill-Smith =======
library(ade4)
library(adegraphics)

hs <- dudi.hillsmith(
  traits_log,
  scannf = FALSE,
  nf = 5
)

#select 2 
#screeplot(hs)
summary(hs) #23 eicxos (cada categoria de var quant é um eixo)
# até o eixo 5, acumula 43.29% da explicacao
hs$eig #variancia generalizada explicada por eixo expalhada

axes_contribution <- 100*hs$eig/sum(hs$eig)

#plot flowers in morpho space
plot(hs$li[,1],
     hs$li[,2],
     xlab = "pc1",
     ylab = "pc2",
     pch = 19) +
  
  text(hs$li[,1],hs$li[,2],labels = traits_log$specie, pos = 1)

hs$cr # quanto cada variavel esta associada a cada eixo
hs$index # tipos de cada var

scatter(hs)
s.label(hs$li, labels = NULL) 
#essa estrutura mais achatada pode tar
# refletindo a baixa explicacao por eixo

hs$cr[1] #importancia de cada variavel para cada eixo. 
# tamanho floral com maior peso. podemos padronizar as var. continuas
# pela media geometrica 

hs$c1 #laodings (autovetores) 
hs$co[1] #c1 reescalonado pelos autovalores.
# direcao (quantitativas) ou qais categorias puxam pra qual lado
# do eixo
hs_df <- hs$li %>%
  as.data.frame() %>%
  rownames_to_column("species") %>%
  left_join(
    traits_log %>%
      select(species, clade),
    by = "species"
  )

ggplot(hs_df, aes(x = Axis1, y = Axis2, color = clade)) +
  geom_point(size = 3) +
  theme_classic() +
  labs(
    x = "PC1",
    y = "PC2",
    color = "Clade"
  )
