### REDUÇÃO DE DIMENSIONALIDADE (GOWER) E
# ANALISES DE DISPARIDADE 
library(dplyr)
library(tibble)

#read phylogenetic tree and ecological data
traits_phylo <- read.csv()
biomes <- 
tree <- read.tree("4.trees/mimosoid_calibrated_clean_updated.tre")
## prune phylogeny
traits_phylo <- traits_phylo %>% rownames_to_column("species") 

tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, traits_phylo$species))

traits_phylo <- traits_phylo %>% column_to_rownames("species")

##Gower distance x PCoA =====------ 
library(cluster)
gower <- as.matrix(daisy(traits_phylo, metric = "gower")) ##cluster package
gower_modemean <- as.matrix(daisy(traits_modemean, metric = "gower"))

mantel(gower_modemean,gower) #as duas matrizes estao bem correlacionadas

## checagem: quantos pares tem NA na distancia (caso alguma linha nao compartilhe
## nenhuma variavel observada com outra - daria distancia NA)
#como fiz a imputação vai dar 0 
sum(is.na(as.matrix(gower_phylo)))

## PCoA com a matriz de gower =====================================
#como a matriz nao é euclidiana, 
# talvez a gente possa considerar uma correção para autovalores negativos
pcoa_res <- pcoa(gower, correction = "cailliez")  #ape package

## garantir ordem identica a arvore 
scores_pcoa <- pcoa_res$vectors
scores_pcoa <- scores_pcoa[match(tree_pruned$tip.label, rownames(scores_pcoa)), ]
stopifnot(identical(rownames(scores_pcoa), tree_pruned$tip.label))

#contribuicao de cada eixo
pcoa_values <- pcoa_res$values
#coordenadas por especie
pcoa_res$vectors

#visualizar em porcentagem a contribuicao dos eixos
percent_explained <- 100 * pcoa_values$Eigenvalues / 
  sum(pcoa_values$Eigenvalues[pcoa_values$Eigenvalues > 0])
percent_explained

#plot pcoa
pcoa_scores <- as.data.frame(pcoa_res$vectors)
pcoa_scores$species <- rownames(pcoa_scores)

library(ggplot2)
ggplot(pcoa_scores, aes(x = Axis.1, y = Axis.2)) +
  geom_point(size = 3) +
  geom_text(aes(label = species), vjust = -0.5) +
  theme_classic()


##incluindo a informaçao de clados
pcoa_df <- pcoa_res$vectors %>%
  as.data.frame() %>%
  rownames_to_column("species") %>%
  left_join(
    traits_final %>%
      select(species, clade),
    by = "species"
  )

##colorindo por clado
ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, color = clade)) +
  geom_point(size = 3) +
  theme_classic() +
  labs(
    x = "PCoA axis 1",
    y = "PCoA axis 2",
    color = "Clade"
  )

#definindo por generos
genus <- data.frame("genus" = sub("_.*" ,"", 
                                  row.names(pcoa_res$vectors)),
                    "species" = row.names(pcoa_res$vectors))

pcoa_data <- as.data.frame(pcoa_res$vectors.cor[, c(1, 2)])
pcoa_data$species <- row.names(pcoa_data)
pcoa_data <- merge(pcoa_data, genus, by = "species")

library(ggrepel)
ggplot(pcoa_data, aes(x = Axis.1, y = Axis.2, color = genus)) +
  geom_point(size = 3) +
  geom_text_repel(
    aes(label = species),
    size = 2.5,
    show.legend = FALSE
  ) +
  theme_classic()

ggplot(pcoa_res$vectors.cor[,c(1,2)], aes(x = Axis.1, y = Axis.2)) +
  geom_point(size = 3) +
  geom_text(aes(label = row.names(pcoa_res$vectors.cor)), 
           vjust = -0.5, size = 2) +
  theme_classic()

##coord fixed to adjust scale 
sov_clade <- pcoa_df %>%
  group_by(clade) %>%
  summarise(
    SOV_Axis1 = max(Axis.1, na.rm = TRUE) -
      min(Axis.1, na.rm = TRUE),
    
    SOV_Axis2 = max(Axis.2, na.rm = TRUE) -
      min(Axis.2, na.rm = TRUE),
    
    SOV_Axis3 = max(Axis.3, na.rm = TRUE) -
      min(Axis.3, na.rm = TRUE),
    
    SOV = SOV_Axis1 +
      SOV_Axis2 +
      SOV_Axis3
  )


##mean pairwise distance
mpd_clade <- pcoa_df %>%
  group_by(clade) %>%
  group_modify(~ {
    
    coords <- .x %>%
      select(Axis.1, Axis.2, Axis.3)
    
    d <- as.matrix(dist(coords))
    
    data.frame(
      MPD = mean(d[upper.tri(d)])
    )
  })

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
    traits_final %>%
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

##calcular metricas de disparidade =======
## SV, SR, MPD
## --- SV (sum of variances) e SR (sum of ranges) via dispRity,
##     calculadas sobre TODOS os eixos retidos do PCoA (nao so os 2 primeiros,
##     para nao perder disparidade que esta em eixos de ordem maior) ---
# objeto dispRity a partir dos eixos da PCoA


##using claddis
## create Claddis cladistic matrix
library(Claddis)
as.matrix(traits_phylo) -> traits

#precisaria substituir valores que sao nomes pra categorias
cladistic_matrix <- build_cladistic_matrix(traits, ordering = ord)

## calculating distance matrix
dist <- calculate_morphological_distances(cladistic_matrix)

## pcoa
pcoa <- ape::pcoa(dist$distance_matrix, correction = "cailliez")

## calculating variance explained explained by each principal component
eigen <- pcoa_res$values$Corr_eig

# PC1
round(eigen[1] / sum(eigen[eigen > 0]) * 100, 2)
#10,98
# PC2
round(eigen[2] / sum(eigen[eigen > 0]) * 100, 2)
#7,35
# PC3
round(eigen[3] / sum(eigen[eigen > 0]) * 100, 2)
#5,37

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

