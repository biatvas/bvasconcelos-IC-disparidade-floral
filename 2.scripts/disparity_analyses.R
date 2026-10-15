### REDUÇÃO DE DIMENSIONALIDADE (GOWER) E
# ANALISES DE DISPARIDADE 
library(dplyr)
library(tibble)
library(purrr)
library(ggrepel)
library(ggplot2)
library(cluster)

#read phylogenetic tree and ecological data
traits <- read.csv("3.outputs/mimoseae_flower_data.csv")
biomes <- read.csv("3.outputs/ecological_data_mimoseae.csv")
tree <- read.tree("4.trees/mimosoid_calibrated_clean_updated.tre")

## prune phylogeny
tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, traits$species))

#salvar info de clados
clade_info <- traits %>%
 select(species, clade)

traits <- traits %>%
  as.data.frame() %>%
  remove_rownames() %>%
  column_to_rownames("species")


continuous_cols <- c( grep("_mean$", colnames(traits), value = TRUE)) 
# Definindo quais são as colunas categóricas 
# categorical_cols <- setdiff(colnames(traits), continuous_cols)

# log nos traços contínuos
traits_log <- traits %>%
mutate(across(all_of(continuous_cols), log))

#transformar character em factor
categorical_cols <- c(
  "perianth_variable",
  "anther_gland_presence",
  "nectary",
  "filament_color",
  "filament_union"
)

traits_log[categorical_cols] <- lapply(
  traits_log[categorical_cols],
  factor
)

#3 é perianth variable e 5 é anther gland presence, 6 é nectary
traits <- traits_log %>%
  select(-clade)

#lidar com polimorfismos
#alguns traços foram corrigidos pra manter como 5&6 por exemplo, 
#quero sortear algum dos dois numeros pra entrar na analise 
# Identificar polimorfismos nas variáveis categóricas
polymorphism_summary <- tibble(
  variable = categorical_cols,
  n_polymorphic = sapply(categorical_cols, function(col) {
    x <- as.character(traits[[col]])
    sum(grepl("&", x, fixed = TRUE), na.rm = TRUE)
  })
) %>%
  filter(n_polymorphic > 0) %>%
  arrange(desc(n_polymorphic))

polymorphism_summary

#funçao para resolver polimorfismos
resolve_polymorphisms <- function(data, columns) {
  
  data_resolved <- data
  
  for (col in columns) {
    
    x <- as.character(data_resolved[[col]])
    
    idx <- which(
      grepl("&", x, fixed = TRUE) & !is.na(x)
    )
    
    for (i in idx) {
      
      states <- strsplit(
        x[i],
        "&",
        fixed = TRUE
      )[[1]]
      
      states <- trimws(states)
      states <- states[states != ""]
      
      if (length(states) > 0) {
        x[i] <- sample(states, size = 1)
      }
    }
    
    data_resolved[[col]] <- factor(
      x,
      levels = unique(c(
        levels(data[[col]]),
        x
      ))
    )
  }
  
  data_resolved
}

traits <- resolve_polymorphisms(data = traits, columns = "flower_merosity")
traits <- resolve_polymorphisms(data = traits, columns = "filament_color")

#resolve stamen count with NA (imputed as mean or according to phylogeny, thus, continuous values)
traits$stamen_count <- round(traits$stamen_count, digits = 0)
traits$stamen_count <- as.integer(traits$stamen_count)
traits$flower_merosity <- as.integer(traits$flower_merosity)

str(traits)
unique(traits$flower_merosity)

##Gower distance x PCoA =====------ 
gower <- as.matrix(daisy(traits, metric = "gower"))

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
ggplot(pcoa_scores, aes(x = Axis.2, y = Axis.3)) +
  geom_point(size = 3) +
  geom_text(aes(label = species), vjust = -0.5) +
  theme_classic()


##incluindo a informaçao de clados
traits_info <- traits %>% rownames_to_column("species") 
traits_info$clade <- clade_info$clade

pcoa_df <- pcoa$vectors %>%
  as.data.frame() %>%
  rownames_to_column("species") %>%
  left_join(
    traits_info %>%
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
  ) + ggplot2::coord_equal()

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
  theme_classic() + 
  ggplot2::coord_equal()

##plot morphospace
ggplot(pcoa$vectors.cor[,c(1,2)], aes(x = Axis.1, y = Axis.2)) +
  geom_point(size = 3) +
  geom_text(aes(label = row.names(pcoa$vectors.cor)), 
           vjust = -0.5, size = 2) +
  theme_classic()

##definindo por biomas 
#eco-morphospace
data_vectors_eco <- pcoa$vectors.cor[, c (1, 2)]

data_vectors_eco <- as.data.frame(data_vectors_eco) %>%
  add_column(as.factor(biomes$biome), .before = T) %>%
  colnames(data_vectors_eco) <- c("biome", "PCoA1", "PCoA2")

# traits$species <- gsub(pattern = "_"," ", traits$species)
# 
# setdiff(traits$species, biomes$species)
# 
# traits_2 <- traits[which(traits$species %in% biomes$species),] #153
# 
# all(traits_2$species %in% biomes$species)

##plot a phylomorphospace 
##phylomorphospace ======
# phylomorphospace exige matrix numerica pura, na mesma ordem da arvore
library(phytools)

space_axis <- pcoa$vectors.cor[, c(1,2)]
space_axis <- space_axis[rownames(space_axis) %in% tree$tip.label, ]

phylomorphospace <- phylomorphospace(
  tree_pruned,
  space_axis,
  xlab = "PCoA1",
  ylab = "PCoA2",
  label = "off"
) + ggplot2::coord_equal()


#PCA Hill-Smith =======
library(ade4)
library(adegraphics)

hs <- dudi.hillsmith(
  traits,
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
     pch = 19)

# + text(hs$li[,1],hs$li[,2],labels = traits$specie, pos = 1)

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
    traits_info %>%
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


# objeto dispRity a partir dos eixos da PCoA
bootstraps <- 1000
rarefaction <- "min" # "min" or FALSE

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

### Check code of Casali et al. (2023)
### available at https://doi.org/10.5281/zenodo.7240006
