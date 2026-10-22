### REDUÇÃO DE DIMENSIONALIDADE (GOWER) E
# ANALISES DE DISPARIDADE 
library(dplyr)
library(tibble)
library(purrr)
library(ggrepel)
library(ggplot2)
library(cluster)
library(ape)
library(vegan)

setwd("~/Documents/GitHub/bvasconcelos-IC-disparidade-floral/")

#read phylogenetic tree and ecological data
traits <- read.csv("3.outputs/mimoseae_flower_data.csv")
biomes <- read.csv("3.outputs/ecological_data_mimoseae.csv", stringsAsFactors = FALSE)
tree <- read.tree("4.trees/mimosoid_calibrated_clean_updated.tre")

## prune phylogeny
tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, traits$species))

#salvar info de clados
# Guardar metadados ANTES de remover clade
meta <- traits %>%
  select(species, clade) %>%
  distinct(species, .keep_all = TRUE)

# Depois, preparar os traits
traits <- traits %>%
  select(-clade) %>%
  as.data.frame() %>%
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
traits <- traits_log 
#%>%
#   select(-clade)

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
#traits$stamen_count <- as.numeric(as.character((traits$stamen_count)))
traits$flower_merosity <- as.numeric((as.character(traits$flower_merosity)))

str(traits)
unique(traits$flower_merosity)

##Gower distance x PCoA =====------ 
gower <- as.matrix(daisy(traits, metric = "gower"))

rownames(gower) <- gsub(pattern = " ","_", rownames(gower))

saveRDS(
  gower,
  "3.outputs/flower_gower_distance.rds")

## checagem: quantos pares tem NA na distancia (caso alguma linha nao compartilhe
## nenhuma variavel observada com outra - daria distancia NA)
#como fiz a imputação vai dar 0 
#sum(is.na(as.matrix(gower)))

###test data with NMDS
set.seed(123)

# gower_dist <- as.dist(gower)

nmds <- metaMDS(gower,
                k=2,
                trymax=100)
#check ordenation
nmds$stress
stressplot(nmds)

#it plots the observed dissimilarities against distances in ordination space,
#and also shows the non-metric monotone regression.
#distances in ordination space???????

#scores
scores_nmds <- as.data.frame(scores(nmds, display = "sites"))%>%
  rownames_to_column("species")

nmds_df <- scores_nmds %>%
  left_join(meta, by = "species")

table(nmds_df$clade, useNA = "ifany")

#definir cores para os clados
clade_levels <- unique(nmds_df$clade)

pal <- setNames(
  hcl.colors(length(clade_levels), palette = "Dark 3"),
  clade_levels
)

##plot morphospace
ggplot(
  nmds_df,
  aes(x = NMDS1, y = NMDS2, color = clade)
) +
  geom_point(size = 2.5, alpha = 0.85) +
  theme_classic() +
  coord_equal() +
  scale_color_manual(values = pal) +
  labs(
    x = "NMDS axis 1",
    y = "NMDS axis 2",
    color = "Clade",
    title = "Floral morphospace of Mimosoideae"
  )


## NMDS com polígonos convexos por clado

# Reordenar os metadados para corresponder aos scores
meta_nmds <- meta[
  match(rownames(nmds$points), meta$species),
]

stopifnot(
  identical(
    as.character(meta_nmds$species),
    rownames(nmds$points)
  )
)

# Garantir que os grupos estejam definidos na mesma ordem
groups <- factor(meta_nmds$clade)

# Paleta de cores
clade_levels <- levels(groups)
pal <- setNames(
  hcl.colors(length(clade_levels), palette = "Dark 3"),
  clade_levels
)

point_cols <- unname(pal[as.character(groups)])

# Plot
plot(
  nmds,
  type = "n",
  xlab = "NMDS axis 1",
  ylab = "NMDS axis 2",
  main = "Floral morphospace by clade"
)

# Polígonos convexos
ordihull(
  nmds,
  groups = groups,
  draw = "polygon",
  col = adjustcolor(pal, alpha.f = 0.20),
  border = NA
)

# Pontos das espécies
points(
  nmds,
  display = "sites",
  pch = 19,
  col = point_cols
)

legend(
  "topright",
  legend = clade_levels,
  col = pal,
  pch = 19,
  bty = "n",
  cex = 0.8
)

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
traits_info$clade <- meta$clade

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


##definindo grupos 
tree_pruned$tip.label <- gsub(pattern = "_", " ", tree_pruned$tip.label)
data_vectors_clade <- as.data.frame(pcoa$vectors.cor[, c (1, 2)])

data_vectors_clade %>% 
  add_column((meta$clade), (meta$species))-> data_vectors_clade

colnames(data_vectors_clade) <- c("PCoA1", "PCoA2", "Clade", "species")

data_vectors_clade <- data_vectors_clade %>% column_to_rownames("species")

data_vectors_clade <-
  data_vectors_clade[row.names(data_vectors_clade) %in% tree_pruned$tip.label, ]

##plot a phylomorphospace 
##phylomorphospace ======
# phylomorphospace exige matrix numerica pura, na mesma ordem da arvore

tree_pruned$tip.label <- gsub(pattern = " ", "_", tree_pruned$tip.label)
space_axis <- pcoa$vectors.cor[, c(1,2)]
space_axis <- space_axis[rownames(space_axis) %in% tree_pruned$tip.label, ]

phylomorphospace <- phylomorphospace(
  tree_pruned,
  space_axis,
  xlab = "PCoA1",
  ylab = "PCoA2",
  label = "off")

##extract config
phylo_data <- data.frame(
  xstart = phylomorphospace$xx[phylomorphospace$edge[, 1]],
  ystart = phylomorphospace$yy[phylomorphospace$edge[, 1]],
  xstop = phylomorphospace$xx[phylomorphospace$edge[, 2]],
  ystop = phylomorphospace$yy[phylomorphospace$edge[, 2]],
  nodestart = phylomorphospace$edge[, 1],
  nodestop = phylomorphospace$edge[, 2]
)

# Nomes que existem nos dados
clados_dados <- sort(unique(as.character(data_vectors_clade$Clade)))

# Nomes definidos na escala manual
clados_escala <- c(
  "Adenanthera clade",
  "Albizia clade",
  "Archidendron clade",
  "Inga clade",
  "Samanea clade",
  "Jupunba clade",
  "Mimosa clade",
  "Pithecellobium clade",
  "Cojoba clade",
  "Calliandra clade",
  "Zapoteca clade",
  "Stryphnodendron clade",
  "Senegalia grade",
  "Parkia clade",
  "Dichrostachys clade",
  "Neltuma clade",
  "Prosopis clade",
  "Newtonia grade",
  "Entada clade"
)

# Clados dos dados que não estão na escala
setdiff(clados_dados, clados_escala)

# Nomes da escala que não existem nos dados
setdiff(clados_escala, clados_dados)

# plot
phylomorphospace <-
  ggplot() +
  geom_segment(
    data = phylo_data,
    aes(
      x = xstart,
      y = ystart,
      xend = xstop,
      yend = ystop
    ),
    linewidth = 0.25,
    color = "grey75"
  ) +
  geom_point(
    data = data_vectors_clade,
    aes(
      x = PCoA1,
      y = PCoA2,
      color = Clade,
      shape = Clade
    ),
    size = 1.8,
    alpha = 0.8
  ) +
  scale_color_manual(
    values = paletteer_d("ggthemes::calc")[
      seq_along(unique(data_vectors_clade$Clade))
    ],
    name = "Clade"
  ) +
  scale_shape_manual(
    values = c(
      "Adenanthera clade" = 17,
      "Albizia clade" = 0,
      "Archidendron clade" = 4,
      "Inga clade" = 2,
      "Samanea clade" = 1,
      "Jupunba clade" = 14,
      "Mimosa clade" = 20,
      "Pithecellobium clade" = 5,
      "Cojoba clade" = 6,
      "Calliandra clade" = 15,
      "Zapoteca clade" = 13,
      "Stryphnodendron clade" = 9,
      "Senegalia grade" = 4,
      "Parkia clade" = 2,
      "Dichrostachys clade" = 8,
      "Neltuma clade" = 10,
      "Prosopis clade" = 19,
      "Newtonia grade" = 3,
      "Entada clade" = 7,
      "outgroup" = 16
    ),
    guide = "none"
  ) +
  coord_fixed(ratio = 1) +
  labs(
    x = paste0(
      "PCoA1 (",
      round(pcoa$values$Rel_corr_eig[1] * 100, 2),
      "%)"
    ),
    y = paste0(
      "PCoA2 (",
      round(pcoa$values$Rel_corr_eig[2] * 100, 2),
      "%)"
    )
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9),
    legend.key.height = grid::unit(0.35, "cm"),
    axis.text = element_text(color = "black")
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        shape = 16,
        size = 3,
        alpha = 1
      )
    )
  )

#check
phylomorphospace

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

# traits$species <- gsub(pattern = "_"," ", traits$species)
# 
# setdiff(traits$species, biomes$species)
# 
# traits_2 <- traits[which(traits$species %in% biomes$species),] #153
# 
# all(traits_2$species %in% biomes$species)


#PCA Hill-Smith =======
library(ade4)
library(adegraphics)

# quant_vars <- c(
#   "flower_merosity",
#   "stamen_count",
#   "calyx_length_mean",
#   "corolla_length_mean",
#   "corolla_lobe_length_mean",
#   "filament_length_mean"
# )
# 
# cat_vars <- c(
#   "perianth_variable",
#   "anther_gland_presence",
#   "nectary",
#   "filament_color",
#   "filament_union"
# )

###
hs <- dudi.hillsmith(
  traits,
  scannf = FALSE,
  nf = 5
)

#screeplot(hs, main = "variancia por eixo")

summary(hs) #23 eicxos (cada categoria de var quant é um eixo)
# até o eixo 5, acumula 43.29% da explicacao

#autovalores
hs$eig #variancia generalizada explicada por eixo expalhada

axes_contribution <- 100*hs$eig/sum(hs$eig)

#cumulative percentual
axes_cumulative <- cumsum(axes_contribution)

data.frame(
  eixo = seq_along(hs$eig),
  variancia = axes_contribution,
  acumulada = axes_cumulative
)
## 56% de variacao acumulada em 4 eixos
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

#plot flowers in morpho space
plot(hs$li[,1],
     hs$li[,2],
     xlab = "pc1",
     ylab = "pc2",
     pch = 19)

# + text(hs$li[,1],hs$li[,2],labels = traits$specie, pos = 1)

hs_df <- hs$li %>%
  as.data.frame() %>%
  rownames_to_column("species") %>%
  left_join(
    traits_info %>% select(species, clade),
    by = "species"
  )

# Conferir se todos os pontos receberam um clado
table(hs_df$clade, useNA = "ifany")

# direcao (quantitativas) ou qais categorias puxam pra qual lado
# do eixo
# hs_df <- hs$li %>%
#   as.data.frame() %>%
#   rownames_to_column("species") %>%
#   left_join(
#     traits_info %>%
#       select(species, clade),
#     by = "species"
#   )

ggplot(hs_df, aes(x = Axis1, y = Axis2, color = clade)) +
  geom_point(size = 2.5, alpha = 0.85) +
  coord_equal() +
  theme_classic() +
  labs(
    x = sprintf("Axis 1 (%.1f%%)", axes_contribution[1]),
    y = sprintf("Axis  2 (%.1f%%)", axes_contribution[2]),
    color = "Clado",
    title = "Floral Morphospace of tribe Mimoseae",
    subtitle = "Hill-Smith ordenation"
  )

##calcular metricas de disparidade =======
## SV, SR, MPD
## --- SV (sum of variances) e SR (sum of ranges) via dispRity,
# objeto dispRity a partir da matriz de distancia
library(dispRity)
traits_info$species <- gsub(pattern = "_"," ", (traits_info$species))

traits_2 <- traits[traits_info$species %in% biomes$species, ]

# mpd <- dispRity(disp_obj_dist, metric = c(mean, pairwise.dist))
# sov <- dispRity(disp_obj_dist, metric = c(sum, variances))
# sor <- dispRity(disp_obj_dist, metric = c(sum, ranges))

### Check code of Casali et al. (2023)
### available at https://doi.org/10.5281/zenodo.7240006

##check code of Clark et al. (2023)
#sum of variances for disparity of a group 
#Sum of range, extensao ocupada pelos grupos
#MPD distancia media entre os pares do grupo 

# Lista de espécies por clado
groups_clade <- split(
  meta$species,
  meta$clade
)

clade_disparity <- custom.subsets(
  gower,
  group = groups_clade)

## Bootstrapping the data
bootstrapped_data <- boot.matrix(clade_disparity, bootstraps = 100)

# Métricas de disparidade
sv_clade <- dispRity(
    bootstrapped_data,
    metric = c(sum, variances))

sr_clade <- dispRity(
    bootstrapped_data,
    metric = c(sum, ranges))

  
mpd_clade <- dispRity(
    bootstrapped_data,
    metric = c(mean, pairwise.dist))

##definir grupos de biomas 

#biomes names
#1 Tropical & Subtropical Moist Broadleaf Forests
#2 Tropical & Subtropical Dry Broadleaf Forests
#3 Tropical & Subtropical Coniferous Forests
#4 Temperate Broadleaf & Mixed Forests
#5 Temperate Conifer Forests
#7 Tropical & Subtropical Grasslands, Savannas & Shrublands
#8 Temperate Grasslands, Savannas & Shrublands
#9 Flooded Grasslands & Savannas
#10 Montane Grasslands & Shrublands
#12 Mediterranean Forests, Woodlands & Scrub
#13 Deserts & Xeric Shrublands
#14 Mangroves

# Conferir estrutura
names(biomes)
table(biomes$biome, useNA = "ifany")


# Remover espécies sem bioma
groups_biomes <- list(
  "Tropical & Subtropical Moist Broadleaf Forests" = filter(biomes, grepl("1", biome)) %>%
                            pull(species) %>% 
                            as.vector(),
  "Tropical & Subtropical Dry Broadleaf Forests" = filter(biomes, grepl("2", biome)) %>%
    pull(species) %>% 
    as.vector(),
  "Tropical & Subtropical Coniferous Forests" = filter(biomes, grepl("3", biome)) %>%
    pull(species) %>% 
    as.vector(),
  "Temperate Broadleaf & Mixed Forests" = filter(biomes, grepl("4", biome)) %>%
    pull(species) %>% 
    as.vector(),
  "Temperate Conifer Forests" = filter(biomes, grepl("5", biome)) %>%
    pull(species) %>% 
    as.vector(),
  "Tropical & Subtropical Grasslands, Savannas & Shrublands" = filter(biomes, grepl("7", biome)) %>%
    pull(species) %>% 
    as.vector(),
  "Temperate Grasslands, Savannas & Shrublands" = filter(biomes, grepl("8", biome)) %>%
    pull(species) %>% 
    as.vector(),
  "Flooded Grasslands & Savannas" = filter(biomes, grepl("9", biome)) %>%
    pull(species) %>% 
    as.vector(),
  "Montane Grasslands & Shrublands" = filter(biomes, grepl("10", biome)) %>%
    pull(species) %>% 
    as.vector(),
  "Mediterranean Forests, Woodlands & Scrub" = filter(biomes, grepl("12", biome)) %>%
    pull(species) %>% 
    as.vector(),
  "Deserts & Xeric Shrublands" = filter(biomes, grepl("13", biome)) %>%
    pull(species) %>% 
    as.vector(),
  "Mangroves" = filter(biomes, grepl("14", biome)) %>%
    pull(species) %>% 
    as.vector())

##preciso manter apenas especies que tem biomas 
# species_biomes <- biomes %>%
#   filter(!is.na(biome), biome != "") %>%
#   pull(species) %>%
#   unique()
# 
# # Manter na matriz somente espécies com bioma
# pcoa_biomes <- pcoa$vectors.cor[
#   rownames(pcoa$vectors.cor) %in% species_biomes,
# ]

### SV - Ecological groups
sub_ecological <- custom.subsets(
  pcoa_biomes,
  groups_biomes
)

bootstrapped_data <- boot.matrix(
  sub_ecological,
  bootstraps = bootstraps,
  rarefaction = rarefaction
)

ecological_sv <- dispRity(
  bootstrapped_data,
  metric = c(sum, variances)
)

ecological_sv_smry <- summary(ecological_sv)

ecological_sv_wilcox <- test.dispRity(
  ecological_sv,
  wilcox.test,
  correction = "holm"
)

ecological_sv_bhatt <- test.dispRity(
  ecological_sv,
  bhatt.coeff
)

### SV - Ecological groups
sub_ecological <- custom.subsets(pcoa_biomes, groups_biomes)
bootstrapped_data <- boot.matrix(sub_ecological, 
                                 bootstraps = bootstraps,
                                 rarefaction = rarefaction)
ecological_sv <- dispRity(bootstrapped_data, metric = c(sum, variances))
ecological_sv_smry <- summary(ecological_sv)
ecological_sv_wilcox <- test.dispRity(ecological_sv, 
                                      wilcox.test, 
                                      correction = "holm")
ecological_sv_bhatt <- test.dispRity(ecological_sv, bhatt.coeff)

### SV - Ecological lineages
sub_ecolineages <- custom.subsets(pcoa$vectors.cor, ecological_lineages)
bootstrapped_data <- boot.matrix(sub_ecolineages, 
                                 bootstraps = bootstraps,
                                 rarefaction = rarefaction)
ecolineages_sv <- dispRity(bootstrapped_data, metric = c(sum, variances))
ecolineages_sv_smry <- summary(ecolineages_sv)
ecolineages_sv_wilcox <- test.dispRity(ecolineages_sv,
                                       wilcox.test, 
                                       correction = "holm")
ecolineages_sv_bhatt <- test.dispRity(ecolineages_sv, bhatt.coeff)

##SALVAR RESULTADOS FINAIS
saveRDS(
  disparity_clade,
  "3.outputs/disparity_flower_clades.rds"
)

saveRDS(
  disparity_biome,
  "3.outputs/disparity_flower_biomes.rds"
)

# Exportar os resumos para consulta
capture.output(
  disparity_clade$SV_summary,
  file = "3.outputs/SV_clades.txt"
)

capture.output(
  disparity_clade$SR_summary,
  file = "3.outputs/SR_clades.txt"
)

capture.output(
  disparity_clade$MPD_summary,
  file = "3.outputs/MPD_clades.txt"
)

capture.output(
  disparity_biome$SV_summary,
  file = "3.outputs/SV_biomes.txt"
)

capture.output(
  disparity_biome$SR_summary,
  file = "3.outputs/SR_biomes.txt"
)

capture.output(
  disparity_biome$MPD_summary,
  file = "3.outputs/MPD_biomes.txt"
)

