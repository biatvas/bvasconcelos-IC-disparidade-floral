### REDUÇÃO DE DIMENSIONALIDADE (GOWER) E
# ANALISES DE DISPARIDADE
# Dimensionality reduction on MIXED data (continuous + categorical traits): 
# (A) Gower distance -> NMDS
# (B) Hill-Smith PCA (ade4::dudi.hillsmith)
#Both ordinate species in a "morphospace"; disparity is analysed from there. 
################################################################################
##nao considerar glandula nas anteras, 
#mas visualizar se existe algum padrao associado a esse carater 
library(dplyr) #data wrangling (select, mutate, left_join...)
library(tibble) #rownames_to_column / column_to_rownames
library(purrr) 
library(ggrepel) #
library(ggplot2) #plots
library(cluster) #daisy(): Gower distance for mixed data
library(ape) #read.tree, drop.tip (phylogenies)
library(vegan) #metaMDS (NMDS), ordihull, stressplot 

setwd("~/Documents/GitHub/bvasconcelos-IC-disparidade-floral/")

#read phylogenetic tree and ecological data
traits <- read.csv("3.outputs/mimoseae_flower_data.csv") #floral traits
biomes <- read.csv("3.outputs/ecological_data_mimoseae.csv", stringsAsFactors = FALSE) #ecological data 
tree <- read.tree("4.trees/mimosoid_calibrated_clean_updated.tre")

## prune phylogeny
tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, traits$species))

# salvar informaçao de clados ANTES de remover clade da tabela dos caracteres
# meta is a species clade lookup table, used later to colour plots and
# define groups for disparity analysis 

meta <- traits %>%
  dplyr::select(species, clade) %>%
  dplyr::distinct(species, .keep_all = TRUE)

# Prepare trait matrix ----------
traits <- traits %>%
  dplyr::select(-clade, -perianth_variable) %>%
  as.data.frame() %>%
  column_to_rownames("species")

continuous_cols <- c( grep("_mean$", colnames(traits), value = TRUE)) 
# Definindo quais são as colunas categóricas 
# categorical_cols <- setdiff(colnames(traits), continuous_cols)
# Log-transform continuous traits. Why: floral sizes are right-skewed and 
# multiplicative, so log makes differences proportional (a 1 mm vs 2 mm)
# difference matters more than 20 mm vs 21 mm)

# log nos traços contínuos
traits_log <- traits %>%
mutate(across(all_of(continuous_cols), log))

#this matters because daisy () and dudi.hillsmith() decide how to treat 
#each column by its R class (numeric = quantitative, factor = categorical)
#transformar character em factor
categorical_cols <- c(
  "filament_color",
  "filament_union",
  "nectary",
  "anther_gland_presence"
)


traits_log[categorical_cols] <- lapply(
  traits_log[categorical_cols],
  factor
)

traits <- traits_log 

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
# for each polymorphic cell ("A&B), sample one state.
resolve_polymorphisms <- function(data, columns) {
  
  data_resolved <- data
  
  for (col in columns) {
    
    x <- as.character(data_resolved[[col]])
    
    #rows with a polymorphism
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

set.seed(123)
traits <- resolve_polymorphisms(data = traits, columns = "flower_merosity")
traits <- resolve_polymorphisms(data = traits, columns = "filament_color")

#resolve stamen count with non-integers values (imputed as mean or according to phylogeny, thus, continuous values)
traits$stamen_count <- round(traits$stamen_count, digits = 0)
#traits$stamen_count <- as.numeric(as.character((traits$stamen_count)))

traits$flower_merosity <- as.numeric((as.character(traits$flower_merosity)))

str(traits) #check classes: numeric and factor
unique(traits$flower_merosity) #check the merosity values

# transforming names of colors into numbers
traits$filament_color <- as.factor(as.integer(factor(traits$filament_color)) - 1)
#yellow: 3; white: 2; pink:1; bicolored: 0
traits <- rename(traits,  color = filament_color)

str(traits)

#################################################
##Gower distance + NMDS =====------ 
###################################################
#Gower distance handles mixed data: each variable gets a dissimilarity in
# [0,1] (range-normalised absolute difference for numeric; 0/1 mismatch for
# factors) and the final distance is the average across variables.

gower <- as.matrix(daisy(traits, metric = "gower"))

saveRDS(
  gower,
  "3.outputs/flower_gower_distance.rds")

## checagem: quantos pares tem NA na distancia (caso alguma linha nao compartilhe
## nenhuma variavel observada com outra - daria distancia NA)
#como fiz a imputação vai dar 0 
#sum(is.na(as.matrix(gower)))

### NMDS =====------ 
set.seed(1456)

nmds <- metaMDS(gower, k=2, trymax=1000)

# Stress = mismatch between ordination distances and original dissimilarities.
# Rough guide: <0.1 good, 0.1-0.2 acceptable, >0.2 use with caution.
nmds$stress
#with trymax = 1000 we obtained value equal  0.169
#value equal 0.1705168

# Shepard plot: observed dissimilarities (x) vs ordination distances (y).
# Points hugging the step-like monotone regression line = good fit.
# ("Ordination distances" = Euclidean distances between species in the 2D plot.)
stressplot(nmds)

stress.vec<-vector()
for (i in 1:10){
  nmds_2d <-metaMDS(gower,k=i)
  p<-stressplot(nmds)
  rsquared<-summary(lm(p$y~p$x))
  stress.vec[i]<-nmds$stress
}

plot(stress.vec,pch=19,ylab="Stress",xlab="Number of Dimensions")

# Extract species scores and attach clade
scores_nmds <- as.data.frame(vegan::scores(nmds, display = "sites")) %>%
  rownames_to_column("species")

#create a dataframe with nmds results
nmds_df <- scores_nmds %>%
  left_join(meta, by = "species")

table(nmds_df$clade, useNA = "ifany") #any species without clade?

#acho que vou excluir tbm o outgroup
#definir cores para os clados
clade_levels <- unique(nmds_df$clade)

pal <- setNames(
  hcl.colors(length(clade_levels), palette = "Dark 3"),
  clade_levels
)

## ---- Morphospace plot (ggplot) ---------
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
    title = "Floral morphospace of Mimoseae"
  )

##---- NMDS morphospace with convex hulls per clade (base R + vegan) ------- 
# Reordenar os metadados para corresponder aos scores
meta_nmds <- meta[match(rownames(nmds$points), meta$species),]

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
  clade_levels)

point_cols <- unname(pal[as.character(groups)])

# Plot
plot(
  nmds,
  type = "n",    #empty plot and layers manually
  xlab = "NMDS axis 1",
  ylab = "NMDS axis 2",
  main = "Floral morphospace by clade"
)

# Hull = smallest convex polygon enclosing all species of a clade
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
  cex = 0.5
)


#it plots the observed dissimilarities against distances in ordination space,
#and also shows the non-metric monotone regression.
#distances in ordination space???????

### save NMDS coordinates + clade as text
# 

# (B) PCA Hill-Smith =======
library(ade4)
library(adegraphics)

# Hill-Smith needs NO distance matrix: it works directly on the data table.
#   - numeric columns   -> treated as quantitative (centred and scaled)
#   - factor columns    -> treated as categorical (dummy-coded, weighted)
# Quantitative and categorical variables are balanced so each variable
# contributes comparably to total inertia (variance). The result is a
# "PCA" where both kinds of variable are represented.
#
# Therefore: every column of `traits` must be correctly typed (numeric vs factor).
#
# scannf = FALSE -> don't ask interactively how many axes to keep
# nf = 5         -> keep the first 5 axes (eigenvalues are still computed for all)
hs <- dudi.hillsmith(
  traits,
  scannf = T
)

# considerando flower_merosity as factor
# traits_2 <- traits
# traits_2$flower_merosity <- as.factor(traits_2$flower_merosity)
# hs_test <- dudi.hillsmith(
#   traits_2,
#   scannf = F,
#   nf = 13
# )
# 
# cor(hs$li$Axis1, hs_test$li$Axis1)

#screeplot(hs, main = "variancia por eixo")
summary(hs)

# The number of axes = (n quantitative vars) + sum over factors of (levels - 1).
# (It is NOT "one axis per category of a quantitative variable": it is the
# categories of the FACTORS that add axes.)

#autovalores
hs$eig #variancia generalizada explicada por eixo expalhada

axes_contribution <-  100*hs$eig/sum(hs$eig)
axes_cumulative <- cumsum(axes_contribution)

data.frame(
  eixo = seq_along(hs$eig),
  variancia = axes_contribution,
  acumulada = axes_cumulative
)

# Contribution of each VARIABLE to each axis.
# For a quantitative variable: squared correlation with the axis.
# For a factor: correlation ratio (eta^2) = how much the factor's categories
# separate species along that axis. Both range from 0 to 1, so they are comparable.
hs$cr # quanto cada variavel esta associada a cada eixo

# Type of each variable as understood by dudi.hillsmith
# ("q" = quantitative, "f" = factor). Use this to confirm that columns
# like flower_merosity / stamen_count were read as intended.
hs$index # tipos de cada var

scatter(hs)  # quick overview plot ade4
screeplot(hs)
s.arrow(hs$co[1:2])
s.label(hs$li, labels = NULL) 
# A flat/elongated cloud usually reflects that variance is spread over many
# axes (each one explains little), i.e. no single dominant gradient.

hs$cr[1]  #importancia no eixo 1
#importancia de cada variavel para cada eixo. 
# tamanho floral com maior peso. podemos padronizar as var. continuas
# pela media geometrica 

hs$c1 #laodings (autovetores) 
hs$co[1] #c1 reescalonado pelos autovalores.

# Quick plot of species scores on axes 1 and 2 (hs$li = row/species coordinates)
plot(hs$li[, 1], hs$li[, 2], xlab = "PC1", ylab = "PC2", pch = 19)
# text(hs$li[,1], hs$li[,2], labels = rownames(traits), pos = 1)


# ---- Attach clade to Hill-Smith scores --------------------------------------
# The species -> clade table is `meta`.
hs_df <- hs$li %>%
  as.data.frame() %>%
  rownames_to_column("species") %>%
  left_join(meta %>% dplyr::select(species, clade), by = "species")

table(hs_df$clade, useNA = "ifany")   # every species should have a clade

## Final morphospace plot (Hill-Smith)
ggplot(hs_df, aes(x = Axis1, y = Axis2, color = clade)) +
  geom_point(size = 2.5, alpha = 0.85) +
  coord_equal() +
  theme_classic() +
  labs(
    x = sprintf("Axis 1 (%.1f%%)", axes_contribution[1]),
    y = sprintf("Axis  2 (%.1f%%)", axes_contribution[2]),
    color = "Clado",
    title = "Floral Morphospace of tribe Mimoseae",
    subtitle = "Hill-Smith ordenation")

## Final morphospace colored by variables state
ggplot(hs_df, aes(x = Axis1, y = Axis2, 
                  color = factor(traits$flower_merosity))) +
  geom_point(size = 2.5, alpha = 0.85) +
  coord_equal() +
  theme_classic()

ggplot(hs_df, aes(x = Axis1, y = Axis2, 
                  color = traits$anther_gland_presence)) +
  geom_point(size = 2.5, alpha = 0.85) +
  coord_equal() +
  theme_classic()

ggplot(hs_df, aes(x = Axis1, y = Axis2, 
                  color = traits$nectary)) +
  geom_point(size = 2.5, alpha = 0.85) +
  coord_equal() +
  theme_classic()

ggplot(hs_df, aes(x = Axis1, y = Axis2, 
                  color = traits$filament_union)) +
  geom_point(size = 2.5, alpha = 0.85) +
  coord_equal() +
  theme_classic()

ggplot(hs_df, aes(x = Axis1, y = Axis2, 
                  color = traits$filament_color)) +
  geom_point(size = 2.5, alpha = 0.85) +
  coord_equal() +
  theme_classic()

### plot morphospace with only one clade colored


##DISPARITY INDICES (SV, SR, MPD) ON HILL-SMITH PCA SCORES ======
# Groups: clades and biomes
#
# Run AFTER the first script. It needs these objects in the environment:
#   traits : data.frame, species in rownames, numeric/factor columns correctly
#            typed and polymorphisms already resolved (input of dudi.hillsmith)
#   meta   : data.frame with columns `species` and `clade`
#   biomes : data.frame with columns `species` and `biome` (codes such as
#            "1", "2", "1&7", "10, 12"...)

set.seed(1234)

## ---- settings -------
bootstrap <- 1000 #defining bootstrap replicates
rarefaction <- 3 #all groups rarefied to this size 
var_threshold <- 0.90   # keep as many PCA axes as needed to reach 90% variance
drop_clades   <- c("outgroup")   # clades to exclude from the clade analysis

# ---- 1. Hill-Smith scores to use as the morphospace -------------------------
# The first script kept only 5 axes (nf = 5), which is only ~43% of variance.
# Disparity should be computed on enough axes to represent the morphospace,
# so we recompute the analysis keeping all axes and then choose a cutoff.
hs_tmp  <- dudi.hillsmith(traits, scannf = FALSE, nf = 2) 
hs_full <- dudi.hillsmith(traits, scannf = FALSE, nf = hs_tmp$rank)

cum_var <- cumsum(hs_full$eig) / sum(hs_full$eig)
n_axes  <- which(cum_var >= var_threshold)[1]
message(n_axes, " axes retained (", round(100 * cum_var[n_axes], 1),
        "% of total variance)")

# Species x axes matrix. THIS (coordinates), not the Gower distance matrix, is
# what variances/ranges must be computed on.
scores_hs <- as.matrix(hs_full$li[, 1:n_axes])

# ---- 2. Standardise species names everywhere (underscore format) ------------
rownames(scores_hs) <- gsub(" ", "_", rownames(scores_hs))
meta$species     <- gsub(" ", "_", meta$species)
biomes$species   <- gsub(" ", "_", biomes$species)

# Quick check of name matching
cat("Species in scores:           ", nrow(scores_hs), "\n") #219
cat("Species in scores w/ clade:  ", sum(rownames(scores_hs) %in% meta$species), "\n") #219
cat("Species in scores w/ biome:  ", sum(rownames(scores_hs) %in% biomes$species), "\n")
#168 species in biome

## define groups 
# Clades
groups_clade <- split(meta$species, meta$clade) 
groups_clade <- groups_clade[!names(groups_clade) %in% drop_clades]

## Biomes
# a species occur in more than one biome, so groups overlap.
#split each cell into its individual numeric codes and match exactly.
## colorir por biomas 

biome_names <- c(
  "1"  = "Tropical & Subtropical Moist Broadleaf Forests",
  "2"  = "Tropical & Subtropical Dry Broadleaf Forests",
  "3"  = "Tropical & Subtropical Coniferous Forests",
  "4"  = "Temperate Broadleaf & Mixed Forests",
  "5"  = "Temperate Conifer Forests",
  "7"  = "Tropical & Subtropical Grasslands, Savannas & Shrublands",
  "8"  = "Temperate Grasslands, Savannas & Shrublands",
  "9"  = "Flooded Grasslands & Savannas",
  "10" = "Montane Grasslands & Shrublands",
  "12" = "Mediterranean Forests, Woodlands & Scrub",
  "13" = "Deserts & Xeric Shrublands",
  "14" = "Mangroves"
)

biome_codes <- lapply(as.character(biomes$biome), function(x) {
  if (is.na(x) || x == "") return(character(0))
  setdiff(strsplit(x, "[^0-9]+")[[1]], "")
})

groups_biome <- lapply(names(biome_names), function(code) {
  unique(biomes$species[vapply(biome_codes, function(x) code %in% x, logical(1))])
})
names(groups_biome) <- unname(biome_names)

# Keep only species present in the score matrix and groups with >= min_n.
# Why rarefy: SV, and especially SR, grow with sample size. Without rarefaction,
# large clades/biomes look more "disparate" simply because they have more
# species. Rarefying every group to the same n makes them comparable.

##calculate also PCoA index
#we can compare after how pcoa and PCA hill smith differ from coordinates values and morphospace distribution 
## PCoA com a matriz de gower =====================================
#como a matriz nao é euclidiana, 
# talvez a gente possa considerar uma correção para autovalores negativos
pcoa <- pcoa(gower, correction = "cailliez")  #ape package

## calculating variance explained explained by each principal component
eigen <- pcoa$values$Corr_eig

#export and summarize PCoA eigenvalues and vectors
flower_eigenvalues<-pcoa$values
flower_eigenvectors<-pcoa$vectors.cor

# wb<-createWorkbook()
# addWorksheet(wb, "Flower Eigenvalues")
# addWorksheet(wb, "Flower Eigenvectors")
# writeData(wb, 1, flower_eigenvalues)
# writeData(wb, 2, flower_eigenvectors)

#saveWorkbook(wb, "3.outputs/PCoA_values_and_vectors.xlsx", overwrite=TRUE)

# PC1
round(eigen[1] / sum(eigen[eigen > 0]) * 100, 2)
# 11.48
# PC2
round(eigen[2] / sum(eigen[eigen > 0]) * 100, 2)
# 9.43
# PC3
round(eigen[3] / sum(eigen[eigen > 0]) * 100, 2)
# 6.2

#salvar coordenadas
coords <- pcoa$vectors.cor

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
# traits_info <- traits %>% rownames_to_column("species") 
# traits_info$clade <- meta$clade

pcoa_df <- pcoa$vectors %>%
  as.data.frame() %>%
  rownames_to_column("species") %>%
  left_join(
    meta %>%
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
# tree_pruned$tip.label <- gsub(pattern = "_", " ", tree_pruned$tip.label)
data_vectors_clade <- as.data.frame(pcoa$vectors.cor[, c (1, 2)])

data_vectors_clade %>% 
  add_column((meta$clade), (meta$species))-> data_vectors_clade

colnames(data_vectors_clade) <- c("PCoA1", "PCoA2", "Clade", "species")

data_vectors_clade <- data_vectors_clade %>% column_to_rownames("species")

data_vectors_clade <-
  data_vectors_clade[row.names(data_vectors_clade) %in% tree_pruned$tip.label, ]

## phylomorphospace ======
# phylomorphospace exige matrix numerica pura, na mesma ordem da arvore

tree_pruned$tip.label <- gsub(pattern = " ", "_", tree_pruned$tip.label)
space_axis <- pcoa$vectors.cor[, c(1,2)]
space_axis <- space_axis[rownames(space_axis) %in% tree_pruned$tip.label, ]

phylomorphospace <- phytools::phylomorphospace(
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

library(paletteer)

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

# pcoa_data <- as.data.frame(pcoa$vectors.cor[, c(1, 2)])
# pcoa_data$species <- row.names(pcoa_data)
# pcoa_data <- merge(pcoa_data, genus, by = "species")


#plotar morfoespaço colorido por genero
ggplot(pcoa_data, aes(x = Axis.1, y = Axis.2, color = genus)) +
  geom_point(size = 3) +
  theme_classic() + 
  ggplot2::coord_equal()

##plot morphospace
ggplot(pcoa$vectors.cor[,c(1,2)], aes(x = Axis.1, y = Axis.2)) +
  geom_point(size = 2) +
  # geom_text(aes(label = row.names(pcoa$vectors.cor)), 
  #           vjust = -0.5, size = 2) +
  theme_classic()


#plot coloring traits
ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, color = traits$anther_gland_presence)) +
  geom_point(size = 2) +
  theme_classic() +
  labs(
    x = "PCoA axis 1",
    y = "PCoA axis 2",
  ) + ggplot2::coord_equal()

# ====== SV, SR, MPD + summaries + tests =============
library(dispRity)
sub_clade <- custom.subsets(gower, groups_clade, dist.data = TRUE) #dataset for mpd 
sub_clade_pcoa <- custom.subsets(pcoa$vectors.cor, groups_clade) #dataset for sv and sr

# podar especies para conter apenas que tem dados de bioma
# gower_biomes <- row.names(gower)[(which(row.names(gower) %in% biomes$species))]
gower_biomes <- gower[rownames(gower) %in% biomes$species, rownames(gower) %in% biomes$species]

# Manter apenas espécies presentes na matriz filtrada
groups_biome <- lapply(groups_biome, function(sp) {
  intersect(sp, rownames(gower_biomes))
})

# Remover grupos vazios
groups_biome <- groups_biome[
  lengths(groups_biome) > 0
]

#conferir tamanho dos grupos
sapply(groups_biome, length)

sub_biomes_pcoa <- custom.subsets(pcoa$vectors.cor, groups_biome)
sub_biomes <- custom.subsets(gower_biomes, groups_biome, dist.data = TRUE)
  
# Metrics:
#   SV  = sum of variances (spread around the centroid, per axis, summed)
#   SR  = sum of ranges    (extent occupied, per axis, summed; outlier-sensitive)
#   MPD = mean pairwise distance between species of the group
#bhatt coef fo probability of overlap between morphological disparity calculated for different groups
##SV for clade
bootstrapped_clades <- boot.matrix(sub_clade_pcoa, 
                                   bootstraps = bootstrap, 
                                   rarefaction = rarefaction)
sv_clade <- dispRity(bootstrapped_clades, metric = c(sum, variances))
clade_sv_smry <- summary(sv_clade)
clade_sv_wilcox <- test.dispRity(sv_clade, 
                                 wilcox.test,
                                 correction = "holm")
clade_sv_bhatt <- test.dispRity(sv_clade, bhatt.coeff)
#salvar sv
#sum of range
sr_clade <- dispRity(bootstrapped_clades, metric = c(sum, ranges))
clade_sr_smry <- summary(sr_clade) 

clade_sr_wilcox <- test.dispRity(sr_clade, 
                                 wilcox.test,
                                 correction = "holm")
clade_sr_bhatt <- test.dispRity(sr_clade, bhatt.coeff)
#mean pairwise distance
#for mpd data, we should use distance matrix
bootstrapped_clade_mpd <- boot.matrix(sub_clade, 
                                      bootstraps = bootstrap, 
                                      boot.by = "dist",
                                      rarefaction = rarefaction)

mpd_clade <- dispRity(bootstrapped_clade_mpd, metric = c(mean, pairwise.dist))

mean_pw <- function(m) { m <- as.matrix(m); mean(m[lower.tri(m)]) }
mpd_clade_2 <- dispRity(bootstrapped_clade_mpd, metric = mean_pw)

clade_mpd_smry <- summary(mpd_clade)
clade_mpd_wilcox <- test.dispRity(mpd_clade, 
                                  wilcox.test,
                                  correction = "holm")

# # observed (no bootstrap) from dispRity
# obs_mpd <- dispRity(sub_clade, metric = c(mean, pairwise.dist))
# summary(obs_mpd)$obs

# manual, straight from the Gower matrix
sapply(groups_clade, function(sp) {
  m <- gower[sp, sp]
  mean(m[lower.tri(m)])
})

## metrics for biomes =====
# o mpd precisa ser pela matriz de distancia e as outras metricas q usam o pcoa!!!
### sum of variances
bootstrapped_biomes <- boot.matrix(sub_biomes_pcoa, bootstraps = bootstrap, 
                                   rarefaction = rarefaction)


biomes_sv <- dispRity(bootstrapped_biomes, metric = c(sum, variances))
biomes_sv_smry <- summary(biomes_sv)
biomes_sv_wilcox <- test.dispRity(biomes_sv, 
                                      wilcox.test, 
                                      correction = "holm")
biomes_sv_bhatt <- test.dispRity(biomes_sv, bhatt.coeff)
#valores entre 0 - 1, o 1 indicando sobreposiçao total 

#sum of range
sr_biomes <- dispRity(bootstrapped_biomes, metric = c(sum, ranges))
biomes_sr_smry <- summary(sr_biomes)
biomes_sr_wilcox <- test.dispRity(sr_biomes, 
                                 wilcox.test,
                                 correction = "holm")

biomes_sr_bhatt <- test.dispRity(sr_biomes, bhatt.coeff) 

#mean pairwise distance
#generate a bootstraped dataset
bootstrapped_biomes_mpd <- boot.matrix(sub_biomes, 
                                      bootstraps = bootstrap, 
                                      boot.by = "dist",
                                      rarefaction = rarefaction)
#calculate mean pairwise distance
mpd_2 <- dispRity(sub_biomes_pcoa, metric = c(mean, pairwise.dist))

mpd_biomes <- dispRity(bootstrapped_biomes_mpd, metric = c(mean, pairwise.dist))
biomes_mpd_smry <- summary(mpd_biomes)
#wilcoxon for significance test
biomes_mpd_wilcox <- test.dispRity(mpd_biomes, 
                                  wilcox.test,
                                  correction = "holm")

# Pairwise tests between groups (Holm-corrected Wilcoxon) and overlap (Bhattacharyya).
# CAUTION: with many bootstraps the Wilcoxon p-values become tiny almost
# automatically (n is set by you). Treat them as descriptive and rely on the
# confidence intervals and the Bhattacharyya overlap too.
 
### Check code of Casali et al. (2023)
### available at https://doi.org/10.5281/zenodo.7240006
##check code of Clark et al. (2023)
#sum of variances for disparity of a group 
library(patchwork)

##STATS - summarize disprity descriptive and test statistics and export
labels<-c(rep("Cranial",7),rep("Postcranial",7),rep(c(rep("Cranial",3),rep("Postcranial",3)),2))
fam_disp<-do.call(rbind,fam_DISP)
clade_disp<-do.call(rbind,clade_DISP)
colnames(fam_disp)[8]<-"97.5%"
colnames(clade_disp)[8]<-"97.5%"
colnames(diet_DISP)[8]<-"97.5%"
colnames(loc_DISP)[8]<-"97.5%"
info_disp<-rbind(fam_disp,clade_disp,diet_DISP,loc_DISP)
options(digits = 2)
summary_DISP<-data.frame(labels,info_disp)
summary_DISP<-summary_DISP[-c(7:8)]
colnames(summary_DISP)<-c("Partition","Group","N","SV_Obs","Median_BS","0.025_BS","0.975_BS")
summary_DISP$Group[summary_DISP$Group=="Mixed_Grazer"]<-"Mixed/Grazer"

statistic_W<-c(unlist(fam_DISP_wilcoxon)[c(1:6,13:18)],unlist(clade_DISP_wilcoxon)[c(1,3)],unlist(diet_DISP_wilcoxon)[1],unlist(loc_DISP_wilcoxon)[1])
p_value<-c(unlist(fam_DISP_wilcoxon)[c(7:12,19:24)],unlist(clade_DISP_wilcoxon)[c(2,4)],unlist(diet_DISP_wilcoxon)[2],unlist(loc_DISP_wilcoxon)[2])
b_coeff<-c(unlist(fam_DISP_bhatt),unlist(clade_DISP_bhatt),unlist(diet_DISP_bhatt),unlist(loc_DISP_bhatt))
gr_lab<-c(rep(c(rownames(fam_DISP_bhatt[[1]])),2),rep("Megatherioidea : Mylodontoidea",2),"Browser : Mixed/Grazer","Scansorial : Terrestrial")
pt_lab<-c(rep("Cranial",6),rep("Postcranial",6),rep(c("Cranial","Postcranial"),2))
summary_stats<-data.frame(pt_lab,gr_lab,round(statistic_W,2),round(p_value,2),round(b_coeff,2))
colnames(summary_stats)<-c("Partition","Groups","Statistic W","P-value","Bhattacharyya Coefficient")
wb<-createWorkbook()
addWorksheet(wb, "Sum of variances")
addWorksheet(wb, "Stats")
writeData(wb, 1, summary_DISP)
writeData(wb, 2, summary_stats)
saveWorkbook(wb, "SV_stats.xlsx", overwrite=TRUE)

##PLOTS - disparity (SV) violin plots

#FAMILY
F1<-c(rep("Megatheriidae",bootstraps))
F2<-c(rep("Megalonychidae",bootstraps))
F3<-c(rep("Mylodontidae",bootstraps))
F4<-c(rep("Scelidotheriidae",bootstraps))
fams_plot_disp<-vector(mode="list",length=length(fam_disparity))
for (i in 1:length(fam_disparity))
{
  sum_Megatheriidae<-as.vector(fam_disparity[[i]]$disparity$Megatheriidae[[3]])
  sum_Megalonychidae<-as.vector(fam_disparity[[i]]$disparity$Megalonychidae[[3]])
  sum_Mylodontidae<-as.vector(fam_disparity[[i]]$disparity$Mylodontidae[[3]])
  sum_Scelidotheriidae<-as.vector(fam_disparity[[i]]$disparity$Scelidotheriidae[[2]])
  
  SV_Megatheriidae<-data.frame(F1,sum_Megatheriidae)
  colnames(SV_Megatheriidae)<-c("Group","Sum_of_variances")
  SV_Megalonychidae<-data.frame(F2,sum_Megalonychidae)
  colnames(SV_Megalonychidae)<-c("Group","Sum_of_variances")
  SV_Mylodontidae<-data.frame(F3,sum_Mylodontidae)
  colnames(SV_Mylodontidae)<-c("Group","Sum_of_variances")
  SV_Scelidotheriidae<-data.frame(F4,sum_Scelidotheriidae)
  colnames(SV_Scelidotheriidae)<-c("Group","Sum_of_variances")
  fams_to_plot<-rbind(SV_Megatheriidae,SV_Megalonychidae,SV_Mylodontidae,SV_Scelidotheriidae)
  ylab_list<-c(rep("Sum of variances",2))
  
  ci <- function(x)
  {
    exit <- quantile(x, probs = c(0.025,0.975))
    names(exit) <- c("ymin","ymax")
    return(exit) 
  }
  
  median <- function(x)
  {
    exit <- quantile(x, probs = c(0.5))
    names(exit) <- c("y")
    return(exit) 
  }	
  
  fams_to_plot$Group<-factor(fams_to_plot$Group,levels=c("Megatheriidae","Megalonychidae","Mylodontidae","Scelidotheriidae"))
  fams_plot_disp[[i]]<-ggplot(fams_to_plot, aes(x=Group, y=Sum_of_variances, colour=Group, fill=Group),alpha=0.9) +
    geom_violin(alpha=0.5, show.legend=T) + labs(x="", y="") +
    xlab("") +
    ylab(ylab_list[i]) +
    scale_fill_manual(values=c("dodgerblue3","limegreen","red","gold2"))+
    scale_colour_manual(values=c("dodgerblue3","limegreen","red","gold2")) +
    stat_summary(fun=ci, geom="line", size=0.6) +
    stat_summary(fun=median, geom="point", size=2.5) +
    theme_bw() + labs("") +
    theme(legend.title=element_blank(), legend.position="bottom", legend.text=element_text(size=12), plot.title=element_text(size=14,face= "bold",hjust=0.5),
          axis.text.x=element_blank(), axis.text.y=element_text(size=10), axis.ticks.x=element_blank(), axis.title.x=element_text(size=12),axis.title.y=element_text(size=12),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_y_continuous(labels=label_number(accuracy=0.1 ,decimal.mark="."))
}

#CLADE
M1<-c(rep("Megatherioidea",bootstraps))
M2<-c(rep("Mylodontoidea",bootstraps))
clade_plot_disp<-vector(mode="list",length=length(clade_disparity))
for (i in 1:length(clade_disparity))
{
  sum_Megatherioidea<-as.vector(clade_disparity[[i]]$disparity$Megatherioidea[[3]])
  sum_Mylodontoidea<-as.vector(clade_disparity[[i]]$disparity$Mylodontoidea[[2]])
  
  SV_Megatherioidea<-data.frame(M1,sum_Megatherioidea)
  colnames(SV_Megatherioidea)<-c("Group","Sum_of_variances")
  SV_Mylodontoidea<-data.frame(M2,sum_Mylodontoidea)
  colnames(SV_Mylodontoidea)<-c("Group","Sum_of_variances")
  clades_to_plot<-rbind(SV_Megatherioidea, SV_Mylodontoidea)
  ylab_list<-c(rep("Sum of variances",2))
  
  ci <- function(x)
  {
    exit <- quantile(x, probs = c(0.025,0.975))
    names(exit) <- c("ymin","ymax")
    return(exit) 
  }
  
  median <- function(x)
  {
    exit <- quantile(x, probs = c(0.5))
    names(exit) <- c("y")
    return(exit) 
  }	
  
  clades_to_plot$Group<-factor(clades_to_plot$Group,levels=c("Megatherioidea","Mylodontoidea"))
  clade_plot_disp[[i]]<-ggplot(clades_to_plot, aes(x=Group, y=Sum_of_variances, colour=Group, fill=Group),alpha=0.9) +
    geom_violin(alpha=0.5, show.legend=T) + labs(x="", y="") +
    xlab("") +
    ylab(ylab_list[i]) +
    scale_fill_manual(values=c("purple", "orange")) +
    scale_colour_manual(values=c("purple", "orange")) +
    stat_summary(fun=ci, geom="line", size=0.6) +
    stat_summary(fun=median, geom="point", size=2.5) +
    theme_bw() + labs("") +
    theme(legend.title=element_blank(), legend.position="bottom", legend.text=element_text(size=12), plot.title=element_text(size=14,face= "bold",hjust=0.5),
          axis.text.x=element_blank(), axis.text.y=element_text(size=10), axis.ticks.x=element_blank(), axis.title.x=element_text(size=12),axis.title.y=element_text(size=12),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_y_continuous(labels=label_number(accuracy=0.1 ,decimal.mark="."))
}

#DIET
D1<-c(rep("Browser",bootstraps))
D2<-c(rep("Mixed/Grazer",bootstraps))	
sum_Browser<-as.vector(diet_disparity$disparity$Browser[[3]])
sum_Mixed_Grazer<-as.vector(diet_disparity$disparity$Mixed_Grazer[[2]])

SV_Browser<-data.frame(D1,sum_Browser)
colnames(SV_Browser)<-c("Group","Sum_of_variances")
SV_Mixed_Grazer<-data.frame(D2,sum_Mixed_Grazer)
colnames(SV_Mixed_Grazer)<-c("Group","Sum_of_variances")
diets_to_plot<-rbind(SV_Browser, SV_Mixed_Grazer)
ylab_list<-"Sum of variances"

ci <- function(x)
{
  exit <- quantile(x, probs = c(0.025,0.975))
  names(exit) <- c("ymin","ymax")
  return(exit) 
}

median <- function(x)
{
  exit <- quantile(x, probs = c(0.5))
  names(exit) <- c("y")
  return(exit) 
}	

diets_to_plot$Group<-factor(diets_to_plot$Group,levels=c("Browser","Mixed/Grazer"))
diet_plot_disp<-ggplot(diets_to_plot, aes(x=Group, y=Sum_of_variances, colour=Group, fill=Group),alpha=0.9) +
  geom_violin(alpha=0.5, show.legend=T) + labs(x="", y="") +
  xlab("") +
  ylab(ylab_list) +
  scale_fill_manual(values=c("darkgreen", "aquamarine3")) +
  scale_colour_manual(values=c("darkgreen", "aquamarine3")) +
  stat_summary(fun=ci, geom="line", size=0.6) +
  stat_summary(fun=median, geom="point", size=2.5) +
  theme_bw() + labs(title=) +
  theme(legend.title=element_blank(), legend.position="bottom", legend.text=element_text(size=12), plot.title=element_text(size=14,face= "bold",hjust=0.5),
        axis.text.x=element_blank(), axis.text.y=element_text(size=10), axis.ticks.x=element_blank(), axis.title.x=element_text(size=12),axis.title.y=element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_y_continuous(labels=label_number(accuracy=0.1 ,decimal.mark="."))


#Final_disparity_plots - merge morphospaces and violin plots

DISP_fam1<-ggarrange(pcoa_fam[[1]],fams_plot_disp[[1]], ncol = 1, nrow = 2, legend="none", align="hv")		
DISP_clade1<-ggarrange(pcoa_clade[[1]],clade_plot_disp[[1]], ncol = 1, nrow = 2, legend="none", align="hv")	
DISP_fam2<-ggarrange(pcoa_fam[[2]],fams_plot_disp[[2]], ncol = 1, nrow = 2, legend="none", align="hv")	
DISP_clade2<-ggarrange(pcoa_clade[[2]],clade_plot_disp[[2]],ncol = 1, nrow = 2, legend="none", align="hv")		
DISP_diet<-ggarrange(pcoa_diet,diet_plot_disp, ncol = 1, nrow = 2, legend="none", align="hv")	
DISP_locs<-ggarrange(pcoa_loc,loc_plot_disp, ncol = 1, nrow = 2, legend="none", align="hv")

DISP_cranial<-ggarrange(DISP_fam1,DISP_clade1,DISP_diet,ncol = 3, nrow = 1, legend="none", align="hv",
                        labels=c("A","B","C"),font.label=list(size = 16,face="bold"))
DISP_cranial
DISP_cranial<-annotate_figure(DISP_cranial,left=text_grob("Cranial",hjust=0, face = "bold", size = 16, rot = 90))	

DISP_postcranial<-ggarrange(DISP_fam2,DISP_clade2,DISP_locs, ncol = 3, nrow = 1, legend="none", align="hv",
                            labels=c("D","E","F"),font.label=list(size = 16,face="bold"))
DISP_postcranial
DISP_postcranial<-annotate_figure(DISP_postcranial,left=text_grob("Postcranial",hjust=0, face = "bold", size = 16, rot = 90))

DISP_final<-ggarrange(DISP_cranial,DISP_postcranial, ncol = 1, nrow = 2, legend.grob=get_legend(pcoa_fam[[1]]), legend="bottom" ,align="hv")
DISP_final	
ggsave("DISP_final.pdf",DISP_final,width=8, height=10)

###################################################################CONVERGENCE#############################################################
##test convergence hypotheses
conv_hypotheses<-ReorderData(time_tree,node_edge[1:41,7:10])
cranial_states<-conv_hypotheses[,1:2]
postcranial_states<-conv_hypotheses[,3:4]
states<-list(cranial_states,postcranial_states)

conv2<-vector(mode="list",length=length(pcoa_asr_list))
for (i in 1:length(pcoa_asr_list))
{	
  tips<-ReorderData(time_tree,pcoa_asr_list[[i]]$vectors.cor[1:41,1:2])
  nodes<-pcoa_asr_list[[i]]$vectors.cor[42:81,1:2]
  
  conv<-vector(mode="list",length=length(colnames(states[[i]])))
  for (j in 1:length(colnames(states[[i]])))
  {
    st<-states[[i]][,j]
    names(st)<-rownames(states[[i]])		
    conv[[j]]<-search.conv(RR=NULL,tree=time_tree,y=tips,state=st,aceV=nodes,declust=TRUE,nsim=1000,rsim=1000,clus=.5)
  }
  names(conv)<-colnames(states[[i]])
  conv2[[i]]<-conv
}
names(conv2)<-c("cranial","postcranial")


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
# names(biomes)
# table(biomes$biome, useNA = "ifany")
# 
#one tidy data frame, built with joins (no order assumptions)
scores_df <- as.data.frame(pcoa$vectors.cor[, 1:3]) %>%
  rownames_to_column("species") %>%
  left_join(meta, by = "species") %>%
  mutate(genus = sub("_.*", "", species))

#generate a dataframe with scores from pcoa and species data 

lab1 <- sprintf("PCoA 1 (%.1f%%)", 100 * pcoa$values$Rel_corr_eig[1]) #calculate each vector percentage
lab2 <- sprintf("PCoA 2 (%.1f%%)", 100 * pcoa$values$Rel_corr_eig[2])

#Rel_corr_eig, que, nesse objeto, representa a proporção relativa corrigida de autovalores associada ao primeiro eixo.

# everything black, the focal group coloured
plot_highlight <- function(df, group_col, level, col,
                           xvar = "Axis.1", yvar = "Axis.2",
                           xlab = lab1, ylab = lab2, hull = TRUE) {
  hl  <- !is.na(df[[group_col]]) & df[[group_col]] == level
  foc <- df[hl, ]
  p <- ggplot() +
    geom_point(data = df[!hl, ], aes(.data[[xvar]], .data[[yvar]]),
               colour = "black", size = 1.4, alpha = 0.6)
  if (hull && nrow(foc) >= 3) {
    h <- foc[chull(foc[[xvar]], foc[[yvar]]), ]
    p <- p + geom_polygon(data = h, aes(.data[[xvar]], .data[[yvar]]),
                          fill = col, alpha = 0.15)
  }
  p +
    geom_point(data = foc, aes(.data[[xvar]], .data[[yvar]]),
               colour = col, size = 3) +
    coord_equal() + theme_classic(base_size = 12) +
    labs(x = xlab, y = ylab, title = level, subtitle = paste("n =", nrow(foc)))
}

#foc is focal
#hull is more than 3, so less then this number of species doesnt generate a polygon
# 3. multi-page PDF, one page per group
save_highlight_pdf <- function(df, group_col, file, min_n = 1, ...) {
  tab <- table(df[[group_col]])
  lv  <- names(tab)[tab >= min_n]
  cols <- setNames(hcl.colors(length(lv), "Dark 3"), lv)
  pdf(file, width = 6, height = 6, onefile = TRUE)
  for (l in lv) print(plot_highlight(df, group_col, l, cols[[l]], ...))
  dev.off()
}

save_highlight_pdf(scores_df, "clade", "3.outputs/morphospace_by_clade.pdf")
save_highlight_pdf(scores_df, "genus", "3.outputs/morphospace_by_genus.pdf", min_n = 2)

# # SV Boxplots =================================================================
# 
# ### This part of the script was based on the code of Casali et al. (2023) & Barduzzi et al (2024)
# ### available at https://doi.org/10.5281/zenodo.7240006
# ### required functions
# 
