#####
# ANALISES DE DISPARIDADE
# (B) Hill-Smith PCA (ade4::dudi.hillsmith)
#Both ordinate species in a "morphospace"; disparity is analysed from there. 

# PCA Hill-Smith =======
library(ade4)
library(adegraphics)
library(dplyr) #data wrangling (select, mutate, left_join...)
library(tibble) #rownames_to_column / column_to_rownames
library(purrr) 
library(ggrepel) #
library(ggplot2) #plots
library(cluster) #daisy(): Gower distance for mixed data
library(ape) #read.tree, drop.tip (phylogenies)
library(vegan) #metaMDS (NMDS), ordihull, stressplot 

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

## read dataset
traits <-readRDS("3.outputs/flower_traits_matrix.rds")
meta <-readRDS("3.outputs/meta.rds")
biomes <- read.csv("3.outputs/ecological_data_mimoseae.csv", stringsAsFactors = FALSE)
gower <-readRDS("3.outputs/flower_gower_distance.rds")

### run pca hill smith analysis
hs <- dudi.hillsmith(
  traits,
  scannf = T
)

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

##aqui, colorir por verticilo?
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
##ajustar os plots pra definir cor de cada variavel
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

### coloring categorical variables in the morphospace =====
## Final morphospace colored by variables state
merosity_plot <- ggplot(hs_df, aes(x = Axis1, y = Axis2, 
                  color = factor(traits$flower_merosity))) +
  geom_point(size = 2.5, alpha = 0.85) +
  labs(x = paste0("Axis 1 (", pct[1], "%)"),
       y = paste0("Axis 2 (", pct[2], "%)"),
       color = "Flower merosity") +
  coord_equal() +
  theme_classic()

anther_gland_plot <- ggplot(hs_df, aes(x = Axis1, y = Axis2, 
                  color = traits$anther_gland_presence)) +
  geom_point(size = 2.5, alpha = 0.85) +
  labs(x = paste0("Axis 1 (", pct[1], "%)"),
       y = paste0("Axis 2 (", pct[2], "%)"),
       color = "Presença de glândula nas anteras") +
  coord_equal() +
  theme_classic()

nectary_plot <- ggplot(hs_df, aes(x = Axis1, y = Axis2, 
                  color = traits$nectary)) +
  geom_point(size = 2.5, alpha = 0.85) +
  labs(x = paste0("Axis 1 (", pct[1], "%)"),
       y = paste0("Axis 2 (", pct[2], "%)"),
       color = "Presença de nectário") +
  coord_equal() +
  theme_classic()

filament_union_plot <- ggplot(hs_df, aes(x = Axis1, y = Axis2, 
                  color = traits$filament_union)) +
  geom_point(size = 2.5, alpha = 0.85) +
  labs(x = paste0("PC1 (", pct[1], "%)"),
       y = paste0("PC2 (", pct[2], "%)"),
       color = "União dos filetes") +
  coord_equal() +
  theme_classic()


#aqui na cor dos filetes, pintar a bolinha da cor deles :)
filament_color_plot <- ggplot(hs_df, aes(x = Axis1, y = Axis2, 
                  color = traits$color)) +
  geom_point(size = 2.5, alpha = 0.85) +
  labs(x = paste0("PC1 (", pct[1], "%)"),
       y = paste0("PC2 (", pct[2], "%)"),
       color = "Cor dos filetes") +
  coord_equal() +
  theme_classic()

### plot morphospace with only one clade colored?


## phylomorphospace ======
# phylomorphospace exige matrix numerica pura, na mesma ordem da arvore
tree <- read.tree("4.trees/mimosoid_calibrated_clean_updated.tre")

## prune phylogeny
tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, row.names(traits)))


tree_pruned$tip.label <- gsub(pattern = " ", "_", tree_pruned$tip.label)
space_axis <- hs$li[, c(1,2)]
space_axis <- space_axis[rownames(space_axis) %in% tree_pruned$tip.label, ]

phylomorphospace <- phytools::phylomorphospace(
  tree_pruned,
  space_axis,
  xlab = "PCA1",
  ylab = "PCA2",
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
clados_dados <- sort(unique(as.character(meta$clade)))

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

#defining data axis 
data_vectors_clade <- as.data.frame(hs$li[, c (1, 2)])

data_vectors_clade %>% 
  add_column((meta$clade), (meta$species))-> data_vectors_clade

colnames(data_vectors_clade) <- c("PCA1", "PCA2", "Clade", "species")

data_vectors_clade <- data_vectors_clade %>% column_to_rownames("species")

data_vectors_clade <-
  data_vectors_clade[row.names(data_vectors_clade) %in% tree_pruned$tip.label, ]

# Clados dos dados que não estão na escala
setdiff(clados_dados, clados_escala)
#excluir o outgroup

# Nomes da escala que não existem nos dados
setdiff(clados_escala, clados_dados)

library(paletteer)
colors <- paletteer_dynamic("cartography::multi.pal", 20)
shapes <- setNames(c(16, 17, 15, 18, 3, 4, 8, 0, 1, 2, 5, 6, 7, 9, 10, 11, 12, 13, 14),
                   clados_escala)
pct <- round(100 * hs$eig / sum(hs$eig), 1)

# plot
phylomorphospace <-
  ggplot() +
  geom_segment(data = phylo_data,
               aes(x = xstart, y = ystart, xend = xstop, yend = ystop),
               linewidth = 0.25, colour = "grey75") +
  geom_point(data = data_vectors_clade,
             aes(x = PCA1, y = PCA2, colour = Clade, shape = Clade),
             size = 2, alpha = 0.85) +
  scale_colour_manual(values = colors, name = "Clade", drop = FALSE) +
  scale_shape_manual(values = shapes, name = "Clade", drop = FALSE) +
  coord_fixed(ratio = 1) +
  labs(x = paste0("Axis 1 (", pct[1], "%)"),
       y = paste0("Axis 2 (", pct[2], "%)")) +
  theme_classic(base_size = 12) +
  theme(legend.position = "right",
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 9),
        legend.key.height = grid::unit(0.4, "cm"),
        axis.text = element_text(colour = "black"))
#check
phylomorphospace


###plot arrow morphospace
# fill in with variable names
# verticil_tab <- data.frame(
#   variable = c("flower_merosity", "stamen_count", "anther_gland_presence", "nectary", "n
# calyx_length_mean"),   # ... add all variables
#   verticil = c("Perianth",        "Androecium")      # ... Calyx, Corolla, Androecium, Gynoecium, etc.
# )
# 
# verticil_cols <- c(Calyx = "#2E7D32", Corolla = "#C2185B",
#                    Perianth = "#F9A825", Androecium = "#1565C0",
#                    Gynoecium = "#6A1B9A", Other = "grey50")

co_df <- data.frame(hs$co[, 1:4])
colnames(co_df) <- c("Axis1", "Axis2","Axis3", "Axis4")
co_df$column   <- rownames(hs$co)                 # variable or category name
co_df$variable <- rownames(hs$cr)[hs$assign]      # parent variable of each column
co_df$type     <- as.character(hs$index)[hs$assign]

# co_df <- co_df %>%
#   left_join(verticil_tab, by = "variable") %>%
#   mutate(verticil = ifelse(is.na(verticil), "Other", verticil))

pct <- round(100 * hs$eig / sum(hs$eig), 1)

# optional: show only the strongest arrows to avoid clutter
co_df$len <- sqrt(co_df$Axis1^2 + co_df$Axis2^2)
co_plot <- co_df %>% filter(len > quantile(len, 0.5))   # top 50%; adjust or remove

p_arrows_12 <- ggplot(co_plot) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_segment(aes(x = 0, y = 0, xend = Axis1, yend = Axis2),
               arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.7) +
  geom_text_repel(aes(x = Axis1, y = Axis2, label = column),
                  size = 3, max.overlaps = Inf, show.legend = FALSE) +
  labs(x = paste0("Axis 1 (", pct[1], "%)"),
       y = paste0("Axis 2 (", pct[2], "%)")) +
  coord_equal() +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom")


p_arrows_34 <- ggplot(co_plot) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_segment(aes(x = 0, y = 0, xend = Axis3, yend = Axis4),
               arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.7) +
  geom_text_repel(aes(x = Axis3, y = Axis4, label = column),
                  size = 3, max.overlaps = Inf, show.legend = FALSE) +
  labs(x = paste0("Axis 3 (", pct[3], "%)"),
       y = paste0("Axis 4 (", pct[4], "%)")) +
  coord_equal() +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom")

## plot arrows and morphospace togheter
li_df <- data.frame(hs$li[, 1:2]); colnames(li_df) <- c("Axis1", "Axis2")

mult <- 0.8 * max(abs(li_df)) / max(abs(co_plot[, c("Axis1", "Axis2")]))

p_biplot <- ggplot() +
  geom_point(data = li_df, aes(Axis1, Axis2), colour = "grey70", size = 1.5) +
  geom_segment(data = co_plot,
               aes(x = 0, y = 0, xend = Axis1 * mult, yend = Axis2 * mult),
               arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.7) +
  geom_text_repel(data = co_plot,
                  aes(x = Axis1 * mult, y = Axis2 * mult, label = column),
                  size = 3, max.overlaps = Inf, show.legend = FALSE) +
  labs(x = paste0("Axis 1 (", pct[1], "%)"), y = paste0("Axis 2 (", pct[2], "%)")) +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "bottom")

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
rarefaction <- "min" #all groups rarefied to this size 
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

# ====== SV, SR, MPD + summaries + tests =============
# Metrics:
#   SV  = sum of variances (spread around the centroid, per axis, summed)
#   SR  = sum of ranges    (extent occupied, per axis, summed; outlier-sensitive)
#   MPD = mean pairwise distance between species of the group
#bhatt coef fo probability of overlap between morphological disparity calculated for different groups
##SV for clade
library(dispRity)
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


#create subset for sv and sr analysis
sub_clade_hs <- custom.subsets(hs_full$li, groups_clade) #dataset for sv and sr
sub_biomes_hs <- custom.subsets(hs_full$li, groups_biome)


bootstrapped_clades <- boot.matrix(sub_clade_hs, 
                                   bootstraps = bootstrap, 
                                   rarefaction = rarefaction)

sv_clade <- dispRity(bootstrapped_clades, metric = c(sum, variances))

clade_sv_smry <- summary(sv_clade)

clade_sv_wilcox <- test.dispRity(sv_clade, wilcox.test, correction = "holm")

clade_sv_bhatt <- test.dispRity(sv_clade, bhatt.coeff)
#save sv analysis and summary

#calculate sum of range
sr_clade <- dispRity(bootstrapped_clades, metric = c(sum, ranges))

clade_sr_smry <- summary(sr_clade)

clade_sr_wilcox <- test.dispRity(sr_clade, 
                                 wilcox.test,
                                 correction = "holm")

clade_sr_bhatt <- test.dispRity(sr_clade, bhatt.coeff)

##calculate sv for biomes
bootstrapped_biomes <- boot.matrix(sub_biomes_hs, bootstraps = bootstrap, 
                                   rarefaction = rarefaction)

sv_biomes <- dispRity(bootstrapped_biomes, metric = c(sum, variances))

biomes_sv_smry <- summary(sv_biomes)

biomes_sv_wilcox <- test.dispRity(sv_biomes, 
                                 wilcox.test, 
                                 correction = "holm")

biomes_sv_bhatt <- test.dispRity(sv_biomes, bhatt.coeff)
#calculate sum of range
sr_biomes <- dispRity(bootstrapped_biomes, metric = c(sum, ranges))

biomes_sr_smry <- summary(sr_biomes)

biomes_sr_wilcox <- test.dispRity(sr_biomes, 
                                 wilcox.test,
                                 correction = "holm")

biomes_sr_bhatt <- test.dispRity(sr_biomes, bhatt.coeff)



#subset for mpd analysis
sub_clade <- custom.subsets(gower, groups_clade, dist.data = TRUE) #dataset for mpd
sub_biomes <- custom.subsets(gower_biomes, groups_biome, dist.data = TRUE)

#conferir tamanho dos grupos
sapply(groups_biome, length)

#adapted from casali et al 
sv_clade$disparity$`Adenanthera clade` #valores de sv do bootstrap
sv_clade$disparity$`Albizia clade` #valores de sv do bootstrap
# exemplos do violin pra mimosa e adenanthera
sv_clade$disparity$`Adenanthera clade`[[3]] #dados rarefeiros para usar no violin plot
sv_clade$disparity$`Cojoba clade`[[2]] #dados observados pois n é igual ao numero de rarefacao

#CLADE - sum of variances
sum_1 <- as.vector(sv_clade$disparity$`Adenanthera clade`[[3]])
C1 <- c(rep("Adenanthera", bootstrap))
sv_1 <- data.frame(C1, sum_1)
colnames(sv_1) <- c("Group", "Sum_of_variance")

sum_2 <- as.vector(sv_clade$disparity$`Albizia clade`[[3]])
C2 <- c(rep("Albizia", bootstrap))
sv_2 <- data.frame(C2, sum_2)
colnames(sv_2) <- c("Group", "Sum_of_variance")

sum_3 <- as.vector(sv_clade$disparity$`Archidendron clade`[[3]])
C3 <- c(rep("Archidendron", bootstrap))
sv_3 <- data.frame(C3, sum_3)
colnames(sv_3) <- c("Group", "Sum_of_variance")

sum_4 <- as.vector(sv_clade$disparity$`Calliandra clade`[[3]])
C4 <- c(rep("Calliandra", bootstrap))
sv_4 <- data.frame(C4, sum_4)
colnames(sv_4) <- c("Group", "Sum_of_variance")

sum_5 <- as.vector(sv_clade$disparity$`Cojoba clade`[[2]]) #use the observed data
C5 <- c(rep("Cojoba", bootstrap))
sv_5 <- data.frame(C5, sum_5)
colnames(sv_5) <- c("Group", "Sum_of_variance")

sum_6 <- as.vector(sv_clade$disparity$`Dichrostachys clade`[[3]])
C6 <- c(rep("Dichrostachys", bootstrap))
sv_6 <- data.frame(C6, sum_6)
colnames(sv_6) <- c("Group", "Sum_of_variance")

sum_7 <- as.vector(sv_clade$disparity$`Entada clade`[[3]])
C7 <- c(rep("Entada", bootstrap))
sv_7 <- data.frame(C7, sum_7)
colnames(sv_7) <- c("Group", "Sum_of_variance")

sum_8 <- as.vector(sv_clade$disparity$`Inga clade`[[3]])
C8 <- c(rep("Inga", bootstrap))
sv_8 <- data.frame(C8, sum_8)
colnames(sv_8) <- c("Group", "Sum_of_variance")

sum_9 <- as.vector(sv_clade$disparity$`Jupunba clade`[[3]])
C9 <- c(rep("Jupunba", bootstrap))
sv_9 <- data.frame(C9, sum_9)
colnames(sv_9) <- c("Group", "Sum_of_variance")

sum_10 <- as.vector(sv_clade$disparity$`Mimosa clade`[[3]])
C10 <- c(rep("Mimosa", bootstrap))
sv_10 <- data.frame(C10, sum_10)
colnames(sv_10) <- c("Group", "Sum_of_variance")

sum_11 <- as.vector(sv_clade$disparity$`Neltuma clade`[[3]])
C11 <- c(rep("Neltuma", bootstrap))
sv_11 <- data.frame(C11, sum_11)
colnames(sv_11) <- c("Group", "Sum_of_variance")

sum_12 <- as.vector(sv_clade$disparity$`Newtonia grade`[[3]])
C12 <- c(rep("Newtonia", bootstrap))
sv_12 <- data.frame(C12, sum_12)
colnames(sv_12) <- c("Group", "Sum_of_variance")

sum_13 <- as.vector(sv_clade$disparity$`Parkia clade`[[3]])
C13 <- c(rep("Parkia", bootstrap))
sv_13 <- data.frame(C13, sum_13)
colnames(sv_13) <- c("Group", "Sum_of_variance")

sum_14 <- as.vector(sv_clade$disparity$`Pithecellobium clade`[[3]])
C14 <- c(rep("Pithecellobium", bootstrap))
sv_14 <- data.frame(C14, sum_14)
colnames(sv_14) <- c("Group", "Sum_of_variance")

# sum_15 <- as.vector(sv_clade$disparity$`Prosopis clade`[[3]])
# C15 <- c(rep("Prosopis", bootstrap))
# sv_15 <- data.frame(C15, sum_15)
# colnames(sv_15) <- c("Group", "Sum_of_variance")

sum_16 <- as.vector(sv_clade$disparity$`Samanea clade`[[2]]) #use the observed data
C16 <- c(rep("Samanea", bootstrap))
sv_16 <- data.frame(C16, sum_16)
colnames(sv_16) <- c("Group", "Sum_of_variance")

sum_17 <- as.vector(sv_clade$disparity$`Senegalia grade`[[3]])
C17 <- c(rep("Senegalia", bootstrap))
sv_17 <- data.frame(C17, sum_17)
colnames(sv_17) <- c("Group", "Sum_of_variance")

sum_18 <- as.vector(sv_clade$disparity$`Stryphnodendron clade`[[3]])
C18 <- c(rep("Stryphnodendron", bootstrap))
sv_18 <- data.frame(C18, sum_18)
colnames(sv_18) <- c("Group", "Sum_of_variance")

sum_19 <- as.vector(sv_clade$disparity$`Zapoteca clade`[[3]])
C19 <- c(rep("Zapoteca", bootstrap))
sv_19 <- data.frame(C19, sum_19)
colnames(sv_19) <- c("Group", "Sum_of_variance")

clades_to_plot <- rbind(sv_1, sv_2, sv_3, sv_4, sv_5, sv_6, sv_7, sv_8, sv_9, sv_10,
                        sv_11, sv_12, sv_13, sv_14, sv_16, sv_17, sv_18, sv_19)
ylab_list <- "Sum of variances"

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

clades_to_plot$Group <- factor(clades_to_plot$Group,
                              levels = c("Adenanthera","Albizia","Archidendron","Calliandra",
                                         "Cojoba","Dichrostachys","Entada","Inga","Jupunba",
                                         "Mimosa","Neltuma","Newtonia","Parkia","Pithecellobium",
                                         "Samanea","Senegalia","Stryphnodendron","Zapoteca"))

library(paletteer)
cores <-paletteer_d("ggsci::category20_d3")

##plot violin
ggplot(clades_to_plot, aes(x=Group, y=Sum_of_variance,
                           colour=Group, fill=Group)) +
  geom_violin(alpha=0.6, scale="width", width=0.85, trim=TRUE, linewidth=0.4) +
  xlab("") +
  ylab(ylab_list) +
  scale_fill_manual(values=cores) +
  scale_colour_manual(values=cores) +
  stat_summary(fun.data=ci, geom="linerange", linewidth=0.5) +
  stat_summary(fun=median, geom="point", size=2.5) +
  theme_bw()  +
  theme(legend.title=element_blank(), legend.position="bottom",
        legend.text=element_text(size=12), plot.title=element_text(size=14,face= "bold",hjust=0.5),
        axis.text.x=element_blank(), axis.text.y=element_text(size=10), axis.ticks.x=element_blank(), axis.title.x=element_text(size=12),axis.title.y=element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# ordenar clados pela mediana da disparidade (maior no topo)
ordem <- tapply(clades_to_plot$Sum_of_variance, clades_to_plot$Group, median)
clades_to_plot$Group <- factor(clades_to_plot$Group,
                               levels = names(sort(ordem)))   # menor embaixo, maior em cima

# cores fixas por clado (continuam as mesmas independente da ordem)
# cores_nomeadas <- setNames(colorRampPalette(as.vector(cores))(19),
#                            c("Adenanthera","Albizia","Archidendron","Calliandra","Cojoba",
#                              "Dichrostachys","Entada","Inga","Jupunba","Mimosa","Neltuma",
#                              "Newtonia","Parkia","Pithecellobium","Prosopis","Samanea",
#                              "Senegalia","Stryphnodendron","Zapoteca"))

clade_sv_violin <- ggplot(clades_to_plot, aes(x=Group, y=Sum_of_variance, colour=Group, fill=Group)) +
  geom_violin(alpha=0.6, scale="width", width=0.85, trim=TRUE, linewidth=0.4) +
  stat_summary(fun.data=ci, geom="linerange", linewidth=0.5, colour="black") +
  stat_summary(fun=median, geom="point", size=2.5, shape=21, colour="black") +
  scale_fill_manual(values=cores) +
  scale_colour_manual(values=cores) +
  labs(x=NULL, y="Sum of variances") +
  coord_flip() +
  theme_bw() +
  theme(legend.position="none",
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=10),
        axis.title.x=element_text(size=12),
        panel.grid.major.x=element_line(colour="grey90"),
        panel.grid.major.y=element_blank(),
        panel.grid.minor=element_blank())

ggsave("5.figuras/clade_sv_plot.pdf", clade_sv_violin,width=8, height=10)


#BIOMe
cor_bioma <- 
  
#BIOMES - sum of variances
sum_1 <- as.vector(sv_biomes$disparity$`Tropical & Subtropical Moist Broadleaf Forests`[[3]])
C1 <- c(rep("1", length(sum_1)))
sv_1 <- data.frame(C1, sum_1)
colnames(sv_1) <- c("Group", "Sum_of_variance")

sum_2 <- as.vector(sv_biomes$disparity$`Tropical & Subtropical Dry Broadleaf Forests`[[3]])
C2 <- c(rep("2", length(sum_2)))
sv_2 <- data.frame(C2, sum_2)
colnames(sv_2) <- c("Group", "Sum_of_variance")

sum_3 <- as.vector(sv_biomes$disparity$`Tropical & Subtropical Coniferous Forests`[[3]]) #use [[2]] if n = rarefaction value
C3 <- c(rep("3", length(sum_3)))
sv_3 <- data.frame(C3, sum_3)
colnames(sv_3) <- c("Group", "Sum_of_variance")

sum_4 <- as.vector(sv_biomes$disparity$`Temperate Broadleaf & Mixed Forests`[[3]])
C4 <- c(rep("4", length(sum_4)))
sv_4 <- data.frame(C4, sum_4)
colnames(sv_4) <- c("Group", "Sum_of_variance")

sum_5 <- as.vector(sv_biomes$disparity$`Temperate Conifer Forests`[[2]])
C5 <- c(rep("5", length(sum_5)))
sv_5 <- data.frame(C5, sum_5)
colnames(sv_5) <- c("Group", "Sum_of_variance")

sum_7 <- as.vector(sv_biomes$disparity$`Tropical & Subtropical Grasslands, Savannas & Shrublands`[[3]])
C7 <- c(rep("7", length(sum_7)))
sv_7 <- data.frame(C7, sum_7)
colnames(sv_7) <- c("Group", "Sum_of_variance")

sum_8 <- as.vector(sv_biomes$disparity$`Temperate Grasslands, Savannas & Shrublands`[[3]])
C8 <- c(rep("8", length(sum_8)))
sv_8 <- data.frame(C8, sum_8)
colnames(sv_8) <- c("Group", "Sum_of_variance")

sum_9 <- as.vector(sv_biomes$disparity$`Flooded Grasslands & Savannas`[[3]])
C9 <- c(rep("9", length(sum_9)))
sv_9 <- data.frame(C9, sum_9)
colnames(sv_9) <- c("Group", "Sum_of_variance")

sum_10 <- as.vector(sv_biomes$disparity$`Montane Grasslands & Shrublands`[[3]])
C10 <- c(rep("10", length(sum_10)))
sv_10 <- data.frame(C10, sum_10)
colnames(sv_10) <- c("Group", "Sum_of_variance")

sum_12 <- as.vector(sv_biomes$disparity$`Mediterranean Forests, Woodlands & Scrub`[[3]])
C12 <- c(rep("12", length(sum_12)))
sv_12 <- data.frame(C12, sum_12)
colnames(sv_12) <- c("Group", "Sum_of_variance")

sum_13 <- as.vector(sv_biomes$disparity$`Deserts & Xeric Shrublands`[[3]])
C13 <- c(rep("13", length(sum_13)))
sv_13 <- data.frame(C13, sum_13)
colnames(sv_13) <- c("Group", "Sum_of_variance")

sum_14 <- as.vector(sv_biomes$disparity$Mangroves[[3]])
C14 <- c(rep("14", length(sum_14)))
sv_14 <- data.frame(C14, sum_14)
colnames(sv_14) <- c("Group", "Sum_of_variance")

biomes_to_plot <- rbind(sv_1, sv_2, sv_3, sv_4, sv_5, sv_7, sv_8,
                        sv_9, sv_10, sv_12, sv_13, sv_14)

# replace biome codes with biome names
biomes_to_plot$Group <- biome_names[as.character(biomes_to_plot$Group)]

ylab_list <- "Sum of variances"

ci <- function(x)
{
  exit <- quantile(x, probs = c(0.025, 0.5, 0.975))
  exit <- data.frame(y = exit[[2]], ymin = exit[[1]], ymax = exit[[3]])
  return(exit)
}

median <- function(x)
{
  exit <- quantile(x, probs = c(0.5))
  names(exit) <- c("y")
  return(exit)
}

# order biomes by median disparity (highest on top after coord_flip)
ordem <- tapply(biomes_to_plot$Sum_of_variance, biomes_to_plot$Group, median)
biomes_to_plot$Group <- factor(biomes_to_plot$Group, levels = names(sort(ordem)))

# fixed colour per biome
#pal <- paletteer_d("tvthemes::gravityFalls")
pal <- paletteer_d("RColorBrewer::Set3")
cores_biomes <- setNames(colorRampPalette(as.vector(pal))(12), unname(biome_names))

# wrap long biome names on the y axis
wrap_lab <- function(x) vapply(strwrap(x, width = 28, simplify = FALSE),
                               paste, collapse = "\n", character(1))

biome_sv_violin <-ggplot(biomes_to_plot, aes(x=Group, y=Sum_of_variance,
                           colour=Group, fill=Group)) +
  geom_violin(alpha=0.6, scale="width", width=0.85, trim=TRUE, linewidth=0.4) +
  stat_summary(fun.data=ci, geom="linerange", linewidth=0.5, colour="black") +
  stat_summary(fun=median, geom="point", size=2.5, shape=21, colour="black") +
  scale_fill_manual(values=cores_biomes) +
  scale_colour_manual(values=cores_biomes) +
  scale_x_discrete(labels=wrap_lab) +
  labs(x=NULL, y=ylab_list) +
  coord_flip() +
  theme_bw() +
  theme(legend.position="none",
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10),
        axis.title.x=element_text(size=12),
        panel.grid.major.x=element_line(colour="grey90"),
        panel.grid.major.y=element_blank(),
        panel.grid.minor=element_blank())


#save plot 
ggsave("5.figuras/biomes_sv_plot.pdf", biome_sv_violin,width=8, height=10)

