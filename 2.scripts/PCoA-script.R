################################################################################
# PCoA (Gower) + FILOMORFOESPACO + DISPARIDADE (SV, SR, MPD)
# a matriz de Gower nao e euclidiana -> correcao de Cailliez
#      para autovalores negativos.
################################################################################

library(dplyr)
library(tibble)
library(ggplot2)
library(ape)       # pcoa, read.tree, drop.tip
library(phytools)  # phylomorphospace
library(dispRity)
library(paletteer)

setwd("~/Documents/GitHub/bvasconcelos-IC-disparidade-floral/")

# ---- ler dados --------------------------------------------------------------
traits <- readRDS("3.outputs/flower_traits_matrix.rds")
meta   <- readRDS("3.outputs/meta.rds")
gower  <- readRDS("3.outputs/flower_gower_distance.rds")
biomes <- read.csv("3.outputs/ecological_data_mimoseae.csv", stringsAsFactors = FALSE)
tree   <- read.tree("4.trees/mimosoid_calibrated_clean_updated.tre")

# podar filogenia
tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, rownames(traits)))
tree_pruned$tip.label <- gsub(" ", "_", tree_pruned$tip.label)

# ---- 1. PCoA ----------------------------------------------------------------
pcoa_res <- ape::pcoa(gower, correction = "cailliez")
eigen <- pcoa_res$values$Corr_eig

# variancia explicada pelos 3 primeiros eixos
round(eigen[1:3] / sum(eigen[eigen > 0]) * 100, 2)
# PC1 11.48 | PC2 9.43 | PC3 6.2

flower_eigenvalues  <- pcoa_res$values
flower_eigenvectors <- pcoa_res$vectors.cor

# coordenadas corrigidas 
coords <- pcoa_res$vectors.cor

lab1 <- sprintf("PCoA 1 (%.1f%%)", 100 * pcoa_res$values$Rel_corr_eig[1])
lab2 <- sprintf("PCoA 2 (%.1f%%)", 100 * pcoa_res$values$Rel_corr_eig[2])

# ordem identica a da arvore (checagem)
scores_pcoa <- coords[match(tree_pruned$tip.label, rownames(coords)), ]
stopifnot(identical(rownames(scores_pcoa), tree_pruned$tip.label))

saveRDS(pcoa_res, "3.outputs/pcoa_gower.rds")

# ---- 2. data frame unico (joins, sem depender de ordem) ---------------------
scores_df <- as.data.frame(coords[, 1:3]) %>%
  rownames_to_column("species") %>%
  left_join(meta, by = "species") %>%
  mutate(genus = sub("_.*", "", species)) %>%
  left_join(rownames_to_column(traits, "species") %>%
              dplyr::select(species, anther_gland_presence, nectary),
            by = "species")

# ---- 3. morfoespaco simples -------------------------------------------------
# por clado
ggplot(scores_df, aes(x = Axis.1, y = Axis.2, color = clade)) +
  geom_point(size = 3) +
  theme_classic() +
  labs(x = lab1, y = lab2, color = "Clade") +
  coord_equal()

# por genero
# ggplot(scores_df, aes(x = Axis.1, y = Axis.2, color = genus)) +
#   geom_point(size = 3) +
#   theme_classic() +
#   labs(x = lab1, y = lab2) +
#   coord_equal()

# por glandula nas anteras (padrao associado ao carater)
ggplot(scores_df, aes(x = Axis.1, y = Axis.2, color = anther_gland_presence)) +
  geom_point(size = 2) +
  theme_classic() +
  labs(x = lab1, y = lab2, color = "Anther gland") +
  coord_equal()

# ---- 4. filomorfoespaco -----------------------------------------------------
# phylomorphospace exige matriz numerica pura, na mesma ordem da arvore
space_axis <- coords[tree_pruned$tip.label, c(1, 2)]

phylo_obj <- phytools::phylomorphospace(
  tree_pruned,
  space_axis,
  xlab = "PCoA1",
  ylab = "PCoA2",
  label = "off"
)

phylo_data <- data.frame(
  xstart = phylo_obj$xx[phylo_obj$edge[, 1]],
  ystart = phylo_obj$yy[phylo_obj$edge[, 1]],
  xstop  = phylo_obj$xx[phylo_obj$edge[, 2]],
  ystop  = phylo_obj$yy[phylo_obj$edge[, 2]],
  nodestart = phylo_obj$edge[, 1],
  nodestop  = phylo_obj$edge[, 2]
)

# clados (sem outgroup) para cores/formas
clados_escala <- c(
  "Adenanthera clade", "Albizia clade", "Archidendron clade", "Inga clade",
  "Samanea clade", "Jupunba clade", "Mimosa clade", "Pithecellobium clade",
  "Cojoba clade", "Calliandra clade", "Zapoteca clade", "Stryphnodendron clade",
  "Senegalia grade", "Parkia clade", "Dichrostachys clade", "Neltuma clade",
  "Prosopis clade", "Newtonia grade", "Entada clade"
)

clados_dados <- sort(unique(as.character(meta$clade)))
setdiff(clados_dados, clados_escala) # clados nos dados fora da escala (outgroup)
setdiff(clados_escala, clados_dados) # clados da escala ausentes nos dados

colors <- paletteer_dynamic("cartography::multi.pal", 20)
shapes <- setNames(
  c(16, 17, 15, 18, 3, 4, 8, 0, 1, 2, 5, 6, 7, 9, 10, 11, 12, 13, 14),
  clados_escala
)

data_plot <- scores_df %>%
  filter(species %in% tree_pruned$tip.label, clade != "outgroup") %>%
  rename(PCoA1 = Axis.1, PCoA2 = Axis.2, Clade = clade)

phylomorphospace_plot <-
  ggplot() +
  geom_segment(
    data = phylo_data,
    aes(x = xstart, y = ystart, xend = xstop, yend = ystop),
    linewidth = 0.25, colour = "grey75"
  ) +
  geom_point(
    data = data_plot,
    aes(x = PCoA1, y = PCoA2, colour = Clade, shape = Clade),
    size = 2, alpha = 0.85
  ) +
  scale_colour_manual(values = colors, name = "Clade", drop = FALSE) +
  scale_shape_manual(values = shapes, name = "Clade", drop = FALSE) +
  coord_fixed(ratio = 1) +
  labs(x = lab1, y = lab2) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9),
    legend.key.height = grid::unit(0.4, "cm"),
    axis.text = element_text(colour = "black")
  )

phylomorphospace_plot

# ---- grupos (clados e biomas) --------------------------------------------
set.seed(1234)
bootstrap     <- 1000
rarefaction   <- "min"          # (no script de Hill-Smith voce usou "min")
drop_clades   <- c("outgroup")

# padronizar nomes (underscore)
meta$species   <- gsub(" ", "_", meta$species)
biomes$species <- gsub(" ", "_", biomes$species)
rownames(coords) <- gsub(" ", "_", rownames(coords))
rownames(gower)  <- gsub(" ", "_", rownames(gower))
colnames(gower)  <- gsub(" ", "_", colnames(gower))

# cat("Species in scores:          ", nrow(coords), "\n")
# cat("Species in scores w/ clade: ", sum(rownames(coords) %in% meta$species), "\n")
# cat("Species in scores w/ biome: ", sum(rownames(coords) %in% biomes$species), "\n")

# clados
groups_clade <- split(meta$species, meta$clade)
groups_clade <- groups_clade[!names(groups_clade) %in% drop_clades]

# biomas: uma especie pode ocorrer em mais de um bioma (grupos se sobrepoem)
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

# podar para especies com dados de bioma
gower_biomes <- gower[rownames(gower) %in% biomes$species,
                      rownames(gower) %in% biomes$species]

groups_biome <- lapply(groups_biome, function(sp) intersect(sp, rownames(gower_biomes)))
groups_biome <- groups_biome[lengths(groups_biome) > 0]

sapply(groups_biome, length) # tamanho dos grupos

# ---- subsets dispRity ----------------------------------------------------
library(dispRity)
# SV e SR -> coordenadas da PCoA | MPD -> matriz de distancia (Gower)

sub_clade_pcoa  <- custom.subsets(coords, groups_clade)
sub_clade       <- custom.subsets(gower, groups_clade, dist.data = TRUE)

sub_biomes_pcoa <- custom.subsets(coords, groups_biome)
sub_biomes      <- custom.subsets(gower_biomes, groups_biome, dist.data = TRUE)

# Metricas:
#   SV  = soma das variancias (dispersao em torno do centroide)
#   SR  = soma das amplitudes (extensao ocupada; sensivel a outliers)
#   MPD = distancia media par-a-par entre especies do grupo
# Bhattacharyya: sobreposicao entre distribuicoes de disparidade (0-1; 1 = total)
# Por que rarefazer: SV e principalmente SR crescem com n; rarefazer todos os
# grupos ao mesmo n torna-os comparaveis.
# CUIDADO: com muitos bootstraps os p-valores de Wilcoxon ficam minusculos
# quase automaticamente; use como descritivo e olhe IC + Bhattacharyya.

# ---- 7. CLADOS --------------------------------------------------------------
bootstrapped_clades <- boot.matrix(sub_clade_pcoa,
                                   bootstraps = bootstrap,
                                   rarefaction = rarefaction)
#prosopis and samanea clade has less than 3 elemtns

# SV
sv_clade         <- dispRity(bootstrapped_clades, metric = c(sum, variances))
clade_sv_smry    <- summary(sv_clade)
clade_sv_wilcox  <- test.dispRity(sv_clade, wilcox.test, correction = "holm")
clade_sv_bhatt   <- test.dispRity(sv_clade, bhatt.coeff)

# SR
sr_clade         <- dispRity(bootstrapped_clades, metric = c(sum, ranges))
clade_sr_smry    <- summary(sr_clade)
clade_sr_wilcox  <- test.dispRity(sr_clade, wilcox.test, correction = "holm")
clade_sr_bhatt   <- test.dispRity(sr_clade, bhatt.coeff)

# MPD (bootstrap sobre a matriz de distancia)
bootstrapped_clade_mpd <- boot.matrix(sub_clade,
                                      bootstraps = bootstrap,
                                      boot.by = "dist",
                                      rarefaction = rarefaction)

mpd_clade <- dispRity(bootstrapped_clade_mpd, metric = c(mean, pairwise.dist))
clade_mpd_smry   <- summary(mpd_clade)
clade_mpd_wilcox <- test.dispRity(mpd_clade, wilcox.test, correction = "holm")

# MPD direto da matriz de Gower (obs)
sapply(groups_clade, function(sp) {
  m <- gower[sp, sp]
  mean(m[lower.tri(m)])
})

# ---- 8. BIOMAS --------------------------------------------------------------
bootstrapped_biomes <- boot.matrix(sub_biomes_pcoa,
                                   bootstraps = bootstrap,
                                   rarefaction = rarefaction)

# SV
biomes_sv         <- dispRity(bootstrapped_biomes, metric = c(sum, variances))
biomes_sv_smry    <- summary(biomes_sv)
biomes_sv_wilcox  <- test.dispRity(biomes_sv, wilcox.test, correction = "holm")
biomes_sv_bhatt   <- test.dispRity(biomes_sv, bhatt.coeff)

# SR
sr_biomes         <- dispRity(bootstrapped_biomes, metric = c(sum, ranges))
biomes_sr_smry    <- summary(sr_biomes)
biomes_sr_wilcox  <- test.dispRity(sr_biomes, wilcox.test, correction = "holm")
biomes_sr_bhatt   <- test.dispRity(sr_biomes, bhatt.coeff)

# MPD
bootstrapped_biomes_mpd <- boot.matrix(sub_biomes,
                                       bootstraps = bootstrap,
                                       boot.by = "dist",
                                       rarefaction = rarefaction)

mpd_biomes        <- dispRity(bootstrapped_biomes_mpd, metric = c(mean, pairwise.dist))
biomes_mpd_smry   <- summary(mpd_biomes)
biomes_mpd_wilcox <- test.dispRity(mpd_biomes, wilcox.test, correction = "holm")

# ---- salvar resultados ---------------------------------------------------
saveRDS(
  list(sv_clade = sv_clade, sr_clade = sr_clade, mpd_clade = mpd_clade,
       sv_biomes = biomes_sv, sr_biomes = sr_biomes, mpd_biomes = mpd_biomes),
  "3.outputs/disparity_pcoa_results.rds"
)
