################################################################################
# 02 - NMDS sobre a distancia de Gower (dados mistos: continuos + categoricos)
# Requer: rodar 01_preparo_dados.R antes
################################################################################

library(dplyr)
library(tibble)
library(ggplot2)
library(vegan)   # metaMDS, ordihull, stressplot

setwd("~/Documents/GitHub/bvasconcelos-IC-disparidade-floral/")

# ---- ler dados --------------------------------------------------------------
traits <- readRDS("3.outputs/flower_traits_matrix.rds")
meta   <- readRDS("3.outputs/meta.rds")
gower  <- readRDS("3.outputs/flower_gower_distance.rds")

# ---- NMDS -------------------------------------------------------------------
set.seed(1456)
nmds <- metaMDS(gower, k = 2, trymax = 1000)

# Stress = desajuste entre distancias na ordenacao e dissimilaridades originais.
# Guia: <0.1 bom, 0.1-0.2 aceitavel, >0.2 usar com cautela.
nmds$stress
# com trymax = 1000: 0.1705168

# Shepard plot: dissimilaridades observadas (x) vs distancias na ordenacao (y).
# Pontos proximos da linha monotona = bom ajuste.
stressplot(nmds)

# ---- stress vs numero de dimensoes -----------------------------------------
stress.vec <- vector()
for (i in 1:10) {
  nmds_k <- metaMDS(gower, k = i)
  stress.vec[i] <- nmds_k$stress
}
plot(stress.vec, pch = 19, type = "b",
     ylab = "Stress", xlab = "Number of Dimensions")

# ---- scores + clados --------------------------------------------------------
scores_nmds <- as.data.frame(vegan::scores(nmds, display = "sites")) %>%
  rownames_to_column("species")

nmds_df <- scores_nmds %>%
  left_join(meta, by = "species")

table(nmds_df$clade, useNA = "ifany") # alguma especie sem clado?

# TODO: considerar excluir o outgroup
# nmds_df <- nmds_df %>% filter(clade != "outgroup")

clade_levels <- unique(nmds_df$clade)
pal <- setNames(
  hcl.colors(length(clade_levels), palette = "Dark 3"),
  clade_levels
)

# ---- morfoespaco (ggplot) ---------------------------------------------------
ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = clade)) +
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

# ---- colorido por carater (ex.: glandula nas anteras) -----------------------
nmds_traits <- nmds_df %>%
  left_join(rownames_to_column(traits, "species"), by = "species")

ggplot(nmds_traits, aes(x = NMDS1, y = NMDS2, color = anther_gland_presence)) +
  geom_point(size = 2.5, alpha = 0.85) +
  theme_classic() +
  coord_equal() +
  labs(color = "Anther gland")

# ---- morfoespaco com hulls por clado (base R + vegan) -----------------------
# reordenar metadados para corresponder aos scores
meta_nmds <- meta[match(rownames(nmds$points), meta$species), ]

stopifnot(identical(as.character(meta_nmds$species), rownames(nmds$points)))

groups <- factor(meta_nmds$clade)

clade_levels <- levels(groups)
pal <- setNames(
  hcl.colors(length(clade_levels), palette = "Dark 3"),
  clade_levels
)
point_cols <- unname(pal[as.character(groups)])

plot(
  nmds,
  type = "n",  # plot vazio; camadas adicionadas manualmente
  xlab = "NMDS axis 1",
  ylab = "NMDS axis 2",
  main = "Floral morphospace by clade"
)

# Hull = menor poligono convexo que contem todas as especies do clado
ordihull(
  nmds,
  groups = groups,
  draw = "polygon",
  col = adjustcolor(pal, alpha.f = 0.20),
  border = NA
)

points(nmds, display = "sites", pch = 19, col = point_cols)

legend(
  "topright",
  legend = clade_levels,
  col = pal,
  pch = 19,
  bty = "n",
  cex = 0.5
)

# ---- salvar coordenadas -----------------------------------------------------
saveRDS(nmds, "3.outputs/nmds_gower.rds")
write.csv(nmds_df, "3.outputs/nmds_scores_clade.csv", row.names = FALSE)