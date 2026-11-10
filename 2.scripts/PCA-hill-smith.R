#####
# ANALISES DE DISPARIDADE
# Hill-Smith PCA (ade4::dudi.hillsmith)
#Both ordinate species in a "morphospace"; disparity is analysed from there. 
# this script calculates a pca hill smith -> plot morphospace -> calculate disparity metrics (mpd, sv and sr),
# -. figures -> supp tables

# PCA Hill-Smith =======
library(ade4) #dudi.hillsmith
library(adegraphics)
library(dplyr) #data wrangling (select, mutate, left_join...)
library(tibble) #rownames_to_column / column_to_rownames
library(purrr) 
library(ggrepel) #
library(ggplot2) #plots
library(cluster) #daisy(): Gower distance for mixed data
library(ape) #read.tree, drop.tip (phylogenies)
library(vegan)  #2.7-3 metaMDS (NMDS), ordihull, stressplot 
library(phytools) #phylomorphospace 
library(paletteer) #colors for plot 
library(flextable)
library(officer)

setwd("Documents/GitHub/bvasconcelos-IC-disparidade-floral/")

# scannf = FALSE -> don't ask interactively how many axes to keep
# nf = 5         -> keep the first 5 axes (eigenvalues are still computed for all)
set.seed(1234)

###settings for next analysis
bootstrap <- 1000 #defining bootstrap replicates
rarefaction <- "min" #all groups rarefied to this size 
var_threshold <- 0.90   # PCA axes kept for disparity = enough for 90% variance
drop_clades   <- c("outgroup")   # clades to exclude from the clade analysis
#drop_clades   <- "outgroup" ?

## read dataset
traits <-readRDS("3.outputs/flower_traits_matrix.rds")
meta <-readRDS("3.outputs/meta.rds")
gower <-readRDS("3.outputs/flower_gower_distance.rds")
biomes <- read.csv("3.outputs/ecological_data_mimoseae.csv", stringsAsFactors = FALSE)
tree   <- read.tree("4.trees/mimosoid_calibrated_clean_updated.tre")


#check name mattings and prune datasets
# Underscore format everywhere (traits, meta, biomes, gower, tree)
rownames(traits) <- gsub(" ", "_", rownames(traits))
meta$species     <- gsub(" ", "_", meta$species)
biomes$species   <- gsub(" ", "_", biomes$species)
gower            <- as.matrix(gower)
rownames(gower)  <- gsub(" ", "_", rownames(gower))
colnames(gower)  <- gsub(" ", "_", colnames(gower))
tree$tip.label   <- gsub(" ", "_", tree$tip.label)

# Name-matching checks
cat("Species in traits:          ", nrow(traits), "\n")
cat("Species in traits w/ clade: ", sum(rownames(traits) %in% meta$species), "\n")
cat("Species in traits w/ biome: ", sum(rownames(traits) %in% biomes$species), "\n")
cat("Species in traits w/ tree:  ", sum(rownames(traits) %in% tree$tip.label), "\n")


### run pca hill smith analysis
#the pca hill smith works directly in the traits dataset
# All axes are kept (nf = rank); plots use axes 1-4, disparity uses
# as many axes as needed to reach `var_threshold` of the variance.
hs0 <- dudi.hillsmith(traits, scannf = FALSE, nf = 2)
hs  <- dudi.hillsmith(traits, scannf = FALSE, nf = hs0$rank)

pct     <- round(100 * hs$eig / sum(hs$eig), 1)
cum_var <- cumsum(hs$eig) / sum(hs$eig)
n_axes  <- which(cum_var >= var_threshold)[1]

eig_table <- data.frame(axis = seq_along(hs$eig),
                        variance = pct,
                        cumulative = round(100 * cum_var, 1))

write.table(eig_table, "3.outputs/axis_contribution.txt")
message(n_axes, " axes retained for disparity (",
        round(100 * cum_var[n_axes], 1), "% of total variance)")

### see some informations
# Type of each variable as understood by dudi.hillsmith
# ("q" = quantitative, "f" = factor). Use this to confirm that columns
# like flower_merosity / stamen_count were read as intended.
hs$index # tipos de cada var
summary(hs)
#autovalores
hs$eig #variancia generalizada explicada por eixo expalhada
hs$cr # quanto cada variavel esta associada a cada eixo

scatter(hs)  # quick overview plot ade4
screeplot(hs)
hs$cr[1]  #importancia no eixo 1

#importancia de cada variavel para cada eixo. 
# tamanho floral com maior peso. 
# podemos padronizar as var. continuas
# pela media geometrica?
hs$c1 #laodings (autovetores) 
hs$co[1] #c1 reescalonado pelos autovalores.

#define coordinates for disparity 
scores_hs <- as.matrix(hs$li[, 1:n_axes])

#make a dataset with clade info to future plots
hs_df <- hs$li[, 1:4] %>%
  as.data.frame() %>%
  rownames_to_column("species") %>%
  left_join(meta %>% dplyr::select(species,clade) %>% distinct(), by = "species")

table(hs_df$clade, useNA = "ifany")   

##aqui, colorir por verticilo?
s.arrow(hs$co[1:2], xax = "Comp1", yax = "Comp2") 
s.label(hs$li, labels = NULL) 

s.arrow(hs$co[3:4])

# A flat/elongated cloud usually reflects that variance is spread over many
# axes (each one explains little), i.e. no single dominant gradient.

## 3. GROUPS ========================================================
# ---- clades ----
groups_clade <- split(meta$species, meta$clade)
groups_clade <- groups_clade[!names(groups_clade) %in% drop_clades]
for (g in names(groups_clade)) {
  groups_clade[[g]] <- intersect(groups_clade[[g]], rownames(scores_hs))
}
groups_clade <- groups_clade[lengths(groups_clade) > 0]

# ---- biomes ----
# A species can occur in >1 biome (e.g. "1&7", "10, 12"), so groups overlap.
# Each code is matched exactly as a whole number inside the cell
# (so "1" does not match "10" or "12").
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

groups_biome <- list()
for (code in names(biome_names)) {
  hit <- grepl(paste0("(^|[^0-9])", code, "($|[^0-9])"), biomes$biome)   # NA -> FALSE
  groups_biome[[biome_names[[code]]]] <- intersect(unique(biomes$species[hit]),
                                                   rownames(scores_hs))
}
groups_biome <- groups_biome[lengths(groups_biome) > 0]

# Group sizes. Why rarefy: SV and especially SR grow with sample size, so
# without rarefaction large groups look more disparate just by having more species.
lengths(groups_clade)
lengths(groups_biome)
message("Rarefaction level (smallest group): ",
        min(lengths(groups_clade), lengths(groups_biome)), " species")



##  Morphospace plot (Hill-Smith) ========
##ajustar os plots pra definir cor de cada variavel
# ---- palettes (one named palette per clade, reused in every figure) ----
clade_names  <- sort(setdiff(unique(as.character(meta$clade)), drop_clades))
pal20        <- as.vector(paletteer_d("ggsci::category20_d3"))
clade_cols   <- setNames(pal20[seq_along(clade_names)], clade_names)
clade_shapes <- setNames(
  rep_len(c(16, 17, 15, 18, 3, 4, 8, 0, 1, 2, 5, 6, 7, 9, 10, 11, 12, 13, 14),
          length(clade_names)),
  clade_names)
okabe <- unname(palette.colors(9, "Okabe-Ito")[-1])


#plot phylomorphospace
tree_pruned <- keep.tip(tree, intersect(tree$tip.label, rownames(hs$li)))
space_axis  <- hs$li[tree_pruned$tip.label, 1:2]

#pdf(NULL)   
# phylomorphospace draws a base plot; we only need its coordinates
phylo_obj <- phytools::phylomorphospace(tree_pruned, space_axis, label = "off")
#dev.off()

phylo_df <- data.frame(
  xstart = phylomorphospace$xx[phylomorphospace$edge[, 1]],
  ystart = phylomorphospace$yy[phylomorphospace$edge[, 1]],
  xstop = phylomorphospace$xx[phylomorphospace$edge[, 2]],
  ystop = phylomorphospace$yy[phylomorphospace$edge[, 2]],
  nodestart = phylomorphospace$edge[, 1],
  nodestop = phylomorphospace$edge[, 2]
)


lim_max <- max(abs(c(phylo_df$xstart, phylo_df$xstop,
                             phylo_df$ystart, phylo_df$ystop, 
                     hs_df$Axis1, hs_df$Axis2 )))
                   

# Round limits to multiples of 2
lim_max <- ceiling(lim_max / 2) * 2
lim <- c(-lim_max, lim_max)

# Common axis breaks
breaks <- seq(-lim_max, lim_max, by = 2)

# plot morphospaces
morpho_scales <- list(
  scale_x_continuous(breaks = breaks, limits = lim, expand = c(0, 0)),
  scale_y_continuous(breaks = breaks, limits = lim, expand = c(0, 0)),
  coord_equal(),
  labs(x = paste0("PC1 (", pct[1], "%)"), y = paste0("PC2 (", pct[2], "%)"))
)

theme_morpho <- theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom")

clade_layers <- list(
  scale_colour_manual(values = clade_cols, name = "Clade", drop = FALSE),
  scale_shape_manual(values = clade_shapes, name = "Clade", drop = FALSE),
  morpho_scales,
  theme_classic(base_size = 12),
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
        legend.position = "right",
        legend.text = element_text(size = 9),
        legend.key.height = grid::unit(0.4, "cm"),
        axis.text = element_text(colour = "black"))
)

# ---- clade morphospace ----
clade_df <- hs_df %>%
  filter(!clade %in% drop_clades) %>%
  mutate(clade = factor(as.character(clade), levels = clade_names))

p_morpho_clade <- ggplot() +
  geom_point(data = clade_df, aes(Axis1, Axis2, colour = clade, shape = clade),
             size = 2.5, alpha = 0.85) +
  clade_layers

# ---- clade phylomorphospace (same plot with the tree underneath) ----
clade_df_tree <- clade_df %>% filter(species %in% tree_pruned$tip.label)

p_phylo_clade <- ggplot() +
  geom_segment(data = phylo_df,
               aes(x = xstart, y = ystart, xend = xstop, yend = ystop),
               linewidth = 0.25, colour = "grey75") +
  geom_point(data = clade_df_tree, aes(Axis1, Axis2, colour = clade, shape = clade),
             size = 2.5, alpha = 0.85) +
  clade_layers

#ggsave("5.figuras/pcahillsmith_morphospace.pdf", p_morpho_clade, width = 8, height = 10)
#ggsave("5.figuras/clades_phylomorphospace.pdf",  p_phylo_clade,  width = 8, height = 10)

# ---- biplot: variable/category loadings over the species cloud ----
co_df <- data.frame(hs$co[, 1:4])
colnames(co_df) <- paste0("Axis", 1:4)
co_df$column   <- rownames(hs$co)                    # variable or category name
co_df$variable <- rownames(hs$cr)[hs$assign]         # parent variable
co_df$type     <- as.character(hs$index)[hs$assign]
co_df$len      <- sqrt(co_df$Axis1^2 + co_df$Axis2^2)

co_plot <- co_df %>% filter(len > quantile(len, 0.5))   # strongest 50% only; adjust as needed
mult    <- 0.8 * max(abs(hs$li[, 1:2])) / max(abs(co_plot[, c("Axis1", "Axis2")]))

p_biplot <- ggplot() +
  geom_point(data = hs_df, aes(Axis1, Axis2), colour = "grey70", size = 1.5) +
  geom_segment(data = co_plot,
               aes(x = 0, y = 0, xend = Axis1 * mult, yend = Axis2 * mult),
               arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.7) +
  geom_text_repel(data = co_plot,
                  aes(x = Axis1 * mult, y = Axis2 * mult, label = column),
                  size = 3, max.overlaps = Inf, show.legend = FALSE) +
  morpho_scales +
  theme_morpho

#ggsave("5.figuras/morphospace_contribution12.pdf", p_biplot, width = 8, height = 10, dpi = 300)

#============ caracteres categoricos no morfoespaco =============
####===========================================================##
### coloring categorical variables in the morphospace =====
## Final morphospace colored by variables state
##flower merosity 
merosity_plot <- ggplot(hs_df, aes(x = Axis1, y = Axis2, fill = factor(traits$flower_merosity))) +
  geom_point(shape = 21,size = 2.5, colour = "grey20", stroke = 0.3,alpha = 0.9) +
  scale_fill_manual(values = setNames(
    okabe[seq_along(unique(na.omit(traits$flower_merosity)))],
    sort(unique(na.omit(traits$flower_merosity)))),name = "Número de partes do perianto") +
  labs(x = paste0("PC1 (", pct[1], "%)"), y = paste0("PC2 (", pct[2], "%)")) +
  coord_equal() + scale_x_continuous(
    breaks = breaks,
    limits = lim,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = breaks,
    limits = lim,
    expand = c(0, 0)
  ) +
  theme_classic() + theme_bw() + theme(panel.grid = element_blank(), legend.position = "bottom")

ggsave("5.figuras/morphospace_with_merosity.pdf", merosity_plot, width=8, height=10)

##anther gland
anther_gland_plot <- ggplot(hs_df, aes(x = Axis1,y = Axis2,fill = factor(traits$anther_gland_presence))) +
  geom_point(shape = 21, size = 2.5,colour = "grey20", stroke = 0.3, alpha = 0.9) + 
  scale_fill_manual(
    values = setNames(
      okabe[seq_along(unique(na.omit(traits$anther_gland_presence)))],
      sort(unique(na.omit(traits$anther_gland_presence)))),name = "glândula nas anteras") + 
  labs(x = paste0("PC1 (", pct[1], "%)"), y = paste0("PC2 (", pct[2], "%)")) + coord_equal() +
  scale_x_continuous(
    breaks = breaks,
    limits = lim,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = breaks,
    limits = lim,
    expand = c(0, 0)
  ) +
  theme_classic() + theme_bw() + theme(panel.grid = element_blank(), legend.position = "bottom")

ggsave("5.figuras/morphospace_anther_gland.pdf", anther_gland_plot, width=8, height=10)

#nectary
nectary_plot <- ggplot(hs_df, aes(
  x = Axis1,
  y = Axis2,
  fill = factor(traits$nectary))) + geom_point(shape = 21, size = 2.5, colour = "grey20", stroke = 0.3, alpha = 0.9) +
  scale_fill_manual(
    values = setNames(
      okabe[seq_along(unique(na.omit(traits$nectary)))],
      sort(unique(na.omit(traits$nectary)))),
    name = "Nectário") + labs(
      x = paste0("PC1 (", pct[1], "%)"),
      y = paste0("PC2 (", pct[2], "%)")
    ) + coord_equal() +
  scale_x_continuous(
    breaks = breaks,
    limits = lim,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = breaks,
    limits = lim,
    expand = c(0, 0)
  ) + theme_classic() + theme_bw() + theme(panel.grid = element_blank(), legend.position = "bottom")

ggsave("5.figuras/morphospace_nectary_plot.pdf", nectary_plot, width=8, height=10)

#filament union
filament_union_plot <- ggplot(hs_df, aes(
  x = Axis1,
  y = Axis2,
  fill = factor(traits$filament_union))) +
  geom_point(shape = 21, size = 2.5, colour = "grey20", stroke = 0.3, alpha = 0.9) +
  scale_fill_manual(
    values = setNames(
      okabe[seq_along(unique(na.omit(traits$filament_union)))],
      sort(unique(na.omit(traits$filament_union)))
    ), name = "Filament union") +
  labs(x = paste0("PC1 (", pct[1], "%)"), y = paste0("PC2 (", pct[2], "%)")) +
  coord_equal() + scale_x_continuous(
    breaks = breaks,
    limits = lim,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = breaks,
    limits = lim,
    expand = c(0, 0)
  ) +
  theme_classic() + theme_bw() + theme(panel.grid = element_blank(), legend.position = "bottom")

ggsave("5.figuras/morphospace_filament_union_plot.pdf", filament_union_plot, width=8, height=10)

#aqui na cor dos filetes, pintar a bolinha da cor deles :)
# a coluna `color` virou 0-3 no script 01
#   "0" = "bicolored" "1" = "pink", 2 = white, 3 = yellow
pals <- list(
  color = c("0" = "#1D7634", "1" = "#F4A6C8", "2"= "#C1CDCD", "3" = "#F2D023"))

filament_color_plot <- ggplot(hs_df, aes(x = Axis1, y = Axis2, fill = factor(traits$color))) +
  geom_point(shape = 21, size = 2, colour = "grey", stroke = 0.2, alpha = 0.85) +
  scale_fill_manual(values = pals$color, name = "cor do filete") +
  labs(x = paste0("PC1 (", pct[1], "%)"),
       y = paste0("PC2 (", pct[2], "%)")) +
  coord_equal() + scale_x_continuous(
    breaks = breaks,
    limits = lim,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = breaks,
    limits = lim,
    expand = c(0, 0)
  ) +
  theme_classic() + theme_bw() + theme(panel.grid = element_blank(), legend.position = "bottom")


ggsave("5.figuras/morphospace_filament_color.pdf", filament_color, width=8, height=10)


# ====== SV, SR, MPD + summaries + tests =============
# Metrics:
#   SV  = sum of variances (spread around the centroid, per axis, summed)
#   SR  = sum of ranges    (extent occupied, per axis, summed; outlier-sensitive)
#   MPD = mean pairwise distance between species of the group
library(dispRity)
#bhatt coef fo probability of overlap between morphological disparity calculated for different groups
##SV for clade
metrics     <- list(SV = c(sum, variances),
                    SR = c(sum, ranges))
groups_list <- list(clade = groups_clade, biome = groups_biome)

# Results end up in disp_all[[group type]][[metric]]$disp / $smry / $wilcox / $bhatt
disp_all <- list()

for (grp in names(groups_list)) {
  
  subsets <- custom.subsets(scores_hs, groups_list[[grp]])
  boot    <- boot.matrix(subsets, bootstraps = bootstrap, rarefaction = rarefaction)
  
  for (m in names(metrics)) {
    d <- dispRity(boot, metric = metrics[[m]])
    disp_all[[grp]][[m]] <- list(
      disp   = d,
      smry   = summary(d),
      wilcox = test.dispRity(d, wilcox.test, correction = "holm"),
      bhatt  = test.dispRity(d, bhatt.coeff)
    )
  }
}

#check results
str(disp_all$clade$SV$bhatt)
str(disp_all$clade$SV$wilcox)

str(disp_all$biome$SV$smry)

##===== MPD ================================================
# Mean pairwise distance among species of a group (mean of the lower triangle
# of the Gower matrix). Needs >= 2 species.

sub_labels <- c(clade = "Clades", biome = "Biomes")
mpd_rows   <- list()

for (grp in names(groups_list)) {
  for (g in names(groups_list[[grp]])) {
    sp <- intersect(groups_list[[grp]][[g]], rownames(gower))
    d  <- gower[sp, sp, drop = FALSE]
    mpd_rows[[paste(grp, g)]] <- data.frame(
      Subdivision = sub_labels[[grp]],
      Group       = g,
      Sample_size = length(sp),
      MPD         = if (length(sp) >= 2) round(mean(d[upper.tri(d)], na.rm = TRUE), 4) else NA_real_
    )
  }
}

mpd_all <- bind_rows(mpd_rows)
mpd_all

#check some significance for mpd values
#conferir se o mpd é maior ou menos do que o esperado ao acaso!!!
#SES mean standardized effect size

n_null <- 1000
pool   <- setdiff(rownames(scores_hs), meta$species[meta$clade %in% drop_clades])  # todas as espécies, sem outgroup

mpd_all$MPD_null <- NA_real_
mpd_all$SES      <- NA_real_
mpd_all$p_value  <- NA_real_

for (i in seq_len(nrow(mpd_all))) {
  n <- mpd_all$Sample_size[i]
  if (n < 2) next
  
  null_mpd <- numeric(n_null)
  for (r in seq_len(n_null)) {
    sp <- sample(pool, n)                     # grupo aleatório do mesmo tamanho
    d  <- gower[sp, sp]
    null_mpd[r] <- mean(d[upper.tri(d)], na.rm = TRUE)
  }
  
  mpd_all$MPD_null[i] <- round(mean(null_mpd), 4)
  mpd_all$SES[i]      <- round((mpd_all$MPD[i] - mean(null_mpd)) / sd(null_mpd), 2)
  mpd_all$p_value[i]  <- min(1, 2 * min(mean(null_mpd <= mpd_all$MPD[i]),
                                        mean(null_mpd >= mpd_all$MPD[i])))
}
mpd_all$p_adj <- round(p.adjust(mpd_all$p_value, method = "holm"), 3)

### plot violin for 
##Prosopis and samanea has less than 3 elemnts
#adapted from casali et al 
str(disp_all$clade$SV$bhatt)


disp_all$clade$SV$disp$disparity$`Adenanthera clade`[[3]]
sv_clade$disparity$`Adenanthera clade` #valores de sv do bootstrap
sv_clade$disparity$`Albizia clade` #valores de sv do bootstrap
# exemplos do violin pra mimosa e adenanthera
sv_clade$disparity$`Adenanthera clade`[[3]] #dados rarefeiros para usar no violin plot
sv_clade$disparity$`Cojoba clade`[[2]] #dados observados pois n é igual ao numero de rarefacao

#CLADE - sum of variances ===============
####===========================================================##
sum_1 <- as.vector(disp_all$clade$SV$disp$disparity$`Adenanthera clade`[[3]])
C1 <- c(rep("Adenanthera", bootstrap))
sv_1 <- data.frame(C1, sum_1)
colnames(sv_1) <- c("Group", "Sum_of_variance")

sum_2 <- as.vector(disp_all$clade$SV$disp$disparity$`Albizia clade`[[3]])
C2 <- c(rep("Albizia", bootstrap))
sv_2 <- data.frame(C2, sum_2)
colnames(sv_2) <- c("Group", "Sum_of_variance")

sum_3 <- as.vector(disp_all$clade$SV$disp$disparity$`Archidendron clade`[[3]])
C3 <- c(rep("Archidendron", bootstrap))
sv_3 <- data.frame(C3, sum_3)
colnames(sv_3) <- c("Group", "Sum_of_variance")

sum_4 <- as.vector(disp_all$clade$SV$disp$disparity$`Calliandra clade`[[3]])
C4 <- c(rep("Calliandra", bootstrap))
sv_4 <- data.frame(C4, sum_4)
colnames(sv_4) <- c("Group", "Sum_of_variance")

sum_5 <- as.vector(disp_all$clade$SV$disp$disparity$`Cojoba clade`[[2]]) #use the observed data
C5 <- c(rep("Cojoba", bootstrap))
sv_5 <- data.frame(C5, sum_5)
colnames(sv_5) <- c("Group", "Sum_of_variance")

sum_6 <- as.vector(disp_all$clade$SV$disp$disparity$`Dichrostachys clade`[[3]])
C6 <- c(rep("Dichrostachys", bootstrap))
sv_6 <- data.frame(C6, sum_6)
colnames(sv_6) <- c("Group", "Sum_of_variance")

sum_7 <- as.vector(disp_all$clade$SV$disp$disparity$`Entada clade`[[3]])
C7 <- c(rep("Entada", bootstrap))
sv_7 <- data.frame(C7, sum_7)
colnames(sv_7) <- c("Group", "Sum_of_variance")

sum_8 <- as.vector(disp_all$clade$SV$disp$disparity$`Inga clade`[[3]])
C8 <- c(rep("Inga", bootstrap))
sv_8 <- data.frame(C8, sum_8)
colnames(sv_8) <- c("Group", "Sum_of_variance")

sum_9 <- as.vector(disp_all$clade$SV$disp$disparity$`Jupunba clade`[[3]])
C9 <- c(rep("Jupunba", bootstrap))
sv_9 <- data.frame(C9, sum_9)
colnames(sv_9) <- c("Group", "Sum_of_variance")

sum_10 <- as.vector(disp_all$clade$SV$disp$disparity$`Mimosa clade`[[3]])
C10 <- c(rep("Mimosa", bootstrap))
sv_10 <- data.frame(C10, sum_10)
colnames(sv_10) <- c("Group", "Sum_of_variance")

sum_11 <- as.vector(disp_all$clade$SV$disp$disparity$`Neltuma clade`[[3]])
C11 <- c(rep("Neltuma", bootstrap))
sv_11 <- data.frame(C11, sum_11)
colnames(sv_11) <- c("Group", "Sum_of_variance")

sum_12 <- as.vector(disp_all$clade$SV$disp$disparity$`Newtonia grade`[[3]])
C12 <- c(rep("Newtonia", bootstrap))
sv_12 <- data.frame(C12, sum_12)
colnames(sv_12) <- c("Group", "Sum_of_variance")

sum_13 <- as.vector(disp_all$clade$SV$disp$disparity$`Parkia clade`[[3]])
C13 <- c(rep("Parkia", bootstrap))
sv_13 <- data.frame(C13, sum_13)
colnames(sv_13) <- c("Group", "Sum_of_variance")

sum_14 <- as.vector(disp_all$clade$SV$disp$disparity$`Pithecellobium clade`[[3]])
C14 <- c(rep("Pithecellobium", bootstrap))
sv_14 <- data.frame(C14, sum_14)
colnames(sv_14) <- c("Group", "Sum_of_variance")

sum_15 <- as.vector(disp_all$clade$SV$disp$disparity$`Prosopis clade`[[2]])
C15 <- c(rep("Prosopis", bootstrap))
sv_15 <- data.frame(C15, sum_15)
colnames(sv_15) <- c("Group", "Sum_of_variance")

sum_16 <- as.vector(disp_all$clade$SV$disp$disparity$`Samanea clade`[[2]]) #use the observed data
C16 <- c(rep("Samanea", bootstrap))
sv_16 <- data.frame(C16, sum_16)
colnames(sv_16) <- c("Group", "Sum_of_variance")

sum_17 <- as.vector(disp_all$clade$SV$disp$disparity$`Senegalia grade`[[3]])
C17 <- c(rep("Senegalia", bootstrap))
sv_17 <- data.frame(C17, sum_17)
colnames(sv_17) <- c("Group", "Sum_of_variance")

sum_18 <- as.vector(disp_all$clade$SV$disp$disparity$`Stryphnodendron clade`[[3]])
C18 <- c(rep("Stryphnodendron", bootstrap))
sv_18 <- data.frame(C18, sum_18)
colnames(sv_18) <- c("Group", "Sum_of_variance")

sum_19 <- as.vector(disp_all$clade$SV$disp$disparity$`Zapoteca clade`[[3]])
C19 <- c(rep("Zapoteca", bootstrap))
sv_19 <- data.frame(C19, sum_19)
colnames(sv_19) <- c("Group", "Sum_of_variance")

clades_to_plot <- rbind(sv_1, sv_2, sv_3, sv_4, sv_5, sv_6, sv_7, sv_8, sv_9, sv_10,
                        sv_11, sv_12, sv_13, sv_14, sv_15, sv_16, sv_17, sv_18, sv_19)
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
                                         "Prosopis" ,"Samanea","Senegalia","Stryphnodendron","Zapoteca"))


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

#ggsave("5.figuras/clade_sv_plot.pdf", clade_sv_violin,width=8, height=10)


#BIOMES - sum of variances =====================================
####===========================================================##
disp_all$biome$SV$disp$disparity$`Tropical & Subtropical Moist Broadleaf Forests`[[3]]

sum_1 <- as.vector(disp_all$biome$SV$disp$disparity$`Tropical & Subtropical Moist Broadleaf Forests`[[3]])
C1 <- c(rep("1", length(sum_1)))
sv_1 <- data.frame(C1, sum_1)
colnames(sv_1) <- c("Group", "Sum_of_variance")

sum_2 <- as.vector(disp_all$biome$SV$disp$disparity$`Tropical & Subtropical Dry Broadleaf Forests`[[3]])
C2 <- c(rep("2", length(sum_2)))
sv_2 <- data.frame(C2, sum_2)
colnames(sv_2) <- c("Group", "Sum_of_variance")

sum_3 <- as.vector(disp_all$biome$SV$disp$disparity$`Tropical & Subtropical Coniferous Forests`[[3]]) #use [[2]] if n = rarefaction value
C3 <- c(rep("3", length(sum_3)))
sv_3 <- data.frame(C3, sum_3)
colnames(sv_3) <- c("Group", "Sum_of_variance")

sum_4 <- as.vector(disp_all$biome$SV$disp$disparity$`Temperate Broadleaf & Mixed Forests`[[3]])
C4 <- c(rep("4", length(sum_4)))
sv_4 <- data.frame(C4, sum_4)
colnames(sv_4) <- c("Group", "Sum_of_variance")

sum_5 <- as.vector(disp_all$biome$SV$disp$disparity$`Temperate Conifer Forests`[[2]])
C5 <- c(rep("5", length(sum_5)))
sv_5 <- data.frame(C5, sum_5)
colnames(sv_5) <- c("Group", "Sum_of_variance")

sum_7 <- as.vector(disp_all$biome$SV$disp$disparity$`Tropical & Subtropical Grasslands, Savannas & Shrublands`[[3]])
C7 <- c(rep("7", length(sum_7)))
sv_7 <- data.frame(C7, sum_7)
colnames(sv_7) <- c("Group", "Sum_of_variance")

sum_8 <- as.vector(disp_all$biome$SV$disp$disparity$`Temperate Grasslands, Savannas & Shrublands`[[3]])
C8 <- c(rep("8", length(sum_8)))
sv_8 <- data.frame(C8, sum_8)
colnames(sv_8) <- c("Group", "Sum_of_variance")

sum_9 <- as.vector(disp_all$biome$SV$disp$disparity$`Flooded Grasslands & Savannas`[[3]])
C9 <- c(rep("9", length(sum_9)))
sv_9 <- data.frame(C9, sum_9)
colnames(sv_9) <- c("Group", "Sum_of_variance")

sum_10 <- as.vector(disp_all$biome$SV$disp$disparity$`Montane Grasslands & Shrublands`[[3]])
C10 <- c(rep("10", length(sum_10)))
sv_10 <- data.frame(C10, sum_10)
colnames(sv_10) <- c("Group", "Sum_of_variance")

sum_12 <- as.vector(disp_all$biome$SV$disp$disparity$`Mediterranean Forests, Woodlands & Scrub`[[3]])
C12 <- c(rep("12", length(sum_12)))
sv_12 <- data.frame(C12, sum_12)
colnames(sv_12) <- c("Group", "Sum_of_variance")

sum_13 <- as.vector(disp_all$biome$SV$disp$disparity$`Deserts & Xeric Shrublands`[[3]])
C13 <- c(rep("13", length(sum_13)))
sv_13 <- data.frame(C13, sum_13)
colnames(sv_13) <- c("Group", "Sum_of_variance")

sum_14 <- as.vector(disp_all$biome$SV$disp$disparity$Mangroves[[3]])
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
#ggsave("5.figuras/biomes_sv_plot.pdf", biome_sv_violin,width=8, height=10)

# ---- Bhattacharyya coefficient heatmaps (SV) ----
### heatmap for clades =====
# library(pheatmap)
# library(ggplot2)
# library(dplyr)

# Extract the clade pairs and coefficients
heatmap_data <- disp_all$clade$SV$bhatt %>%
  mutate(
    clade1 = sub(":.*", "", rownames(disp_all$clade$SV$bhatt)),
    clade2 = sub(".*:", "", rownames(disp_all$clade$SV$bhatt))
  )

ggplot(heatmap_data, aes(x = clade2, y = clade1, fill = bhatt.coeff)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient(
    low = "white",
    high = "red",
    name = "Bhatt coefficient"
  ) +
  theme_minimal() +
  labs(
    x = "Clade",
    y = "Clade",
    title = "Bhatt coefficient heatmap"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  coord_fixed()



##### plot heatmap for biomes ========
### heatmap for biomes =====
heatmap_data_biomes <- disp_all$biome$SV$bhatt %>%
  mutate(biome1 = sub(":.*", "", rownames(disp_all$biome$SV$bhatt)),
         biome2 = sub(".*:", "", rownames(disp_all$biome$SV$bhatt)))

p_heat_biomes <- ggplot(heatmap_data_biomes, aes(x = biome2, y = biome1, fill = bhatt.coeff)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient(
    low = "white",
    high = "red",
    name = "Bhatt coefficient"
  ) +
  theme_minimal() +
  labs(
    x = "biome",
    y = "biome",
    title = "Bhatt coefficient heatmap"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  coord_fixed()

ggsave("5.figuras/heatmap_bhatt_biome.pdf", p_heat_biomes, width = 12, height = 10)


# Save results so the bootstrap does not need to be re-run
saveRDS(list(disp_all = disp_all, mpd = mpd_all, n_axes = n_axes, eig = eig_table),
        "3.outputs/disparity_results.rds")


#save results 
as.data.frame(disp_all$clade$SV$smry)

 