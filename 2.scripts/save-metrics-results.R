###make table of all results
##run this after pca-hill-smith script

library(dplyr)
library(tidyr)
library(tibble)
library(flextable)
library(officer)

## SR - Sum of ranges

clade_sr_table <- clade_sr_smry %>%
  mutate(Subdivision = "Clades") %>%
  transmute(
    Subdivision,
    Group = subsets,
    `Sample size` = n,
    `SR (Obs.)` = round(obs, 2),
    `Median (BS)` = round(bs.median, 2),
    `CI (BS)` = paste0(
      round(`2.5%`, 2), "–",
      round(`97.5%`, 2)
    )
  )

biomes_sr_table <- biomes_sr_smry %>%
  mutate(Subdivision = "Biomes") %>%
  transmute(
    Subdivision,
    Group = subsets,
    `Sample size` = n,
    `SR (Obs.)` = round(obs, 2),
    `Median (BS)` = round(bs.median, 2),
    `CI (BS)` = paste0(
      round(`2.5%`, 2), "–",
      round(`97.5%`, 2)
    )
  )

all_sr_summary <- bind_rows(
  clade_sr_table,
  biomes_sr_table
)


## SV - Sum of variances
clade_sv_table <- clade_sv_smry %>%
  mutate(Subdivision = "Clades") %>%
  transmute(
    Subdivision,
    Group = subsets,
    `Sample size` = n,
    `SV (Obs.)` = round(obs, 2),
    `Median (BS)` = round(bs.median, 2),
    `CI (BS)` = paste0(
      round(`2.5%`, 2), "–",
      round(`97.5%`, 2)
    )
  )

biomes_sv_table <- biomes_sv_smry %>%
  mutate(Subdivision = "Biomes") %>%
  transmute(
    Subdivision,
    Group = subsets,
    `Sample size` = n,
    `SV (Obs.)` = round(obs, 2),
    `Median (BS)` = round(bs.median, 2),
    `CI (BS)` = paste0(
      round(`2.5%`, 2), "–",
      round(`97.5%`, 2)
    )
  )

all_sv_summary <- bind_rows(
  clade_sv_table,
  biomes_sv_table
)


## SV - Clades
bc_sv_clade <- clade_sv_bhatt %>%
  rownames_to_column("Groups") %>%
  mutate(
    Subdivision = "Clades",
    BC = round(bhatt.coeff, 3)
  ) %>%
  select(Subdivision, Groups, BC)


## SV - Biomes
bc_sv_biomes <- biomes_sv_bhatt %>%
  rownames_to_column("Groups") %>%
  mutate(
    Subdivision = "Biomes",
    BC = round(bhatt.coeff, 3)
  ) %>%
  select(Subdivision, Groups, BC)


## Combine SV
bc_sv <- bind_rows(
  bc_sv_clade,
  bc_sv_biomes
)

## SR - Clades
bc_sr_clade <- clade_sr_bhatt %>%
  rownames_to_column("Groups") %>%
  mutate(
    Subdivision = "Clades",
    BC = round(bhatt.coeff, 3)
  ) %>%
  select(Subdivision, Groups, BC)


## SR - Biomes
bc_sr_biomes <- biomes_sr_bhatt %>%
  rownames_to_column("Groups") %>%
  mutate(
    Subdivision = "Biomes",
    BC = round(bhatt.coeff, 3)
  ) %>%
  select(Subdivision, Groups, BC)


## Combine SR
bc_sr <- bind_rows(
  bc_sr_clade,
  bc_sr_biomes
)

str(clade_sr_wilcox)

## SR - Clades
wilcox_sr_clade <- data.frame(
  Groups = rownames(clade_sr_wilcox[[2]]),
  Statistic_W = clade_sr_wilcox[[1]]$statistic,
  P_value = as.numeric(clade_sr_wilcox[[2]][, 1])
) %>%
  mutate(
    Subdivision = "Clades"
  ) %>%
  select(Subdivision, Groups, Statistic_W, P_value)


## SR - Biomes
wilcox_sr_biomes <- data.frame(
  Groups = rownames(biomes_sr_wilcox[[2]]),
  Statistic_W = biomes_sr_wilcox[[1]]$statistic,
  P_value = as.numeric(biomes_sr_wilcox[[2]][, 1])
) %>%
  mutate(
    Subdivision = "Biomes"
  ) %>%
  select(Subdivision, Groups, Statistic_W, P_value)


## Combine SR
wilcox_sr <- bind_rows(
  wilcox_sr_clade,
  wilcox_sr_biomes)

## SV - Clades
wilcox_sv_clade <- data.frame(
  Groups = rownames(clade_sv_wilcox[[2]]),
  Statistic_W = clade_sv_wilcox[[1]]$statistic,
  P_value = as.numeric(clade_sv_wilcox[[2]][, 1])
) %>%
  mutate(
    Subdivision = "Clades"
  ) %>%
  select(Subdivision, Groups, Statistic_W, P_value)


## SV - Biomes
wilcox_sv_biomes <- data.frame(
  Groups = rownames(biomes_sv_wilcox[[2]]),
  Statistic_W = biomes_sv_wilcox[[1]]$statistic,
  P_value = as.numeric(biomes_sv_wilcox[[2]][, 1])
) %>%
  mutate(
    Subdivision = "Biomes"
  ) %>%
  select(Subdivision, Groups, Statistic_W, P_value)


## Combine SV
wilcox_sv <- bind_rows(
  wilcox_sv_clade,
  wilcox_sv_biomes
)

## SV pairwise comparisons
table_sv_tests <- left_join(
  wilcox_sv,
  bc_sv,
  by = c("Subdivision", "Groups")
) %>%
  mutate(
    Statistic_W = round(Statistic_W, 2),
    P_value = signif(P_value, 3)
  ) %>%
  rename(
    `Statistic W` = Statistic_W,
    `P-value` = P_value
  )


## SR pairwise comparisons
table_sr_tests <- left_join(
  wilcox_sr,
  bc_sr,
  by = c("Subdivision", "Groups")
) %>%
  mutate(
    Statistic_W = round(Statistic_W, 2),
    P_value = signif(P_value, 3)
  ) %>%
  rename(
    `Statistic W` = Statistic_W,
    `P-value` = P_value
  )

## Summary tables
ft_sr_summary <- flextable(all_sr_summary) %>%
  theme_booktabs() %>%
  fontsize(size = 9, part = "all") %>%
  bold(part = "header") %>%
  align(align = "center", part = "all") %>%
  autofit()


ft_sv_summary <- flextable(all_sv_summary) %>%
  theme_booktabs() %>%
  fontsize(size = 9, part = "all") %>%
  bold(part = "header") %>%
  align(align = "center", part = "all") %>%
  autofit()


## Pairwise comparisons
ft_sv <- flextable(table_sv_tests) %>%
  theme_booktabs() %>%
  fontsize(size = 9, part = "all") %>%
  bold(part = "header") %>%
  align(align = "center", part = "all") %>%
  autofit()


ft_sr <- flextable(table_sr_tests) %>%
  theme_booktabs() %>%
  fontsize(size = 9, part = "all") %>%
  bold(part = "header") %>%
  align(align = "center", part = "all") %>%
  autofit()

##calculate mpd for each clade
## Calculate MPD for each clade
mpd_clade <- data.frame(
  Clade = names(groups_clade),
  Sample_size = NA_integer_,
  MPD = NA_real_
)

for (i in seq_along(groups_clade)) {
  
  sp <- groups_clade[[i]]
  
  d <- as.matrix(gower[sp, sp])
  
  mpd_clade$Sample_size[i] <- length(sp)
  
  mpd_clade$MPD[i] <- if (length(sp) >= 2) {
    mean(d[upper.tri(d)], na.rm = TRUE)
  } else {
    NA_real_
  }
}

mpd_clade$MPD <- round(mpd_clade$MPD, 4)

mpd_clade

##mpd for each biome 
## Filter biome groups

groups_biome_mpd <- lapply(
  groups_biome,
  intersect,
  y = rownames(gower_biomes)
)


## Calculate MPD

mpd_biomes <- data.frame(
  Biome = names(groups_biome),
  Sample_size = NA_integer_,
  MPD = NA_real_
)

for (i in seq_along(groups_biome)) {
  
  sp <- groups_biome[[i]]
  
  d <- as.matrix(gower_biomes[sp, sp])
  
  mpd_biomes$Sample_size[i] <- length(sp)
  
  mpd_biomes$MPD[i] <- if (length(sp) >= 2) {
    mean(d[upper.tri(d)], na.rm = TRUE)
  } else {
    NA_real_
  }
}

mpd_biomes$MPD <- round(mpd_biomes$MPD, 4)

mpd_biomes

## MPD
mpd_all <- bind_rows(
  mpd_clade %>%
    rename(Group = Clade) %>%
    mutate(Subdivision = "Clades"),
  
  mpd_biomes %>%
    rename(Group = Biome) %>%
    mutate(Subdivision = "Biomes")
) %>%
  select(Subdivision, Group, Sample_size, MPD)

ft_mpd <- flextable(mpd_all) %>%
  theme_booktabs() %>%
  fontsize(size = 9, part = "all") %>%
  bold(part = "header") %>%
  align(align = "center", part = "all") %>%
  autofit()



##save word document
doc <- read_docx()

## Table S1
doc <- body_add_par(
  doc,
  "Table S1. Morphospace occupancy (SR)",
  style = "heading 1"
)

doc <- body_add_flextable(doc, ft_sr_summary)


## Table S2
doc <- body_add_par(
  doc,
  "Table S2. Morphospace occupancy (SV)",
  style = "heading 1"
)

doc <- body_add_flextable(doc, ft_sv_summary)


## Table S3
doc <- body_add_par(
  doc,
  "Table S3. Pairwise comparisons of SR",
  style = "heading 1"
)

doc <- body_add_flextable(doc, ft_sr)


## Table S4
doc <- body_add_par(
  doc,
  "Table S4. Pairwise comparisons of SV",
  style = "heading 1"
)

doc <- body_add_flextable(doc, ft_sv)


## Table S5
doc <- body_add_par(
  doc,
  "Table S5. Mean pairwise morphological distances",
  style = "heading 1"
)

doc <- body_add_flextable(doc, ft_mpd)


## Export
print(
  doc,
  target = "3.outputs/Disparity_results.docx"
)
