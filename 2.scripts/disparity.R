## DISPARITY ANALYSES ======
## calcular distancia de Gower -> PCoA -> disparidade (SV, SR, MPD) -> plots
setwd("~/Documents/Labis/Dados/") #defining work directory
if (!require(librarian)) install.packages("librarian"); library("librarian")
librarian::shelf(phytools, dplyr, tidyr, purrr, factoextra, vegan, tidyverse,ape, stringr, readr,
                 tidytree, ggtree, cluster, dispRity) #installing and/or loading packages

##=== Limpeza e organização dos dados =======
##testing with a small dataset of mimosa
dat_flower <- read.csv("4.outputs/morpho_10perc_genus.csv", )
dat_flower <- dat_flower %>%
  filter(str_detect(taxon, "Mimosa"))

#we have only 62 species from Mimosa here 

tree <- read.tree("3.trees/mimosoid_calibrated_clean_updated.tre")
dat_eco <- #extract from ringelberg?
  
  ## get traits dataframe
  traits <- cbind("taxon" = dat_flower$taxon, dat_flower[, 4:33])

set.seed(7)
## sinonimos: staminal_tube == stamen_tube
traits <- traits %>%
  mutate(stamen_tube_final = coalesce(stamen_tube_presence, staminal_tube_presence)) %>%
  select(-stamen_tube_presence, -staminal_tube_presence)

## prune phylogeny
tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, 
                                      dat_flower$taxon))

##inflorescence shape, calyx shape, 
#inflorescence width, pedicel width, 
#calyx lobe length, corolla lobe length serão todos excluídos. 
#dados de sexualidade das flores precisam ser pensados como tratar
traits_mixed <- traits %>%
  column_to_rownames("taxon") %>%
  select(where(~ !all(is.na(.))))   # remove colunas 100% NA (ex: inflorescence_width_mean)

traits_mixed <- traits_mixed %>%
  mutate(
    flower_merosity = factor(
      flower_merosity,
      levels = c("3-merous", "4-merous", "5-merous", "variable"),
      ordered = FALSE #pra mim nao faz sentido ordenar
    )
  ) 
###mudar esse script aqui pra selecionar 
#salvar um csv do unique dos traços pra fazer uma validação dos dados

##tratamento dos dados validados

## stamen_count 
traits_mixed <- traits_mixed %>%
  mutate(
    stamen_count_state = case_when(
      stamen_count == 4  ~ "haplostemonous",
      stamen_count %in% c(8, 10) ~ "diplostemonous",
      stamen_count > 10  ~ "polystemonous",
      TRUE ~ NA_character_
    ),
    stamen_count_state = factor(
      stamen_count_state,
      levels = c("haplostemonous", "diplostemonous", "polystemonous"),
      ordered = FALSE
    )
  ) %>%
  select(-stamen_count)   # evita duplicar a mesma informacao

## as outras colunas vou transformar em factor 
traits_mixed <- traits_mixed %>%
  mutate(across(where(is.character), as.factor))

#vou selecionar umas variaveis pra testar as analises de disparidade 
vars_selecionadas <- c(
  "inflorescence_type",
  "flower_merosity",
  "filament_color",
  "inflorescence_length_mean",
  "inflorescence_peduncle_length_mean",
  "height_mean",
  "calyx_length_mean",
  "stamen_count_state",
  "filament_length_mean",
  "corolla_lobe_width_mean",
  "corolla_length_mean")

## confere se algum nome foi digitado errado antes de filtrar
setdiff(vars_selecionadas, colnames(traits_mixed))

traits_mixed <- traits_mixed %>%
  select(all_of(vars_selecionadas))

## conferencia
sapply(traits_mixed, class)

## Input de dados com RPhylopars ============ 
install.packages("Rphylopars")
#downloades Rphylopars_0.3.10
library(Rphylopars)

#Rphylopars so funciona para variaveis CONTINUAS (assume um modelo de
## evolucao do caracter - BM ou OU - correlacionado com a filogenia)
cont_cols <- traits_mixed %>%
  select(where(is.numeric)) %>%
  colnames()

## Teste de normalidade nos dados iniciais ======
cont_raw <- traits_mixed %>%
  rownames_to_column("species") %>%
  select(species, all_of(cont_cols))

shapiro_results <- lapply(cont_raw[cont_cols], function(x) {
  x <- na.omit(x)
  if (length(x) >= 3) shapiro.test(x) else NA
})

data.frame(
  trait   = names(shapiro_results),
  W       = sapply(shapiro_results, function(r) if (is.list(r)) r$statistic else NA),
  p.value = sapply(shapiro_results, function(r) if (is.list(r)) r$p.value else NA)
)

## log no dado bruto, antes de imputar (NAs continuam ali, log(NA) = NA,
##    entao nao precisa remover NA antes)
cont_log <- cont_raw
cont_log[cont_cols] <- lapply(cont_log[cont_cols], log)

## Checagem de alinhamento entre taxons e arvore antes do phylopars =====
setdiff(cont_log$species, tree_pruned$tip.label)
setdiff(tree_pruned$tip.label, cont_log$species)

#o setdiff precisa vir (0)

## phylopars() agora roda na escala log ======
phylopars_fit <- phylopars(
  trait_data       = cont_log,
  tree             = tree_pruned,
  model            = "BM",
  pheno_error      = TRUE,
  phylo_correlated = TRUE,
  pheno_correlated = TRUE
)

## Extrair valores das pontas ==== 
n_tip <- length(tree_pruned$tip.label)
imputed_cont <- phylopars_fit$anc_recon[1:n_tip, , drop = FALSE]

## checagem defensiva: a ordem das pontas em anc_recon bate com a arvore?
name_order <- identical(rownames(imputed_cont), tree_pruned$tip.label)
name_order 

if (!ordem_ok) {
  imputed_cont <- imputed_cont[match(tree_pruned$tip.label, rownames(imputed_cont)), ]
  stopifnot(identical(rownames(imputed_cont), tree_pruned$tip.label))
}

## agora sim, se necessario, garanta o rotulo (so depois de confirmar a ordem)
rownames(imputed_cont) <- tree_pruned$tip.label

## quantos NA existiam antes vs depois (em escala log)
sum(is.na(cont_log[cont_cols]))
sum(is.na(imputed_cont))   # deve ser 0


## Substituir em traits_mixed - agora traits_mixed fica DIRETO em escala log
imputed_df <- as.data.frame(imputed_cont) %>%
  rownames_to_column("species")

traits_mixed <- traits_mixed %>%
  rownames_to_column("species") %>%
  select(-all_of(cont_cols)) %>%
  left_join(imputed_df, by = "species") %>%
  column_to_rownames("species")

sapply(traits_mixed, class)
sum(sapply(traits_mixed, function(x) sum(is.na(x))))   # so as categoricas devem ter NA

## NAO tem mais etapa de "log nas continuas" depois disso - ja esta feito

## DISTANCIA DE GOWER 
## type = list(...) permite controlar como cada bloco de variaveis e tratado
## symm = fatores nominais tratados como "simple matching" (0/1)
library(cluster)
gower_d <- daisy(traits_mixed, metric = "gower") ##cluster package

## checagem: quantos pares tem NA na distancia (caso alguma linha nao compartilhe
## nenhuma variavel observada com outra - daria distancia NA)
sum(is.na(as.matrix(gower_d)))

## PCoA com a matriz de gower =====================================
pcoa_res <- pcoa(gower_d).  #ape package

## % de variancia explicada por eixo
pcoa_res$values$Relative_eig[1:5]

scores_pcoa <- pcoa_res$vectors[, 1:2]
colnames(scores_pcoa) <- c("PCo1", "PCo2")

## garantir ordem identica a arvore 
scores_pcoa <- scores_pcoa[match(tree_pruned$tip.label, rownames(scores_pcoa)), ]
stopifnot(identical(rownames(scores_pcoa), tree_pruned$tip.label))
## ============================================ #

## SV, SR, MPD ======================================
## --- SV (sum of variances) e SR (sum of ranges) via dispRity,
##     calculadas sobre TODOS os eixos retidos do PCoA (nao so os 2 primeiros,
##     para nao perder disparidade que esta em eixos de ordem maior) ---

disparity_sv <- dispRity(pcoa_res$vectors, metric = c(sum, variances))
disparity_sr <- dispRity(pcoa_res$vectors, metric = c(sum, ranges))

summary(disparity_sv)   # Sum of Variance
summary(disparity_sr)   # Sum of Range

## --- MPD (mean pairwise distance) direto da matriz de Gower (nao do PCoA) ---
gower_mat <- as.matrix(gower_d)
mpd <- mean(gower_mat[lower.tri(gower_mat)], na.rm = TRUE)
mpd

## aqui posso incluir selecao em grupo pra definir disparidade
#por clado ou por grupos ecologicos,por exemplo
# grupos <- traits_mixed$clade
# names(grupos) <- rownames(traits_mixed)

## dispRity aceita uma lista de subconjuntos de taxons por grupo
## clados <- split(names(grupos), grupos)

# disparity_por_grupo <- dispRity(
#  all_axes,
#  metric = c(sum, variances),
#  group  = grupos_lista)
# summary(disparity_por_grupo)

## morphospace (PCoA simples, sem arvore)
df_plot <- data.frame(
  taxon = rownames(scores_pcoa),
  PCo1  = scores_pcoa[, "PCo1"],
  PCo2  = scores_pcoa[, "PCo2"])

var_pco1 <- round(pcoa_res$values$Relative_eig[1] * 100, 1)
var_pco2 <- round(pcoa_res$values$Relative_eig[2] * 100, 1)

p_morphospace <- ggplot(df_plot, aes(x = PCo1, y = PCo2)) +
  geom_point(size = 3, alpha = 0.8, color = "steelblue") +
  labs(
    title = "Morfoespaco floral - Mimosa (PCoA)",
    x = paste0("PCo1 (", var_pco1, "%)"),
    y = paste0("PCo2 (", var_pco2, "%)")
  ) +
  theme_minimal()

print(p_morphospace)

ggsave("5.figures/morphospace_mimosa.png", p_morphospace, width = 7, height = 5, dpi = 300)

## phylomorphospace ============================================================

## phylomorphospace exige matrix numerica pura, na mesma ordem da arvore
scores_mat <- as.matrix(scores_pcoa)

phylomorphospace(
  tree_pruned,
  scores_mat,
  xlab = paste0("PCo1 (", var_pco1, "%)"),
  ylab = paste0("PCo2 (", var_pco2, "%)"),
  label = "off"
)

#pra salvar o plot 

png("5.figures/phylomorphospace_mimosa.png", width = 1600, height = 1200, res = 200)
phylomorphospace(
  tree_pruned,
  scores_mat,
  xlab = paste0("PCo1 (", var_pco1, "%)"),
  ylab = paste0("PCo2 (", var_pco2, "%)"),
  label = "off",
  node.size = c(0, 1.2)
)
title("Phylomorphospace Mimosa")
dev.off()

