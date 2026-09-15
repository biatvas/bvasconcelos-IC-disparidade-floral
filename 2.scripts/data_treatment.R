if (!require(librarian)) install.packages("librarian")
librarian::shelf(dplyr, purrr, readr, stringr, tidyr, tibble,
                 cluster, ape, vegan, ggplot2, readr, ade4, FactoMineR, 
                 tibble, stats)

#limpeza e normalização dos dados iniciais
morpho_data<- read.csv("~/Documents/GitHub/bvasconcelos-IC-disparidade-floral/1.datasets/mimoseae_subset_clean.csv")
validated_data <- morpho_data %>% filter(Check == "1") #221 obs

traits <- cbind("taxon" = validated_data$taxon, validated_data[, 6:83])

#agora vou trabalhar com o arquivo traits pra limpeza e gerar unique csv dos dados 
#a etapa consiste em 
#aplicar um unique csv e gerar arquivo de substituição dos dados qualitativos +
#transformar dados continuos pra mesma escala de medida e calcular média
# depois dessas etapas, conferir NAs e realizar imputação e também calcular Gower sem input
setwd("Documents/GitHub/bvasconcelos-IC-disparidade-floral/")

get_range_traits <- function(cols) {
  unique(sub("(.+)_min(_|$).*$", "\\1", cols[grepl("_min(_|$)", cols)]))}

get_range_cols <- function(cols, range_traits, suffix) {
  unlist(lapply(range_traits, function(x)
    grep(paste0("^", x, "_", suffix, "(_|$)"), cols, value = TRUE)))
}

cols <- names(traits)
range_traits <- get_range_traits(cols)

continuous_col <- c(
  get_range_cols(cols, range_traits, "min"),
  get_range_cols(cols, range_traits, "low"),
  get_range_cols(cols, range_traits, "high"),
  get_range_cols(cols, range_traits, "max")
)
unit_col <- grep("_unit$", cols, value = TRUE)

# tudo que não é taxon, contínuo (min/low/high/max) ou unit = qualitativo
qual_cols <- setdiff(cols, c("taxon", continuous_col, unit_col))

export_qualitative_lookup <- function(traits, qual_cols, path) {
  lookup <- purrr::map_dfr(qual_cols, function(col) {
    x <- as.character(traits[[col]])
    x <- x[!is.na(x) & x != ""]
    if (length(x) == 0) return(NULL)
    tab <- sort(table(x), decreasing = TRUE)
    tibble::tibble(
      variable = col,
      original_value = names(tab),
      n_occurrences = as.integer(tab),
      standardized_value = names(tab)  # coluna pra editar manualmente
    )
  })
  readr::write_csv(lookup, path)
  message(nrow(lookup), " valores únicos exportados em: ", path)
  lookup
}

#evitar que existam diferenças como Branco, branco, branco ,
traits[qual_cols] <- lapply(traits[qual_cols], function(x) trimws(tolower(as.character(x))))

qual_lookup <- export_qualitative_lookup(traits, qual_cols, "3.outputs/qualitative_lookup_mimoseae.csv")

#checar manualmente e gerar um novo arquivo editado pra substituir os nomes 
#quali data to check are inflorescence type, sex type of flower/inflo, habit,
#flower_merosity, filament_color, stamen_count, nectary presence, anther gland presence (i dont now if i will keep stemonozone & stamen tube)
qual_lookup_check <- read.csv("~/Documents/GitHub/bvasconcelos-IC-disparidade-floral/3.outputs/qualitative_lookup_mimoseae_check.csv")

apply_qualitative_lookup <- function(traits, lookup) {
  for (v in unique(lookup$variable)) {
    map <- lookup[lookup$variable == v, ]
    idx <- match(as.character(traits[[v]]), map$original_value)
    traits[[v]][!is.na(idx)] <- map$standardized_value[idx[!is.na(idx)]]
  }
  traits
}

traits <- apply_qualitative_lookup(traits, qual_lookup_check)
## script rodando normalmente até aqui 

### CORRECTING CONTINUOUS DATA ####
# Dados contínuos em traits nem sempre estão corretamente organizados nos seus respectivos 
# min, low, high e max. São as colunas contínuas:
#conferir se todas as colunas continuas de traits estao aqui
continuous_col[!continuous_col %in% colnames(traits)] 

#conferindo os dados 
unique(unlist(unname(traits[continuous_col]))) #tem um cm no meio dos dados, conferir onde que teve esse erro de digitacao 
##uns dados com (0.7) e um cm 
#conferindo onde ta esse cm
traits[apply(traits[continuous_col], 1, function(x) any(x == "cm", na.rm = TRUE)), ]
#198 tetrapleura tetraptera

# função para incluir valores de min e max em low e high, respectivamente
update_trait_values <- function(traits, min_col, low_col, high_col, max_col) {
  
  min_to_low <- which(is.na(traits[[low_col]]) & !is.na(traits[[min_col]])) #se é NA em low e não é em min, min é transferido pra low 
  max_to_high <- which(is.na(traits[[high_col]]) & !is.na(traits[[max_col]])) #se é NA em max e não é em high, high é transferido pra low
  
  if (length(min_to_low)) {
    traits[[low_col]][min_to_low] <- traits[[min_col]][min_to_low] #transfere min para low
    traits[[min_col]][min_to_low] <- NA #exclui em min
  }
  if (length(max_to_high)) {
    traits[[high_col]][max_to_high] <- traits[[max_col]][max_to_high] #transfere max para min
    traits[[max_col]][max_to_high] <- NA #exclui em max
  }
  
  return(traits)
}

cols <- names(traits)
range_traits <- unique(sub("(.+)_min(_|$).*$", "\\1", cols[grepl("_min(_|$)", cols)]))

traits_2 <- traits
#vamos continuar usando o traits_2
for (root in range_traits) {
  min_col  <- grep(paste0("^", root, "_min(_|$)"),  cols, value = TRUE) 
  low_col  <- grep(paste0("^", root, "_low(_|$)"),  cols, value = TRUE)
  high_col <- grep(paste0("^", root, "_high(_|$)"), cols, value = TRUE)
  max_col  <- grep(paste0("^", root, "_max(_|$)"),  cols, value = TRUE)
  
  n <- min(length(min_col), length(low_col), length(high_col), length(max_col))
  
  if (n > 0) {
    for (i in seq_len(n)) { #looping para fazer para todas as ocorrencias
      traits_2 <- update_trait_values(
        traits_2,
        min_col[i],
        low_col[i],
        high_col[i],
        max_col[i] 
      )
    }
  }
}

sum(is.na(traits_2)) #8403
sum(is.na(traits)) #8403

#===========================================#
## Mean values for continuous traits ####
# Verificando se existe algo que tem informação em min mas não tem em low ou se tem em max mas não tem em high
min_cols <- unlist(lapply(range_traits, function(x) grep(paste0("^", x, 
                                                                "_min(_|$)"),  cols, value = TRUE)))

low_cols <- unlist(lapply(range_traits, function(x) grep(paste0("^", x, 
                                                                "_low(_|$)"),  cols, value = TRUE)))

high_cols <- unlist(lapply(range_traits, function(x) grep(paste0("^", x, 
                                                                 "_high(_|$)"), cols, value = TRUE)))

max_cols  <- unlist(lapply(range_traits,  function(x) grep(paste0("^", x, 
                                                                  "_max(_|$)"), cols, value = TRUE)))

any(sapply(seq_along(min_cols), function(i) {
  idx <- !is.na(traits_2[[min_cols[i]]]) & traits_2[[min_cols[i]]] != "" &
    (is.na(traits_2[[low_cols[i]]]) | traits_2[[low_cols[i]]] == "")
  any(idx)
})) #verifica se há algum caso que tem NA ou é vazio na coluna low, mas tem dados na coluna min

any(sapply(seq_along(max_cols), function(i) {
  idx <- !is.na(traits_2[[max_cols[i]]]) & traits_2[[max_cols[i]]] != "" &
    (is.na(traits_2[[high_cols[i]]]) | traits_2[[high_cols[i]]] == "")
  any(idx)
})) #verifica se há algum caso que tem NA ou é vazio na coluna high, mas tem dados na coluna max

# já que não há nada que tenha em min e max que não tenha valores em low e high (ou seja, FALSE foi retornado), 
# vou remover as colunas min e max e tirar a média entre low e high (se só houver apenas um valor, 
# ele será usado)

traits_2 #221 obs, 79 variables
sum(is.na(traits_2)) #8396

for (i in seq_along(low_cols)) {
  
  low  <- as.numeric(traits_2[[low_cols[i]]])
  high <- as.numeric(traits_2[[high_cols[i]]])
  
  mean_col <- sub("_low$", "_mean", low_cols[i])
  
  traits_2[[mean_col]] <- rowMeans(
    cbind(low, high),
    na.rm = TRUE
  )
}

#Verificando
sum(is.na(traits_2[colnames(traits)])) #8396, o mesmo que antes, então não foram gerados NAs ao estimar a média

all.equal(traits_2$height_mean, rowMeans(cbind(as.numeric(traits_2$height_low),
                                               as.numeric(traits_2$height_high)),na.rm = TRUE))

#===================================#
## Keeping only mean columns ####
# removendo colunas com low, min, max, high (manter só mean)
continuous_col <- c(min_cols, low_cols, high_cols, max_cols)
traits_3 <- traits_2

traits_3 <- traits_3[!colnames(traits_3) %in% continuous_col]
#43 variaveis e 221obs
sum(is.na(traits_3)) #2305

### Unit standardization ####
traits_2 <- traits_3
remove(traits_3)
sum(is.na(traits_2)) #2305

## Correcting typos ####
traits_3 <- traits_2

#corrigindo um erro na escrita
colnames(traits_3) <- sub("calyx_lobe_length_unit.", "calyx_lobe_length_unit", colnames(traits_3))
cols <- names(traits_3)
range_traits <- unique(sub("(.+)_mean(_|$).*$", "\\1", cols[grepl("_mean(_|$)", cols)]))

unit_col <- paste(range_traits, "unit", sep = "_") #colunas com unit
all(unit_col %in% colnames(traits_3)) #todas colunas de unit_col está em traits

unique(unname(unlist(lapply(traits_3[unit_col], function (x) unique(x)))))
#ajuste pq tem um M e um ""

traits_3[unit_col] <- lapply(traits_3[unit_col], function(x) {
  x <- tolower(trimws(x))
  x[x == ""] <- NA
  x
})

unique(unname(unlist(lapply(traits_3[unit_col], function (x) unique(x)))))
#parece estar tudo certo "m"  NA "cm" "mm"
traits_2 <- traits_3
remove(traits_3)

# Applying unit standardization ####
traits_3 <- traits_2

mean_cols <- names(traits_3)[grepl("_mean$", names(traits_3))]

for (var in mean_cols) {
  
  unit_var <- sub("_mean$", "_unit", var)
  
  if (!unit_var %in% names(traits_3)) next
  
  traits_3[[unit_var]] <- as.character(traits_3[[unit_var]])
  traits_3[[unit_var]] <- trimws(tolower(traits_3[[unit_var]]))
  
  idx_m  <- traits_3[[unit_var]] %in% "m"
  idx_dm <- traits_3[[unit_var]] %in% "dm"
  idx_mm <- traits_3[[unit_var]] %in% "mm"
  
  traits_3[[var]][idx_m]  <- traits_3[[var]][idx_m] * 100
  traits_3[[var]][idx_dm] <- traits_3[[var]][idx_dm] * 10
  traits_3[[var]][idx_mm] <- traits_3[[var]][idx_mm] / 10
  
  idx_valid <- !is.na(traits_3[[unit_var]]) & traits_3[[unit_var]] != ""
  traits_3[[unit_var]][idx_valid] <- "cm"
}

#tem que ter a mesma soma de NA
sum(is.na(traits_2)) #3238
sum(is.na(traits_3)) #3238


which(traits_2$inflorescence_length_unit == "mm")[3]
traits_2[3,"inflorescence_length_unit"]
traits_2[3,"inflorescence_length_mean"] #16.5
#precisa ser 16.5/10
traits_3[3,"inflorescence_length_mean"] #1.65

which(traits_2$inflorescence_length_unit == "cm")[7]
traits_2[5,"inflorescence_length_unit"]
traits_2[5,"inflorescence_length_mean"] #3.25
#precisa ser 3.25
traits_3[5,"inflorescence_length_mean"]

#checando se algum NA foi introduzido
which(is.na(traits_3[mean_cols]) & !is.na(traits_2[mean_cols]) == T)

#excluir coluna com unidade
traits_3 <- traits_3[, !names(traits_3) %in% unit_col]

#pelos testes, parece estar tudo ok. podemos remover as colunas descrevendo as unidades
traits_2 <- traits_3
remove(traits_3)

cleaned_traits <- traits_2

#tem alguns NaN, vou limpar pra virar NA
# limpar NaN -> NA
cleaned_traits[] <- lapply(cleaned_traits, function(x) {
  if (is.numeric(x)) {
    x[is.nan(x)] <- NA
  }
  x
})
##checando o funcionamento, script funcional até aqui 07/09/2026

#check dataset completeness
#checking traits with less than 15% completeness
traits_percent <- colMeans(!is.na(cleaned_traits)) * 100
names(traits_percent[traits_percent < 15])

#pedicel width mean (vou excluir esse)
traits_percent_original <- colMeans(!is.na(traits)) * 100
names(traits_percent_original[traits_percent_original < 15])
trait_completeness <- colMeans(!is.na(cleaned_traits[,-1])) * 100

cleaned_traits <- cleaned_traits %>%
  select(-pedicel_width_mean)

##to montando o script em ordem entao nao vou salvar o dataset agora
write.csv(cleaned_traits, "3.outputs/morphological_dataset_clean.csv", row.names = F)
# read.csv("3.outputs/morphological_dataset_clean.csv)
traits <- cleaned_traits

#221 obs & 30 variables

## Traits selected for disparity analyses
traits_selected <- traits %>%
  select(
    taxon,
    inflorescence_type,
    flower_merosity,
    stamen_count,
    anther_gland_presence,
    nectary_presence,
    inflorescence_length_mean,
    inflorescence_peduncle_length_mean,
    calyx_length_mean,
    corolla_length_mean,
    corolla_lobe_length_mean,
    pedicel_length_mean,
    filament_length_mean)

dim(traits_selected)
#221 13 
str(traits_selected)
                            
#optei por usar max e min do numero de estames
# Separa os valores usando "-"
stamen_range <- strsplit(traits_selected$stamen_count, "-")

# Extrai o mínimo
traits_selected$stamen_min <- sapply(
  stamen_range,
  function(x) as.numeric(x[1])
)

# Extrai o máximo
traits_selected$stamen_max <- sapply(
  stamen_range,
  function(x) {
    if (length(x) == 2) {
      as.numeric(x[2])
    } else {
      as.numeric(x[1])
    }
  }
)

## calcular média do número de estames
traits_selected$stamen_count <- rowMeans(
  traits_selected[, c("stamen_min", "stamen_max")],
  na.rm = TRUE
)

#excluir maximo e minimo
traits_selected <- traits_selected %>%
  select(-stamen_min, -stamen_max)

#apareceram alguns NaN
traits_selected[] <- lapply(traits_selected, function(x) {
  if (is.numeric(x)) {
    x[is.nan(x)] <- NA
  }
  x
})

#221 obs e 13 variaveis!!
#11/09

##definir os caracteres pra analise
inflo_traits <- c(
  "inflorescence_type",
  "inflorescence_length_mean",
  "inflorescence_peduncle_mean"
)

flower_traits <- c(
  "flower_merosity",
  "stamen_count",
  "anther_gland_presence",
  "nectary_presence",
  "calyx_length_mean",
  "corolla_length_mean",
  "corolla_lobe_length_mean",
  "pedicel_length_mean",
  "filament_length_mean"
)

#turn species names as a rownames
traits_selected <- traits_selected %>%
  rename(species = taxon)

traits_matrix <- traits_selected %>%
  tibble::column_to_rownames("species")

# Transformando células vazias em NA
traits_matrix[] <- lapply(traits_matrix, function(x) {
  if (is.character(x)) {
    x <- trimws(x)
    x[x == ""] <- NA
  }
  x
})

#definindo quais sao as colunas continuas
continuous_cols <- c( grep("_mean$", colnames(traits_matrix), value = TRUE), "stamen_count" ) 
# Definindo quais são as colunas categóricas 
categorical_cols <- setdiff(colnames(traits_matrix), continuous_cols)

#read phylogenetic tree and ecological data
tree <- read.tree("4.trees/mimosoid_calibrated_clean_updated.tre")
## prune phylogeny
tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, traits_selected$species))
            
## set seed for replicability
set.seed(7) 

##Input de dados com Rphylopars e Moda/Media
#tem alguns traços que tao com valor ausente mas nao é NA, ai nao ta imputando 

### A. IMPUTAÇÃO COM MODA / MÉDIA
traits_sub <- traits_matrix %>%
  mutate(across(all_of(categorical_cols), as.factor)) %>%
  mutate(across(all_of(continuous_cols), as.numeric))

get_moda <- function(x) {
  ux <- unique(x[!is.na(x)])
  ux[which.max(tabulate(match(x, ux)))]
}
#criando um dataset pra aplicar moda e media
traits_modemean <- traits_sub

traits_modemean[continuous_cols] <- lapply(traits_modemean[continuous_cols], function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
})

traits_modemean[categorical_cols] <- lapply(traits_modemean[categorical_cols], function(x) {
  x[is.na(x)] <- get_moda(x)
  x
})

stopifnot(sum(sapply(traits_modemean, function(x) sum(is.na(x)))) == 0)


# log nos traços contínuos que fugirem de normalidade (ajuste conforme shapiro_tab)
traits_sub_log <- traits_modemean %>%
  mutate(across(all_of(continuous_cols), log))

#substituindo as colunas
## Continuas: preenchidas com a MEDIA da coluna.
## Categoricas (factor): preenchidas com a MODA da coluna.
## checagem
sum(sapply(traits_matrix, function(x) sum(is.na(x)))) #947 NAs no começo
sum(sapply(traits_modemean,function(x) sum(is.na(x)))) #aqui é 0

# B. IMPUTAÇÃO COM Rphylopars (só traços contínuos)
# ============================================================
# phylopars precisa de data.frame com coluna "species" + as contínuas,
# NA onde faltar dado — não apenas os nomes das colunas
library(Rphylopars)
library(tibble)

(setdiff(traits_selected$species, tree$tip.label))
#Senegalia catechu/Senegalia chundra 
#Senegalia caesia era p ser Senegalia intsia

traits_selected <- traits_selected %>%
  filter(species %in% tree$tip.label)

phylopars_input <- traits_selected %>%
  select(species, all_of(continuous_cols))

phylopars_fit <- phylopars(
  trait_data       = phylopars_input,
  tree             = tree_pruned,
  model            = "BM",
  pheno_error      = TRUE,
  phylo_correlated = TRUE,
  pheno_correlated = TRUE
)

#excluir senegalia catechu e caesiaa?? 
n_tip <- length(tree_pruned$tip.label)
imputed_cont <- phylopars_fit$anc_recon[1:n_tip, continuous_cols, drop = FALSE]

# checagem defensiva de ordem antes de rotular
name_order_ok <- identical(rownames(imputed_cont), tree_pruned$tip.label)
if (!name_order_ok) {
  imputed_cont <- imputed_cont[match(tree_pruned$tip.label, rownames(imputed_cont)), ]
  stopifnot(identical(rownames(imputed_cont), tree_pruned$tip.label))
}

stopifnot(sum(is.na(imputed_cont)) == 0)

# categóricas entram pela moda (Rphylopars não modela discreto/categórico)
traits_phylo <- as.data.frame(imputed_cont) %>%
  rownames_to_column("species") %>%
  left_join(
    traits_modemean %>% rownames_to_column("species") %>% select(species, all_of(categorical_cols)),
    by = "species"
  ) %>%
  column_to_rownames("species")

stopifnot(sum(sapply(traits_phylo, function(x) sum(is.na(x)))) == 0)
#traits_phylo are dataset with phylogenetic input
### ====================== ######

##Gower distance x PCoA =====------ 
library(cluster)
gower_phylo <- daisy(traits_phylo, metric = "gower") ##cluster package
gower_modemean <- daisy(traits_modemean, metric = "gower")

## checagem: quantos pares tem NA na distancia (caso alguma linha nao compartilhe
## nenhuma variavel observada com outra - daria distancia NA)
#como fiz a imputação vai dar 0 
sum(is.na(as.matrix(gower_phylo)))

## PCoA com a matriz de gower =====================================
pcoa_res <- pcoa(gower_phylo)  #ape package

## garantir ordem identica a arvore 
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

##coord fixed to adjust scale 
# objeto dispRity a partir dos eixos da PCoA
library(dispRity)
pcoa_axes <- pcoa_res$vectors

disp_obj <- custom.subsets(pcoa_axes, group = list(all = rownames(pcoa_axes)))

sov <- dispRity(disp_obj, metric = c(sum, variances))
sor <- dispRity(disp_obj, metric = c(sum, ranges))
mpd <- dispRity(disp_obj, metric = c(mean, pairwise.dist))

#PCA Hill-Smith =======
library(ade4)
library(adegraphics)

hs <- dudi.hillsmith(traits_phylo,
               scannf = TRUE, nf = 2)
#select 2 

hs$eig
axes_contribution <- 100*hs$eig/sum(hs$eig)

#plot flowers in morpho space
plot(hs$li[,1],
     hs$li[,2],
     xlab = "pc1",
     ylab = "pc2",
     pch = 19)

text(hs_phylo$li[,1],hs_phylo$li[,2],labels = traits_phylo$species, pos = 1)



##calcular metricas de disparidade =======
## SV, SR, MPD
## --- SV (sum of variances) e SR (sum of ranges) via dispRity,
##     calculadas sobre TODOS os eixos retidos do PCoA (nao so os 2 primeiros,
##     para nao perder disparidade que esta em eixos de ordem maior) ---

gower_dist <- cluster::daisy(trait_df_clean, metric = "gower")

# PCoA para obter um espaço multivariado contínuo (necessário p/ SOV e SOR)
pcoa_res <- ape::pcoa(gower_dist)
pcoa_axes <- pcoa_res$vectors

# objeto dispRity a partir dos eixos da PCoA
disp_obj <- custom.subsets(pcoa_axes, group = list(all = rownames(pcoa_axes)))

sov <- dispRity(disp_obj, metric = c(sum, variances))
sor <- dispRity(disp_obj, metric = c(sum, ranges))

gower_mat <- as.matrix(gower_dist)
disp_obj_dist <- custom.subsets(gower_mat, group = list(all = rownames(gower_mat)))
mpd <- dispRity(disp_obj_dist, metric = c(mean, pairwise.dist))


##plot clades ======


#plot ecomorphospace ======

### sinal filogenetico?
# δ (delta) statistic Borges
                     
##phylomorphospace ======
# phylomorphospace exige matrix numerica pura, na mesma ordem da arvore
library(phytools)
scores_mat <- as.matrix(scores_pcoa)

phylomorphospace(
  tree_pruned,
  scores_mat,
  xlab = paste0("PC1 (", var_pc1, "%)"),
  ylab = paste0("PC2 (", var_pc2, "%)"),
  label = "off"
)

