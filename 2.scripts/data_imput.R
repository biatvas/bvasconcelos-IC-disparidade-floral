###IMPUTAÇÃO DOS DADOS
# UM NOVO SCRIPT PARA SELEÇÃO DAS VARIAVEIS E PODEMOS INCLUIR ESTATISTICAS
# GERAIS DOS TRAÇOS
library(phytools)
library(ape)
#aqui quero ter uma tabela com dados imputados com os nomes das espécies ordenados de acordo com a filogenia 
#read morpho data
setwd("Documents/GitHub/bvasconcelos-IC-disparidade-floral/")
traits <- read.csv("3.outputs/morphological_dataset_treatment.csv")
#read phylogenetic tree and ecological data
tree <- read.tree("4.trees/mimosoid_calibrated_clean_updated.tre")
## prune phylogeny
tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, traits$taxon))

traits_ordered <- traits[match(tree_pruned$tip.label, traits$taxon),]

traits <- traits_ordered %>%
  rename(species = taxon)

clade_info <- traits %>%
  select(species, clade)

traits <- traits %>%
  select(
    species,
    inflorescence_type,
    flower_merosity,
    stamen_count,
    filament_color,
    anther_gland_presence,
    nectary,
    sex_type,
    dimorphic_flower,
    inflorescence_length_mean,
    inflorescence_peduncle_length_mean,
    calyx_length_mean,
    corolla_length_mean,
    corolla_lobe_length_mean,
    filament_length_mean)

dim(traits)
#219 15

#turn species names as a rownames
traits_matrix <- traits %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames("species")

#definindo quais sao as colunas continuas
continuous_cols <- c( grep("_mean$", colnames(traits_matrix), value = TRUE), "stamen_count" ) 
# Definindo quais são as colunas categóricas 
categorical_cols <- setdiff(colnames(traits_matrix), continuous_cols)

## set seed for replicability
set.seed(7) 

##Input de dados com Rphylopars e Moda/Media
#tem alguns traços que tao com valor ausente mas nao é NA, ai nao ta imputando 
### A. IMPUTAÇÃO COM MODA / MÉDIA
traits <- traits_matrix %>%
  mutate(across(all_of(categorical_cols), as.factor)) %>%
  mutate(across(all_of(continuous_cols), as.numeric))

traits_modemean <- traits

get_moda <- function(x) {
  ux <- unique(x[!is.na(x)])
  ux[which.max(tabulate(match(x, ux)))]
}

traits_modemean[continuous_cols] <- lapply(traits_modemean[continuous_cols], function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
})

traits_modemean[categorical_cols] <- lapply(traits_modemean[categorical_cols], function(x) {
  x[is.na(x)] <- get_moda(x)
  x
})

stopifnot(sum(sapply(traits_modemean, function(x) sum(is.na(x)))) == 0)

#traits modemean here are with data imputed
# podemos considerar a moda de linhagens filogeneticamente proximas
# como genero.
# vamos testar fazer moda e media pra grupos dos clados e/ou por generos
genus <- data.frame("genus" = sub("_.*", "", row.names(traits)))

traits_modemean_2 <- traits %>%
 tibble::rownames_to_column("species") %>% #cria a col species
 dplyr::mutate(genus = sub("_.*", "", species)) %>% #cria a col genus, selecionando apenas o primeiro nome antes de _ de species
 dplyr::group_by(genus) %>% #agrupa por genero
 dplyr::mutate(
   dplyr::across(
     dplyr::all_of(categorical_cols), #considera apenas as variaveis em categorical_cols
     ~ {
       x <- .x
       x[is.na(x)] <- get_moda(x) #usando a funcao criada acima
       x
     }
   )
 ) %>%
 dplyr::ungroup() %>%
 tibble::column_to_rownames("species")

#all(genus$genus %in% sub("_.*","", tree_pruned$tip.label))

# funciona, mas retorna NA pros generos com apenas uma especie no
# dataset e que eh NA pra variavel. Ou seja, nesses casos teriamos que fazer o
# mesmo procedimento mas considerando generos proximos

#talvez fazer algo como: 
# se a ocorrencia de um nome em genus é 1, entao selecionar o 
# genero que ocorre logo antes do nome em sub("_.*","", 
# tree_pruned$tip.label)

traits_modemean_2 <- traits_modemean_2 %>%
  tibble::rownames_to_column("species") %>%
  left_join(clade_info, by = "species")

traits_modemean_3 <- traits_modemean_2 %>%
  dplyr::group_by(clade) %>% #agrupa por clado
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(categorical_cols),
      ~ {
        x <- .x
        x[is.na(x)] <- get_moda(x)
        x
      }
    )
  ) %>%
  dplyr::ungroup() %>%
  tibble::column_to_rownames("species")

#imputar media dos continuos 
traits_modemean_3[continuous_cols] <- lapply(traits_modemean_3[continuous_cols], function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
})

#checar se imputou em tudo 
sum(sapply(traits_modemean_3, function(x) sum(is.na(x))))

# log nos traços contínuos que fugirem de normalidade 
traits_log <- traits_modemean %>%
  mutate(across(all_of(continuous_cols), log))

#substituindo as colunas
## Continuas: preenchidas com a MEDIA da coluna.
## Categoricas (factor): preenchidas com a MODA da coluna.
## checagem
sum(sapply(traits_matrix, function(x) sum(is.na(x)))) #947 NAs no começo
sum(sapply(traits_modemean,function(x) sum(is.na(x)))) #aqui é 0
sum(sapply(traits_log,function(x) sum(is.na(x)))) #aqui é 0 tbm

#isso ja foi feito antes 
# #ordenando traits_modemean pra filogenia
# traits_modemean <- traits_modemean[match(tree_pruned$tip.label, 
#                                          row.names(traits_modemean)),]


#identical(row.names(traits_modemean), tree_pruned$tip.label)
#retorna TRUE 

# B. IMPUTAÇÃO COM Rphylopars (só traços contínuos)
# ============================================================
# phylopars precisa de data.frame com coluna "species" + as contínuas,
# NA onde faltar dado — não apenas os nomes das colunas
library(Rphylopars)
library(tibble)

#quando ordenou baseado na filogenia, isso ja foi resolvido 
# (setdiff(traits_modemean_3$species, tree$tip.label))
#Senegalia catechu/Senegalia chundra 
#Senegalia caesia era p ser Senegalia intsia

phylopars_input <- traits %>%
  tibble::rownames_to_column("species") %>%
  select(species, all_of(continuous_cols)) 

#ordenando as especies para ter a mesma ordem da filogenia
#all(phylopars_input$species %in% tree_pruned$tip.label)
#retorna T
# phylopars_fit <- phylopars(
#   trait_data       = phylopars_input,
#   tree             = tree_pruned,
#   model            = "BM",
#   pheno_error      = TRUE,
#   phylo_correlated = TRUE,
#   pheno_correlated = TRUE
# )

#sem assumir correlacao entre tracos e variacao intraespecifica
phylopars_fit_no_cor <- phylopars(
  trait_data       = phylopars_input,
  tree             = tree_pruned,
  model            = "BM",
  pheno_error      = F,
  phylo_correlated = F,
  pheno_correlated = F
)

#checando
# View(data.frame(phylopars_fit$anc_recon[1:219,1],
#            phylopars_fit_no_cor$anc_recon[1:219,1],
#      phylopars_input_ordered$inflorescence_length_mean[1:219]))


n_tip <- length(tree_pruned$tip.label)
imputed_cont <- phylopars_fit_no_cor$anc_recon[1:n_tip, continuous_cols, drop = FALSE]

# checagem defensiva de ordem antes de rotular
# name_order_ok <- identical(rownames(imputed_cont), tree_pruned$tip.label)
# if (!name_order_ok) {
#   imputed_cont <- imputed_cont[match(tree_pruned$tip.label, rownames(imputed_cont)), ]
#   stopifnot(identical(rownames(imputed_cont), tree_pruned$tip.label))
# }
# 
# stopifnot(sum(is.na(imputed_cont)) == 0)

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

# log nos traços contínuos que fugirem de normalidade 
log_phylo <- traits_phylo %>%
  mutate(across(all_of(continuous_cols), log))
### ====================== ######
traits_final <- log_phylo %>%
  tibble::rownames_to_column("species") %>%
  left_join(clade_info, by = "species")

#por fim, salvar o dataset final gerado
write.csv(traits_final, "traits_treatmentimput30092026.csv", row.names = F)
