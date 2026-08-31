if (!require(librarian)) install.packages("librarian")
librarian::shelf(dplyr, purrr, readr, stringr, tidyr, tibble,
                 cluster, ape, vegan, ggplot2, readr)

morpho_data<- read_csv("~/Documents/GitHub/bvasconcelos-IC-disparidade-floral/1.datasets/mimoseae_subset_clean.csv")
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
qual_lookup_check <- read.csv("3.outputs/qualitative_lookup_mimoseaee_check.csv")

apply_qualitative_lookup <- function(traits, lookup) {
  for (v in unique(lookup$variable)) {
    map <- lookup[lookup$variable == v, ]
    idx <- match(as.character(traits[[v]]), map$original_value)
    traits[[v]][!is.na(idx)] <- map$standardized_value[idx[!is.na(idx)]]
  }
  traits
}

traits <- apply_qualitative_lookup(traits, qual_lookup_check)

sapply(qual_cols, function(v) length(unique(traits[[v]])))  # comparar antes/depois

### CORRECTING CONTINUOUS DATA ####
# Dados contínuos em traits nem sempre estão corretamente organizados nos seus respectivos 
# min, low, high e max. São as colunas contínuas:
#conferir se todas as colunas continuas de traits estao aqui
continuous_col[!continuous_col %in% colnames(traits)] 

#conferindo os dados 
unique(unlist(unname(traits[continuous_col]))) #tem um cm no meio dos dados, conferir onde que teve esse erro de digitacao 
## +- em um dado, uns dados com (0.7) e acho que só?

#conferindo onde ta esse cm
traits[apply(traits[continuous_col], 1, function(x) any(x == "cm", na.rm = TRUE)), ]

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

sum(is.na(traits_2)) #10188
sum(is.na(traits)) #10188

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

traits_2 #33 obs, 79 variables
sum(is.na(traits_2)) #1428

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
sum(is.na(traits_2[colnames(traits)])) #1428, o mesmo que antes, então não foram gerados NAs ao estimar a média

all.equal(traits_2$height_mean, rowMeans(cbind(as.numeric(traits_2$height_low),
                                               as.numeric(traits_2$height_high)),na.rm = TRUE))

#===================================#
## Keeping only mean columns ####
# removendo colunas com low, min, max, high (manter só mean)
continuous_col <- c(min_cols, low_cols, high_cols, max_cols)
traits_3 <- traits_2

traits_3 <- traits_3[!colnames(traits_3) %in% continuous_col]
#43 variaveis e 33obs
sum(is.na(traits_3)) #584

### Unit standardization ####
traits_2 <- traits_3
remove(traits_3)
sum(is.na(traits_2))

## Correcting typos ####
traits_3 <- traits_2

#corrigindo um erro na escrita
colnames(traits_3) <- sub("calyx_lobe_length_unit.", "calyx_lobe_length_unit", colnames(traits_3))

cols <- names(traits_3)
range_traits <- unique(sub("(.+)_mean(_|$).*$", "\\1", cols[grepl("_mean(_|$)", cols)]))

unit_col <- paste(range_traits, "unit", sep = "_") #colunas com unit
all(unit_col %in% colnames(traits_3)) #todas colunas de unit_col está em traits

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
sum(is.na(traits_2)) #584
sum(is.na(traits_3)) #584


which(traits_2$inflorescence_length_unit == "mm")[3]
traits_2[10,"inflorescence_length_unit"]
traits_2[10,"inflorescence_length_mean"] #2.5
#precisa ser 2.5/10
traits_3[10,"inflorescence_length_mean"]

which(traits_2$inflorescence_length_unit == "cm")[7]
traits_2[20,"inflorescence_length_unit"]
traits_2[20,"inflorescence_length_mean"] #1.7
#precisa ser 1.7
traits_3[20,"inflorescence_length_mean"]

#checando se algum NA foi introduzido
which(is.na(traits_3[mean_cols]) & !is.na(traits_2[mean_cols]) == T)

#excluir coluna com unidade
traits_3 <- traits_3[, !names(traits_3) %in% unit_col]

#pelos testes, parece estar tudo ok. podemos remover as colunas descrevendo as unidades
traits_2 <- traits_3
remove(traits_3)

cleaned_traits <- traits_2
#check dataset completeness
#checking traits with less than 15% completeness
traits_percent <- colMeans(!is.na(cleaned_traits)) * 100
names(traits_percent[traits_percent < 15])
traits_percent_original <- colMeans(!is.na(traits)) * 100
names(traits_percent_original[traits_percent_original < 15])
trait_completeness <- colMeans(!is.na(cleaned_traits[,-1])) * 100

head(cleaned_traits)

write.csv(cleaned_traits, "3.outputs/morphological_dataset_clean.csv", row.names = F)

