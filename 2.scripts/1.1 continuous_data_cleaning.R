#============#
# 1. HEAD ####
#============#

if (!require(librarian)) install.packages("librarian"); library("librarian")
librarian::shelf(phytools, dplyr, tidyr, purrr, vegan, tidyverse,ape, stringr, readr)

morphodata <- read.csv("1.datasets/mimoseae_species_clean.csv")

#======================#
# 2. SELECTING COLS ####
#======================#

traits <- cbind("taxon" = morphodata$taxon, morphodata[,6:83])
#3172 obs and 79 variables

#==================================#
# 3. CORRECTING CONTINUOUS DATA ####
#==================================#

# Dados contínuos em traits nem sempre estão corretamente organizados nos seus respectivos 
# min, low, high e max. São as colunas contínuas:

cols <- names(traits)
range_traits <- unique(sub("(.+)_min(_|$).*$", "\\1", cols[grepl("_min(_|$)", cols)]))

min_col <- unlist(lapply(range_traits, function(x) grep(paste0("^", x, 
                                                               "_min(_|$)"),  cols, value = TRUE)))
low_col <- unlist(lapply(range_traits, function(x) grep(paste0("^", x, 
                                                               "_low(_|$)"),  cols, value = TRUE)))
high_col <- unlist(lapply(range_traits, function(x) grep(paste0("^", x, 
                                                                "_high(_|$)"), cols, value = TRUE)))
max_col  <- unlist(lapply(range_traits,  function(x) grep(paste0("^", x, 
                                                                 "_max(_|$)"), cols, value = TRUE)))

all.equal(length(min_col), length(low_col), length(high_col), length(max_col))

continuous_col <- c(min_col, low_col, high_col, max_col)
continuous_col[!continuous_col %in% colnames(traits)] #checando

# Podem estar ter valores imputados como:
unique(unlist(unname(traits[continuous_col]))) #logo, antes de transformações precisam ser corrigidos

# Gerando um outro banco de dados para trabalharmos sem alterar o original
traits_2 <- traits

#=========================#
## 3.1 Removing spaces ####
#=========================#

# Primeiro, vamos remover espaços das colunas contínuas. 

no_space <- as.data.frame(sapply(traits_2[continuous_col], function(x) gsub("\\s+", "", x)))

#verificando se nada foi alterado no nome das colunas ou se NA foram gerados
all(names(no_space) == names(traits_2[continuous_col])) 
sum(is.na(traits_2[continuous_col])) #121684
sum(is.na(no_space)) #121684

#já que tudo parece estar certo, trocando as colunas contínuas para as sem espaços em traits_2
traits_2[continuous_col] <- no_space
sum(is.na(traits)) #166507
sum(is.na(traits_2)) #166507

#============================#
## 3.2 Removing ≤ symbols ####
#============================#

# Abaixo, faremos com que todas as ocorrencias de números iniciadas com ≤ , por exemplo "≤8", 
# sejam copiadas nas colunas de nome high e o simbolo sejá removido

#Novamente, criando um traits_3 para não modificar o que já estava certo antes
traits_3 <- traits_2

for (col in continuous_col) {

  root <- sub("_(min|low|high|max)$", "", col)
  
  col_high <- paste0(root, "_high") #isso aqui é porque os valores serão adicionados em high
  
  if (!(col_high %in% names(traits_3))) next
  
  vals <- traits_3[[col]]
  high_vals <- traits_3[[col_high]]
  
  idx <- grepl("^≤", vals) & (is.na(high_vals) | high_vals == "") #entao se tem esse simbolo em alguma coluna e os valores high são NA, então transfere pra high
  
  traits_3[[col_high]][idx] <- sub("^≤", "", vals[idx])
  traits_3[[col]][idx] <- NA_character_ #a coluna que tinha o símbolo e não era high vira NA
}

#verificando se deu certo:
##antes era:
traits_2$inflorescence_length_min[grepl("^≤", traits_2$inflorescence_length_min)]
which(grepl("^≤", traits_2$inflorescence_length_min) == T)
## agora:
traits_3$inflorescence_length_min[which(grepl("^≤", traits_2$inflorescence_length_min) == T)] #o que pedimos para virar NA
traits_3$inflorescence_length_high[which(grepl("^≤", traits_2$inflorescence_length_min) == T)] #os valores foram adionados em high

sum(is.na(traits_2)) #166507
sum(is.na(traits_3)) #166507

# Já que está tudo certo com traits_3, podemos substituir o dataset
traits_2 <- traits_3
remove(traits_3)

#============================#
## 3.3 Removing ± symbols ####
#============================#

#Se tem o símbolo ±, então o valor é adicionado em low e em high. Onde ocorria antes é removido. 

#criando um dataset para trabalharmos
traits_3 <- traits_2

for (col in continuous_col) {

  root <- sub("_(min|low|high|max)$", "", col)
  
  col_low  <- paste0(root, "_low") #porque incluiremos em low
  col_high <- paste0(root, "_high") #porque incluiremos em high
  
  if (!(col_low %in% names(traits_3)) || !(col_high %in% names(traits_3))) next #para colunas que não são low e high, aplica o definido abaixo
  
  vals <- traits_3[[col]]
  low_vals  <- traits_3[[col_low]]
  high_vals <- traits_3[[col_high]]
  
  #seleciona quem tem ± no início
  idx <- grepl("^±\\s*", vals) &
    (is.na(low_vals) | low_vals == "") &
    (is.na(high_vals) | high_vals == "")
  
  #tira o simbolo
  clean_vals <- sub("^±\\s*", "", vals[idx])
  
  #copia em low e high
  traits_3[[col_low]][idx]  <- clean_vals
  traits_3[[col_high]][idx] <- clean_vals
  
  #apaga onde estavam as ocorrencias de +/-
  traits_3[[col]][idx] <- NA_character_
}

#verificando se deu tudo certo usando inflorescence_length_min como exemplo 
traits_2$inflorescence_length_min[grepl("^±", traits_2$inflorescence_length_min)]
which(grepl("^±", traits_2$inflorescence_length_min) == T)

traits_3$inflorescence_length_min[which(grepl("^±", traits_2$inflorescence_length_min) == T)]
traits_3$inflorescence_length_high[which(grepl("^±", traits_2$inflorescence_length_min) == T)]

sum(is.na(traits_3)) #166506
sum(is.na(traits_2)) #166506

#td certo, então subtituindo os dataset
traits_2 <- traits_3
remove(traits_3)

#===============================#
## 3.4 Removing "()" symbols ####
#===============================#

# agora resolvendo problemas dos que estão entre parenteses, atribuindo esses valores para
# low, min, max ou high, de acordo com a posição dos valores

#dataset para trabalharmos
traits_3 <- traits_2

for (col in continuous_col) {
  
  root <- sub("_(min|low|high|max)$", "", col)
  
  col_min  <- paste0(root, "_min")
  col_low  <- paste0(root, "_low")
  col_high <- paste0(root, "_high")
  col_max  <- paste0(root, "_max")
  col_check <- paste0(root, "_check") #pra checar se deu certo
  
  vals <- traits_3[[col]]
  
  # selecionando quais colunas tem casos de parenteses pra deixar elas na coluna check, caso seja preciso checar
  idx <- grepl("\\(", vals)
  traits_3[[col_check]] <- NA_character_
  traits_3[[col_check]][idx] <- vals[idx]
  
  x <- vals[idx]
  
  #extraindo os numeros dos parenteses no inicio do numero medio, ou seja, valroes que serão incluidos no min
  
  has_min <- grepl("^\\(", x) 
  min_vals <- rep(NA_character_, length(x))
  
  min_vals[has_min] <- sub("^\\(([^\\)]+)\\).*", "\\1", x[has_min])
  min_vals <- sub("-$", "", min_vals)
  
  traits_3[[col_min]][idx] <- min_vals
  
  #extraindo os numeros dos parenteses no final do numero medio, ou seja, valroes que serão incluidos no max
  
  has_max <- grepl("\\)$", x)
  max_vals <- rep(NA_character_, length(x))
  
  max_vals[has_max] <- sub(".*\\(([^\\)]+)\\)$", "\\1", x[has_max])
  max_vals <- sub("^-", "", max_vals)
  
  traits_3[[col_max]][idx] <- max_vals
  
  #removendo parenteses dos que trocamos  
  core <- gsub("\\([^\\)]+\\)", "", x)
  core <- gsub("--+", "-", core)
  core <- gsub("^-|-$", "", core)
  
  #definindo quem é low e high
  low_vals  <- rep(NA_character_, length(core))
  high_vals <- rep(NA_character_, length(core))
  
  has_range <- grepl("-", core) #agora dos dos que tem parenteses, que tbm estao entre traços, ou seja, low e high 
  
  low_vals[has_range]  <- sub("^(\\d+\\.?\\d*).*", "\\1", core[has_range])
  high_vals[has_range] <- sub(".*-(\\d+\\.?\\d*)$", "\\1", core[has_range])
  #outras inconsistencias com hifen
  no_range <- !has_range & core != ""
  low_vals[no_range] <- core[no_range]
  
  traits_3[[col_low]][idx]  <- low_vals
  traits_3[[col_high]][idx] <- high_vals
  
}

#verificando se deu certo:
sum(grepl("\\(", traits_2$inflorescence_length_min)) #tinha 7 ocorrencias com ()
sapply(traits_2[continuous_col], function(x) sum(grepl("\\(", x))) #todas na coluna
#inflorescence_length_min

traits_2$inflorescence_length_min[grepl("\\)", traits_2$inflorescence_length_min)]
which(grepl("\\(", traits_2$inflorescence_length_min) == T)

#valores que antes estavam entre parenteses e antes do número médio, agora estão em min
traits_3$inflorescence_length_min[which(grepl("\\(", traits_2$inflorescence_length_min) == T)]
#valores que antes estavam entre parenteses e depois do número médio, agora estão em max
traits_3$inflorescence_length_max[which(grepl("\\(", traits_2$inflorescence_length_min) == T)]

# valores que não estavam entre parenteses e que se referiam ao low, agora estçai em low
traits_3$inflorescence_length_low[which(grepl("\\(", traits_2$inflorescence_length_min) == T)]

# valores que não estavam entre parenteses e que se referiam ao high, agora estão em high
traits_3$inflorescence_length_high[which(grepl("\\(", traits_2$inflorescence_length_min) == T)]

#de traits_3 selecionando apenas as colunas que tbm estão em traits_2, ou seja, removendo os check
traits_3 <- select(traits_3, -c((names(traits_3[!colnames(traits_3) %in% colnames(traits_2)]))))

#verificando se tá td certo
all(names(traits_3) %in% names(traits_2))
#é pra ter diminuido o numero de NA, já que dividimos uma coluna em outras
sum(is.na(traits_3)) #166489
sum(is.na(traits_2)) #166506

#tudo parece estar certo, substituindo os dataset
traits_2 <- traits_3
remove(traits_3)

#verificando o que falta arrumar
unique(unlist(unname(traits_2[continuous_col])))

#==============================#
## 3.5 Removing "-" symbols ####
#==============================#

#Tratando os numeros entre "-", atribuindo os valores para low ou para high

#dataset para trabalharmos
traits_3 <- traits_2

for (col in continuous_col) {
  
  root <- sub("_(min|low|high|max)$", "", col)
  
  col_low  <- paste0(root, "_low")
  col_high <- paste0(root, "_high")
  
  vals <- traits_3[[col]]
  
  idx <- grepl("^\\d+\\.?\\d*-\\d+\\.?\\d*$", vals)
  
  low_vals  <- sub("^(\\d+\\.?\\d*)-.*", "\\1", vals[idx])
  high_vals <- sub(".*-(\\d+\\.?\\d*)$", "\\1", vals[idx])
  
  traits_3[[col_low]][idx]  <- low_vals
  traits_3[[col_high]][idx] <- high_vals
  
  traits_3[[col]][idx] <- NA_character_
  
}

#verificando
traits_2$inflorescence_length_min[grepl("-", traits_2$inflorescence_length_min)]
which(grepl("-", traits_2$inflorescence_length_min) == T)

traits_3$inflorescence_length_low[which(grepl("-", traits_2$inflorescence_length_min) == T)]
traits_3$inflorescence_length_high[which(grepl("-", traits_2$inflorescence_length_min) == T)]

#novamente, é pra diminuir o número de NA, porque estão dividindo o que tinha em apenas uma coluna em duas
sum(is.na(traits_3)) #166440
sum(is.na(traits_2)) #166489

#porque aparentemente está tudo certo, substituindo os datasets
traits_2 <- traits_3
remove(traits_3)
all(names(traits) %in% names(traits_2))

#=====================#
# 4.SAVIND DATASET ####
#=====================#

write.csv(traits_2, "3.outputs/corrected_continuous_data-20260408.csv", 
          fileEncoding = "UTF-8", row.names = F)
