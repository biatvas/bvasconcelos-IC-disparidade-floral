#============#
# 1. HEAD ####
#============#

if (!require(librarian)) install.packages("librarian"); library("librarian")
librarian::shelf(phytools, dplyr, tidyr, purrr, vegan, tidyverse,ape, stringr, readr)

traits <- read.csv("3.outputs/corrected_continuous_data-20260408.csv")

#========================#
# 2. Merging  columns ####
#========================#

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

#gerando noto dataset para trabalharmos
traits_2 <- traits

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

sum(is.na(traits_2)) #167229
sum(is.na(traits)) #167229

#===========================================#
## 2.1 Mean values for continuous traits ####
#===========================================#

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

traits_2 #3172 obs, 79 variables
sum(is.na(traits_2)) #167229

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
sum(is.na(traits_2[colnames(traits)])) #167229, o mesmo que antes, então não foram gerados NAs ao estimar a média

all.equal(traits_2$height_mean, rowMeans(cbind(as.numeric(traits_2$height_low),
               as.numeric(traits_2$height_high)),na.rm = TRUE))

#===================================#
## 2.2 Keeping only mean columns ####
#===================================#

# removendo colunas com low, min, max, high (manter só mean)
continuous_col <- c(min_cols, low_cols, high_cols, max_cols)
traits_3 <- traits_2

traits_3 <- traits_3[!colnames(traits_3) %in% continuous_col]
#43 variaveis e 3172 obs
sum(is.na(traits_3)) #69858

#============================#
# 3. Unit standardization ####
#============================#
# traits_2_backup <- traits_2

traits_2 <- traits_3
remove(traits_3)
sum(is.na(traits_2)) #69858

#==========================#
## 3.1 Correcting typos ####
#==========================#
traits_3 <- traits_2

#corrigindo um erro na escrita
colnames(traits_3) <- sub("calyx_lobe_length_unit.", "calyx_lobe_length_unit", colnames(traits_3))

cols <- names(traits_3)
range_traits <- unique(sub("(.+)_mean(_|$).*$", "\\1", cols[grepl("_mean(_|$)", cols)]))

unit_col <- paste(range_traits, "unit", sep = "_") #colunas com unit
all(unit_col %in% colnames(traits_3)) #todas colunas de unit_col está em traits

unique(unname(unlist(lapply(traits_3[unit_col], function (x) unique(x))))) #tom? 

which(apply(traits_3[unit_col], 1, function(row) {
  any(row == "tom", na.rm = TRUE)
}) == T)

traits_3[165,]$taxon

#corrigindo
traits_3[165,"pedicel_width_unit"] <- "mm"
#sum(is.na(traits_3)) #69858

unique(unname(unlist(lapply(traits_3[unit_col], function (x) unique(x)))))

traits_2 <- traits_3
remove(traits_3)

#=======================================#
## 3.2 Applying unit standardization ####
#=======================================#
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
sum(is.na(traits_2)) #69858
sum(is.na(traits_3)) #69858

which(traits_2$inflorescence_length_unit == "mm")[3]
traits_2[10,"inflorescence_length_unit"]
traits_2[10,"inflorescence_length_mean"] #1.75
#precisa ser 1.75/10
traits_3[10,"inflorescence_length_mean"]

which(traits_2$inflorescence_length_unit == "cm")[7]
traits_2[95,"inflorescence_length_unit"]
traits_2[95,"inflorescence_length_mean"]
#precisa ser 3
traits_3[95,"inflorescence_length_mean"]

#checando se algum NA foi introduzido
which(is.na(traits_3[mean_cols]) & !is.na(traits_2[mean_cols]) == T)

#pelos testes, parece estar tudo ok. podemos remover as colunas descrevendo as unidades
traits_2 <- traits_3
remove(traits_3)

#traits contain 31 variables
traits_3 <- traits_3[, !names(traits_3) %in% unit_col]
sum(is.na(traits_3)) #54868

#=========================================#
# 4. Saving continuous cleaned dataset ####
#=========================================#

write.csv(traits_3, "3.outputs/continuous_data_cleaned-20260410.csv", row.names = F,
          fileEncoding = "UTF-8")

#============================#
# 5. Dataset completeness ####
#============================#

cleaned_traits <- read.csv("3.outputs/continuous_data_cleaned-20260410.csv")

#checking traits with less than 15% completeness
traits_percent <- colMeans(!is.na(cleaned_traits)) * 100
names(traits_percent[traits_percent < 15])
traits_percent_original <- colMeans(!is.na(traits)) * 100
names(traits_percent_original[traits_percent_original < 15])
trait_completeness <- colMeans(!is.na(cleaned_traits[,-1])) * 100

#checking less than 30% 
names(traits_percent[traits_percent < 30])

# acho que podemos revisar os dados não contínuos tbm e depois a gnt faz a filtragem

# traits_filtered <- cleaned_traits[, c(TRUE, trait_completeness >= 30)]

#now i will keep 30% filtered columns and remove taxa with less than 60% completeness
# species_completeness <- rowMeans(!is.na(traits_filtered[,-1])) * 100
# traits_filtered$species_completeness <- species_completeness
#write.csv(traits_filtered, "4.outputs/traits_filtered.csv", row.names = F)

#remove species less than 60% completeness
# traits_final <- traits_filtered %>%
#   dplyr::filter(species_completeness >= 60) %>%
#   dplyr::select(-species_completeness)

#now we have 1292 species

#write.csv(traits_final,"final_species_filtered_data.csv", row.names = F)

# traits_final <- traits_final %>% 
#   separate(taxon, into = c("genus", "epithet"),
#            sep = "_")
# unique(traits_final$genus)
# traits_final %>% distinct() %>% count(genus, name = "traits_final")
