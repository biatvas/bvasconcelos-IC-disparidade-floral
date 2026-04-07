#======================#
#----- setup ---------- 
#2026.04.07
setwd("~/Documents/Labis/Dados") #defining work directory
if (!require(librarian)) install.packages("librarian"); library("librarian")
librarian::shelf(phytools, dplyr, tidyr, purrr, vegan, tidyverse,ape, stringr, readr) #installing and/or loading packages
#===============================================================#
# 1. read complete dataset
#===============================================================#
morphodata <- read.csv("1.datasets/mimoseae_species_cleaned.csv")
#3174 obs and 84 variables

#i will make a report with species just from phylogeny 
morphodata %>% 
  separate(taxon, into = c("genus", "epithet"),
           sep = "_") %>% 
  distinct() %>% 
  count(genus, name = "morpho_taxon")


colnames(morphodata)
#i will exclude the columns that are not traits 
#continue avaliating the traits completeness
traits <- cbind("taxon" = morphodata$taxon, morphodata[,6:83])

## check data
str(traits)

### some numeric trait columns seems to present non-numeric data...

### transform non-numeric data into numeric for now, to check missing data
cols <- c("inflorescence_length_min", 
          "inflorescence_length_low",
          "inflorescence_width_low",
          "inflorescence_peduncle_length_min",
          "inflorescence_peduncle_length_low",
          "pedicel_length_min",
          "pedicel_length_low",
          "calyx_length_min",
          "calyx_length_low",
          "calyx_lobe_length_high",
          "corolla_length_low",
          "corolla_lobe_length_low",
          "corolla_lobe_width_min",
          "corolla_lobe_width_low",
          "filament_length_min",
          "filament_length_low")
traits[cols] <- lapply(traits[cols], \(v) ifelse(is.na(v), NA, 1))


# Data treatment ==============================================================

## resolve wrong placement of cont trait values -------------------------------

## data extraction had problems dealing with min-low-high-max
update_trait_values <- function(traits, min_col, low_col, high_col, max_col) {
  ### find indices for min→low and max→high separately
  min_to_low <- which(is.na(traits[[low_col]]) & !is.na(traits[[min_col]]))
  max_to_high <- which(is.na(traits[[high_col]]) & !is.na(traits[[max_col]]))
  
  ### only update where needed, then set min/max to NA
  if (length(min_to_low)) {
    traits[[low_col]][min_to_low] <- traits[[min_col]][min_to_low]
    traits[[min_col]][min_to_low] <- NA
  }
  if (length(max_to_high)) {
    traits[[high_col]][max_to_high] <- traits[[max_col]][max_to_high]
    traits[[max_col]][max_to_high] <- NA
  }
  
  return(traits)
}

### find all columns ending with _min_, _low_, _high_, or _max_ (followed by anything)
cols <- names(traits)
range_traits <- unique(sub("(.+)_min(_|$).*$", "\\1", cols[grepl("_min(_|$)", cols)]))

for (root in range_traits) {
  min_col  <- grep(paste0("^", root, "_min(_|$)"),  cols, value = TRUE)
  low_col  <- grep(paste0("^", root, "_low(_|$)"),  cols, value = TRUE)
  high_col <- grep(paste0("^", root, "_high(_|$)"), cols, value = TRUE)
  max_col  <- grep(paste0("^", root, "_max(_|$)"),  cols, value = TRUE)
  
  ### only proceed if all columns exist for this root
  if (length(min_col) && length(low_col) && length(high_col) && length(max_col)) {
    traits <- update_trait_values(traits, min_col[1], low_col[1], high_col[1], max_col[1])
  }
}

## summarise continuous and discrete traits -----------------------------------

### compute midpoints

### first midpoint of low&min and high&max values 
### continuous traits: compute midpoint, use single non-NA value if one is 
### missing, and round to 2 decimal places
midpoint_low <- function(low,min){
  round(
    ifelse(
      is.na(min),low,
      ifelse(is.na(low), min, (low + min) / 2)
    ),
    2
  )
}

midpoint_high <- function(high,max){
  round(
    ifelse(
      is.na(max),high,
      ifelse(is.na(high), max, (high + max) / 2)
    ),
    2
  )
}

### discrete traits: compute midpoint, use single non-NA value if one is 
### missing, and round to nearest integer
midpoint_discrete <- function(low, high) {
  round(
    ifelse(
      is.na(high), low,
      ifelse(is.na(low), high, (low + high) / 2)
    )
  )
}

midpoint_continuous <- function(low, high) {
  round(
    ifelse(
      is.na(high), low,
      ifelse(is.na(low), high, (low + high) / 2)
    ),
    2
  )
}


### identify all low, high, min and max columns
low_cols <- cols[str_detect(cols, "_low")]
min_cols <- cols[str_detect(cols,"_min")]
high_cols <- cols[str_detect(cols,"_high")]
max_cols <- cols[str_detect(cols,"_max")]

### identify continuous (length, width, diameter, etc.) vs 
### discrete (count, merosity)
continuous_pattern <- "(length|lenght|width|diameter|height)"
#discrete_pattern   <- "(count|merosity)"

### find base names for each type
continuous_traits <- str_remove(low_cols[str_detect(low_cols, continuous_pattern)], "_low.*")
#meristic_traits   <- str_remove(low_cols[str_detect(low_cols, discrete_pattern)], "_low.*")

### apply summarization
traits <- traits %>%
  {
    #low traits
    for (base in range_traits) {
      low_col  <- grep(paste0("^", base, "_low"), names(.), value = TRUE)
      min_col <- grep(paste0("^", base, "_min"), names(.), value = TRUE)
      if (length(low_col) == 1 && length(min_col) == 1) {
        .[[base]] <- midpoint_low(.[[low_col]], .[[min_col]])
      }
    }
        for (base in range_traits) {
          high_col  <- grep(paste0("^", base, "_high"), names(.), value = TRUE)
          max_col <- grep(paste0("^", base, "_max"), names(.), value = TRUE)
          if (length(high_col) == 1 && length(max_col) == 1) {
            .[[base]] <- midpoint_high(.[[high_col]], .[[max_col]])
          }
        }
        ### continuous traits
        for (base in continuous_traits) {
          low_col  <- grep(paste0("^", base, "_low"), names(.), value = TRUE)
          high_col <- grep(paste0("^", base, "_high"), names(.), value = TRUE)
          if (length(low_col) == 1 && length(high_col) == 1) {
            .[[base]] <- midpoint_continuous(.[[low_col]], .[[high_col]])
          }
        }
        .
      } %>%
  select(-matches("(_min|_low|_high|_max)"))

#now we have 43 variables 

## unit standardization -------------------------------------------------------

for (var in continuous_traits) {
  unit_col <- paste0(var, "_unit")
  
  ### make sure values are numeric
  traits[[var]] <- as.numeric(as.character(traits[[var]]))
  
  ### meters to cm
  idx_m <- !is.na(traits[[unit_col]]) & traits[[unit_col]] == "m"
  traits[[var]][idx_m] <- traits[[var]][idx_m] * 100
  
  ### decimeters to cm
  idx_dm <- !is.na(traits[[unit_col]]) & traits[[unit_col]] == "dm"
  traits[[var]][idx_dm] <- traits[[var]][idx_dm] * 10
  
  ### millimeters to cm
  idx_mm <- !is.na(traits[[unit_col]]) & traits[[unit_col]] == "mm"
  traits[[var]][idx_mm] <- traits[[var]][idx_mm] / 10
  
}

#now unit col can be deleted 
traits <- traits %>%
  select(-contains("_unit"))

#traits contain 31 variables
colnames(traits)

#checking traits with less than 15% completeness
traits_percent <- colMeans(!is.na(traits)) * 100
names(traits_percent[traits_percent < 15])
trait_completeness <- colMeans(!is.na(traits[,-1])) * 100 
#checking less than 30% 
names(traits_percent[traits_percent < 30])
traits_filtered <- traits[, c(TRUE, trait_completeness >= 30)]

#now i will keep 30% filtered columns and remove taxa with less than 60% completeness
species_completeness <- rowMeans(!is.na(traits_filtered[,-1])) * 100
traits_filtered$species_completeness <- species_completeness
write.csv(traits_filtered, "4.outputs/traits_filtered.csv", row.names = F)

#remove species less than 60% completeness
traits_final <- traits_filtered %>%
  dplyr::filter(species_completeness >= 60) %>%
  dplyr::select(-species_completeness)
#now we have 1773 species

write.csv(traits_final,"final_species_filtered_data.csv", row.names = F)

traits_final <- traits_final %>% 
  separate(taxon, into = c("genus", "epithet"),
           sep = "_")
unique(traits_final$genus)
traits_final %>% distinct() %>% count(genus, name = "traits_final")

