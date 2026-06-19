# 1. Head ----------------------------------------------------------------------
#==============================================================================#
setwd("C:/Users/Eulália/Desktop/bia/Labis/IC/Dados") #defining work directory
if (!require(librarian)) install.packages("librarian"); library("librarian")
librarian::shelf(phytools, dplyr, purrr, factoextra,vegan,tidyverse,ape,stringr,readr) #installing and/or loading packages
mimosa_species <- read.csv("1.datasets/tabelas/mimosa_species_features.csv")

#========================================================================#
#==2.1.1 clean species names w/ tnrs=====================================#
#========================================================================#
std_tip_labels <- function(tip_labels){
  for(i in 1:length(tip_labels)){
    # Replace underscore by space
    tip_labels[i] <- gsub(pattern = "_",
                          x = tip_labels[i],
                          replacement = " ")
    
    # Remove white space
    tip_labels[i] <- trimws(tip_labels[i]) 
    
    # Replace double space with single space
    tip_labels[i] <- gsub(pattern = "  ", 
                          x = tip_labels[i], 
                          replacement = " ") 
    
    # Remove accents
    tip_labels[i] <- stri_trans_general(str = tip_labels[i], 
                                        id = "Latin-ASCII")
  }
  return(tip_labels)
}

# Standardize tip labels for replacement in the phylogenetic tree
std_tip_labels_2 <- function(tip_labels){
  for(i in 1:length(tip_labels)){
    # Remove abbreviations (var.)
    tip_labels[i] <- gsub(pattern = "var.",
                          replacement = "",
                          x = tip_labels[i],
                          fixed = T)
    
    # Removing abbreviations (subsp.)
    tip_labels[i] <- gsub(pattern = "subsp.",
                          replacement = "",
                          x = tip_labels[i],
                          fixed = T)
    
    # Remove white spaces
    tip_labels[i] <- trimws(tip_labels[i]) 
    
    # Replace double spaces by single spaces
    tip_labels[i] <- gsub(pattern = "  ", 
                          x = tip_labels[i], 
                          replacement = " ") 
    
    # Replace spaces by underscores
    tip_labels[i] <- gsub(pattern = " ",
                          replacement = "_",
                          x = tip_labels[i],
                          fixed = T)
    
    # Remove hyphenation 
    tip_labels[i] <- gsub(pattern = "-",
                          replacement = "",
                          x = tip_labels[i],
                          fixed = T)
  }
  return(tip_labels)
}
rm_duplicate_tips <- function(tip_labels, tree){
  # Create a list of tip names associated with their correspondent index/indexes
  ## If a name is not unique, then it will be associated with two or more 
  ## indexes
  tips_list <- tip_labels %>%
    as.data.frame() %>%
    rownames_to_column("number") %>%
    split(.$., .$number)
  
  # Choose unique tips randomly
  set.seed(7)
  unique_tips <- sapply(tips_list, function(df) sample(df$number, 1)) %>%
    as.numeric()
  # Keep unique tips
  keep.tip(species, unique_tips)
}

# Collapse infraspecific names into species-level names
collapse_infra <- function(tip_labels){
  for(i in 1:length(tip_labels)){
    # Split epithets 
    split_tip_labels <- str_split(string = tip_labels[i],
                                  pattern = "_")[[1]]
    
    # Reassign tip label using only generic and specific epithets
    tip_labels[i] <- paste(split_tip_labels[1], split_tip_labels[2],
                           sep = "_")
  }
  return(tip_labels)
}

## Clean tip labels with TNRS -----------------------------------------------

# Extract tips 
tips <- mimosa_species$taxon

# Standardize tip names for TNRS queries
tips <- std_tip_labels(tips)
#collapse infra specific names

# Run TNRS using WCVP as source
tnrs <- TNRS(tips, sources = "wcvp", classification = "wfo", mode = "resolve", matches = "best")
setdiff(tnrs$Name_submitted, tips)
which(table(tips)>1, T)

write.csv(tnrs,"4.outputs/mimosa/mimosa_queries.csv", row.names = F)
# Read TNRS queries
tnrs <- read.csv("4.outputs/mimosa/mimosa_queries.csv",na.strings = c("", NA), stringsAsFactors = F, encoding = "UTF-8")
# Filter names with an overall score (matching index) equal to 1 (exact match),
# and with all possible taxonomic statuses except "No opinion"
resolved <- tnrs %>% filter(Overall_score == 1 &
                              Taxonomic_status != "No opinion") %>%
  select(Name_submitted,
         Overall_score,
         Name_matched,
         Accepted_name)
fuzzy <- tnrs %>% filter(Overall_score < 1 |
                           Taxonomic_status == "No opinion") %>%
  select(Name_submitted,
         Overall_score,
         Name_matched,
         Accepted_name)  
#103 nomes para corrigir
# Write *.csv for manual checking of fuzzy matches
## Check `clean_occurrence.R` for details on this procedure
write.csv(cbind(fuzzy, data.frame(Keep = NA, Altered = NA)), "4.outputs/mimosa/mimosa_fuzzy_checked.csv",row.names = FALSE)
# Read *.csv after manual checking
fuzzy_checked <- read.csv("4.outputs/mimosa/mimosa_fuzzy_checked.csv", na.strings = c("", NA), stringsAsFactors = F,encoding = "UTF-8")
# Merge automatically and manually checked taxonomic names
tips_clean <- rbind(cbind(resolved, Keep = 1, Altered = NA), fuzzy_checked)

# Replace original tips and flag names to be removed
for(i in 1:length(tips)){
  if(tips_clean$Keep[tips_clean$Name_submitted == tips[i]] == 0){
    tips[i] <- "drop"
  } else{
    tips[i] <- 
      tips_clean$Accepted_name[tips_clean$Name_submitted == tips[i]]
  }
}
#como deu erro talvez pela ausencia de alguns nomes, vou rodar um script adaptado 
for(i in seq_along(tips)) {
  
  keep_val <- tips_clean$Keep[tips_clean$Name_submitted == tips[i]]
  
  if(length(keep_val) == 0) {
    next   # or tips[i] <- NA  (your choice)
  }
  
  if(keep_val[1] == 0) {
    tips[i] <- "drop"
  } else {
    tips[i] <- tips_clean$Accepted_name[
      tips_clean$Name_submitted == tips[i]
    ][1]
  }
  
}
# Standardize tips before placing them back into the phylogenetic tree
tips <- std_tip_labels_2(tips)
# Place clean tips back into the phylogenetic tree
mimosa_species$species <- tips
# Drop flagged tips
mimosa_species <- filter(mimosa_species, which(mimosa_species$species == "drop"))
mimosa_species <- mimosea_species |>
  dplyr::filter(species != "drop")

# Remove duplicate tips ---------------------------------------------
## 
# Which are the duplicate tips? 
#collapsed_species$species[duplicated(collapsed_species$species)]
# Drop duplicate tips
# <- rm_duplicate_tips(tip_labels = collapsed_species$species,
# tree = collapsed_species
# Save new data
write.csv(mimosa_species, "4.outputs/mimosa/mimosadata_clean.csv", row.names = FALSE)
