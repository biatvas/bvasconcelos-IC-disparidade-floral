#============#
# 1. Head ####
#============#

library(sf)
library(dplyr)

morpho_data <- read.csv("1.datasets/mimoseae_subset_clean.csv")
mimosoid_occ <- read.csv("1.datasets/MimosoidsFull-Jens_et_al_2023.csv") # Mimosoid occurence data available in Jens et al., 2023 (doi:)
eco_reg <- read_sf("1.datasets/Ecoregions2017/Ecoregions2017.shp") # Ecoregios shapefiles available in Dinerstein et al., 2017 (doi: 10.1093/biosci/bix014)

#===============================#
# 2. Occurence in ecoregions ####
#===============================#

#===========================#
## 2.1 Preparing dataset ####
#===========================#

# Species to extract occurence data
species <- gsub(pattern = "_"," ", morpho_data$taxon)

setdiff(species, unique(mimosoid_occ$species)) #57 species in our dataset not present in jens occurence data

# occurence data
# which(is.na(mimosoid_occ$decimalLatitude))
# which(is.na(mimosoid_occ$decimalLongitude))

occurrence <- mimosoid_occ 
# all(occurence$species %in% species)

# transforming occurence data in sf object and projecting
occurrence_points <- st_as_sf(
  occurence,
  coords = c("decimalLongitude", "decimalLatitude"),
  crs = 4326) %>% 
  st_transform(st_crs(eco_reg))

# checking if there are invalid geometries
table(st_is_valid(eco_reg)) 
which(!st_is_valid(eco_reg))
st_is_valid(eco_reg, reason = TRUE) %>%
  unique() #duplicated vertex

#correcting invalid geometries
eco_reg_valid <- st_make_valid(eco_reg)

#checking if it worked
table(st_is_valid(eco_reg_valid)) #no more invalid geometries

# occurrence_ecoregions <- st_join(
#   occurrence_points,
#   eco_reg_valid)

# occurrence_ecoregions_nearest <- st_nearest_feature(
#   occurrence_points,
#   eco_reg_valid)

intersects <- st_intersects(
  occurrence_points,
  eco_reg_valid)

test <- sapply(intersects, function(x) length(x) > 1) 
which(test == TRUE)

for (i in seq_along(intersects)) {
  if (length(intersects[[i]]) > 0) {
    occurrence$BIOME_NUM[i] <- eco_reg_valid$BIOME_NUM[intersects[[i]]]
    occurrence$BIOME_NAME[i] <- eco_reg_valid$BIOME_NAME[intersects[[i]]]
  } else {
    occurrence$BIOME_NUM[i] <- NA
    occurrence$BIOME_NAME[i] <- NA
  }
}

# selecting only mimosoid species that occur in morpho_data dataset
occurrence_ecoregions_subset <- occurrence_ecoregions[which(occurrence_ecoregions$species 
                                                            %in% species),]

# all(occurrence_ecoregions_subset$species %in% species)
# unique(occurrence_ecoregions_subset$species)

eco_per_species <- occurrence_ecoregions_subset %>% 
  st_drop_geometry() %>%
  count(species, BIOME_NAME, REALM, BIOME_NUM) %>%
  arrange(species, desc(n))

unique(eco_per_species$species)

eco_per_species %>% group_by(species) %>% 
  summarise(biome = paste(BIOME_NUM, collapse = "&"),
            .groups = "drop")

write_sf(occurrence_point, "3.outputs/occurence-data/occurrence_point")
write_sf(eco_reg_valid, "3.outputs/occurence-data/")

#dados florais, 
