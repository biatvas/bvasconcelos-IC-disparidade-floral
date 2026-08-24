#selecting traits for analysis 
morpho_data <- read.csv("1.datasets/species_features.csv", sep = ";") #reading morphological dataset
#==================#
# selecting traits ----
trait_selection <- read.csv2("1.datasets/trait_selection.csv", sep = ",")
trait_selection[which(trait_selection$keep == "1"), "trait"]
trait_keep <- subset(trait_selection, keep == 1)
unique(trait_keep$trait) 
morpho_data_selected <- select(morpho_data, c(unique(trait_keep$trait)))

#after cleaning the preliminar morphodata, I will combine with other two datasets 
write.csv(morpho_data_selected, "morpho_data_selected.csv", row.names = FALSE)
#manualmente, limpei a tabela do ryan e de mimosa, devido a inconsistência do nome de caracteres
