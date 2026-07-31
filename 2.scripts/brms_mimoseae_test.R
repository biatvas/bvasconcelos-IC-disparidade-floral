#============================#
#==== molecular evolution rates ===
#================================#
library(ape)
library(brms)
library(ggplot2)
library(dplyr)
library(tibble)
library(posterior)
library(ggtree)
library(treeio)
library(circlize)
setwd("~/Documents/Labis/Dados/")

traits <- read.csv("~/Documents/Labis/Dados/1.datasets/moniq/continuous_data_cleaned-20260410.csv")
##read tree 系統樹読み込み
time_tree <- read.tree("3.trees/mimosoid_calibrated_clean_updated.tre")
substitution_tree <- read.tree("3.trees/mimosoid_branchesoptimized_clean_updated.tre")

all.equal.phylo(time_tree, substitution_tree,
                use.edge.length = FALSE)

#prune phylogeny 
tips_remove <- setdiff(time_tree$tip.label, traits$taxon)

#calculate rates
time_tree <- drop.tip(time_tree, tips_remove)
substitution_tree <- drop.tip(substitution_tree, tips_remove)
rates <- substitution_tree$edge.length / time_tree$edge.length
rates_df <- data.frame(rates = rates)

#分子進化速度のヒストグラムを作成/histogram 
ggplot(rates_df, aes(x = rates)) +
  geom_histogram(binwidth = 0.0005, fill = "#9f9f98", color = "black") +
  labs(title = "Distribution of Molecular Evolutionary Rates",
       x = "Molecular Evolutionary Rate",
       y = "Frequency") +
  theme_minimal()

#rates variam entre 0.0015 até 0.023

# 系統樹の色を決める selecting color for phylogeny
color_scale <- colorRamp2(c(0.0015, 0.0030, 0.010), c("#00a1e9", "#eae8e1", "#ea5532"))
rate_colors <- color_scale(rates)

# save plot
pdf("4.outputs/Mimoseae_molecular_evolutionary_rates.pdf", width = 7, height = 7)

plot.phylo(time_tree, type = "fan", edge.color = rate_colors, cex = 0.05, 
           edge.width = 0.5, main = "Rates of molecular evolution")

dev.off()

#============================#
#== brms =====================
#=============================#
tip_edges <- time_tree$edge[,2] <= length(time_tree$tip.label)

rates_df <- data.frame(
  taxon = time_tree$tip.label[time_tree$edge[tip_edges,2]],
  rate = rates[tip_edges]
)

#juntar o dataset com as taxas
mimoseae_data <- merge(traits, rates_df, by = "taxon")
