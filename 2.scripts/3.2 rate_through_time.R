## Plot rates of molecular evolution through time 
#rate through time plot 

#Load libraries 
library(ape)
library(tidyverse)
library(phytools)

#get ready
setwd("~/Documents/Labis/Dados")
#cada ramo tem uma taxa e buscamos ver como isso muda ao longo do tempo 

#ler arvores
time_tree <- read.tree("3.trees/mimosoid_calibrated_clean_updated.tre")
substitution_tree <- read.tree("3.trees/mimosoid_branchesoptimized_clean_updated.tre")

# garantir mesma ordem
time_tree <- reorder.phylo(time_tree, "postorder")
substitution_tree <- reorder.phylo(substitution_tree, "postorder")

#check topology
all.equal.phylo(time_tree, substitution_tree,
                use.edge.length = FALSE)
#table for edge length
time_df <- data.frame(
  parent = time_tree$edge[,1],
  child = time_tree$edge[,2],
  time_length = time_tree$edge.length
)

sub_df <- data.frame(
  parent = substitution_tree$edge[,1],
  child = substitution_tree$edge[,2],
  subst_length = substitution_tree$edge.length
)

edge_df <- merge(time_df, sub_df, by = c("parent", "child"))

#rate for each branch 
edge_df$rate <- edge_df$subst_length / edge_df$time_length #taxa instantanea media do ramo 

#extract node ages 
node_ages <- branching.times(time_tree)

#check names 
node_ages
names(node_ages)
#rate means number of substitution per time 
n_tips <- length(time_tree$tip.label)

#ages <- numeric(n_tips + time_tree$Nnode)
#ages[as.numeric(names(node_ages))] <- node_ages
max_node <- max(as.numeric(names(node_ages)))
ages <- numeric(max_node)
ages[as.numeric(names(node_ages))] <- node_ages

# Build edge table with start/end times
edges <- time_tree$edge
edge_df$start <- ages[edges[,1]]
edge_df$end <- ages[edges[,2]]
#aqui cada ramo vai receber uma idade de quando surgiu 

# Ensure that the starting age of a branch is never lower than its ending age
all(edge_df$start >= edge_df$end)

# Define time bins
max_age <- max(edge_df$start, na.rm = TRUE)
bin_width <- 2 #for 2 Mya 
time_breaks <- seq(0, ceiling(max_age / bin_width) * bin_width, 
                   by = bin_width)

bin_df <- data.frame(
  bin_start = time_breaks[-length(time_breaks)],
  bin_end   = time_breaks[-1]
)

#compute overlap per bin 
#overlap vai definir o tempo do ramo dentro do bin (intervalo de tempo)
slice_overlap <- function(start, end, bin_start, bin_end){
  max(0,min(start,bin_end) - max(end, bin_start))
}

#divides events by branch length 
results <- data.frame()

for(i in 1:nrow(bin_df)) {
  bin_start <- bin_df$bin_start[i]
  bin_end <- bin_df$bin_end[i]
  
  edge_df$overlap <- mapply(
    slice_overlap,
    edge_df$start,
    edge_df$end,
    MoreArgs = list(bin_start, bin_end)
  )
  tmp <- subset(edge_df, overlap > 0)
  
#taxa media e tempo de contribuicao 
#aqui faz uma correção extraindo media ponderada
bin_rate <- sum(tmp$rate * tmp$overlap)/ sum(tmp$overlap)
  
  results <- rbind(results, data.frame(
    bin_start = bin_start,
    bin_end = bin_end, 
    rate = bin_rate
  ))
}
# ou seja, o que a gente ta extraindo é a taxa média dos ramos
#que existiam naquele tempo 

#plot 
rate_graph <- ggplot(results, aes(x = bin_end, y = rate)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_x_reverse() +
  theme_classic() +
  labs(
    x = "Time before present (Mya)",
    y = "Mean substitution rate"
  )

# salvar o plot
pdf("4.outputs/molecular_evolutionary_rates_in_time.pdf", width = 7, height = 7)

ggplot(results, aes(x = bin_end, y = rate)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_x_reverse() +
  theme_classic() +
  labs(
    x = "Time before present (Mya)",
    y = "Mean substitution rate"
  )
dev.off()

#plot tree with rates 
rates <- edge_df$rate
names(rates) <- 1:nrow(time_tree$edge)

summary(rates)

#plot
library(ggtree)
library(ggplot2)
library(patchwork)

#check number of bins 
n_branches <- nrow(tmp)

#boxplot of each rate for bin
all_rates <- data.frame()

for(i in 1:nrow(bin_df)) {
  
  bin_start <- bin_df$bin_start[i]
  bin_end <- bin_df$bin_end[i]
  
  edge_df$overlap <- mapply(
    slice_overlap,
    edge_df$start,
    edge_df$end,
    MoreArgs = list(bin_start, bin_end)
  )
  
  tmp <- subset(edge_df, overlap > 0)
  
  tmp$bin_time <- bin_end
  
  all_rates <- rbind(all_rates, tmp)
}

#boxplot
p_box <- ggplot(all_rates,
                aes(x = factor(bin_time), y = rate)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_classic() +
  labs(x = "Time bin (Mya)",
       y = "Substitution rate") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_time <- ggplot(results, aes(x = bin_end, y = rate)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_x_reverse() +
  theme_classic() +
  labs(x = "Time before present (Mya)",
       y = "Mean substitution rate")

#plot tree and rates 
rates <- substitution_tree$edge.length / time_tree$edge.length
rates_df <- data.frame(rates = rates)

library(circlize)
color_scale <- colorRamp2(c(0.001, 0.004, 0.015), c("#00a1e9", "#eae8e1", "#ea5532"))
rate_colors <- color_scale(rates)

#plot tree
tree_rates <- plot.phylo(time_tree, align.tip.label = TRUE,
                         edge.color = rate_colors, cex = 0.03, 
                         edge.width = 0.8, main = "Molecular evolutionary rate")

#plot both analysis toghether 
# criar dataframe das arestas
tree_data <- data.frame(
  node = time_tree$edge[,2],
  rate = rates
)

# construir árvore
p_tree <- ggtree(time_tree)

# juntar dados corretamente
p_tree <- p_tree %<+% tree_data +
  geom_tree(aes(color = rate), linewidth = 0.5) +
  scale_color_gradient(
    low = "#00a1e9",
    high = "#ea5532"
  ) +
  theme_tree2()

# gráfico temporal
p_time <- ggplot(results, aes(x = bin_end, y = rate)) +
  geom_line() +
  geom_point() +
  scale_x_reverse() +
  theme_classic()

# combinar
combined_plot <- cowplot::plot_grid(
  p_tree,
  p_time,
  ncol = 1,
  align = "v",
  rel_heights = c(2,1)
)

combined_plot
