#===Molecular evolution analysis======
#=====2026.3.4======# 
setwd("C:/Users/Eulália/Desktop/bia/Labis/IC/Dados") #defining work directory

#extracting terminal branchs lengths 末端枝の長さを取り出すコード
library(ape)
tree <- ape:: read.tree("3.trees/mimosoid_calibrated_clean_updated.tre")
lengths <- tree$edge.length[tree$edge[,2] <= length(tree$tip.label)]

#==========================================================#
#==========分子進化速度を計算して、表を作る ===============
#=======calculating molecular evolutionary rates at terminal branches===#
time_tree <- read.tree("3.trees/mimosoid_calibrated_clean_updated.tre")
substitution_tree <- read.tree("3.trees/mimosoid_branchesoptimized_clean_updated.tre")

all.equal.phylo(time_tree, substitution_tree, use.edge.length = FALSE) #checking topology 

time_lengths <- time_tree$edge.length[time_tree$edge[,2] <= length(time_tree$tip.label)]
substitution_lengths <- substitution_tree$edge.length[substitution_tree$edge[,2] <= length(substitution_tree$tip.label)]

rate <- substitution_lengths / time_lengths

molecular_evolutionary_rate <- data.frame(
  Tip_in_tree = time_tree$tip.label,
  rate = rate
)

write.table(molecular_evolutionary_rate,
            file = "4.outputs/molecular_evolutionary_rates_mimoseae.txt",
            row.names = FALSE)

#checking rates distribution 
hist(molecular_evolutionary_rate$rate)

median(molecular_evolutionary_rate$rate)
#0.004034203
min(molecular_evolutionary_rate$rate)
#0.0003192668
max(molecular_evolutionary_rate$rate)
#0.02079727

molecular_evolutionary_rate %>%
  mutate(genus = sub("_.*", "", Tip_in_tree)) %>%
  group_by(genus) %>%
  summarise(mean_rate = mean(rate, na.rm = TRUE))
#to see median value for each genus 

#===2026.03.18====#
# 枝長の比率を計算====
library(ape)
library(ggplot2)
library(circlize)

time_tree <- read.tree("3.trees/mimosoid_calibrated_clean_updated.tre")
substitution_tree <- read.tree("3.trees/mimosoid_branchesoptimized_clean_updated.tre")

rates <- substitution_tree$edge.length / time_tree$edge.length
rates_df <- data.frame(rates = rates)

pdf("4.outputs/distribution_molecular_evolutionary_rates.pdf", width = 7, height = 7)

# ①分子進化速度のヒストグラムを作成/histogram
ggplot(rates_df, aes(x = rates)) +
  geom_histogram(binwidth = 0.0005, fill = "#9f9f98", color = "black") +
  labs(title = "Distribution of Molecular Evolutionary Rates",
       x = "Molecular Evolutionary Rate",
       y = "Frequency") +
  theme_minimal()
#first we plotted with width = 0.0001 but was too tiny 
#check for median, min and max for values 

# ②系統樹の色を決める
color_scale <- colorRamp2(c(0.001, 0.0025, 0.010), c("#00a1e9", "#eae8e1", "#ea5532"))
rate_colors <- color_scale(rates)

#plotting with different numbers
color_scale <- colorRamp2(c(0.001, 0.004, 0.015), c("#00a1e9", "#eae8e1", "#ea5532"))
rate_colors <- color_scale(rates)


# ③系統樹をプロット
plot.phylo(time_tree, type = "fan", align.tip.label = TRUE,
           edge.color = rate_colors, cex = 0.03, 
           edge.width = 0.8, main = "Molecular evolutionary rate")


# salvar o plot
pdf("4.outputs/molecular_evolutionary_rates_2.pdf", width = 7, height = 7)

plot.phylo(time_tree, type = "fan", align.tip.label = TRUE,
           edge.color = rate_colors, cex = 0.05, 
           edge.width = 0.5, main = "Molecular evolutionary rate")

dev.off()
#see tip names with high evolutionary rates 
threshold <- 0.01

red_edges <- which(rates >= threshold)
tip_indices <- time_tree$edge[red_edges, 2]
red_tips_idx <- tip_indices[tip_indices <= length(time_tree$tip.label)]
red_tips <- time_tree$tip.label[red_tips_idx]

#Testing for significance between rates?


