# 環境設定
  library(ape)
  library(brms)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(posterior)
  library(ggtree)
  library(treeio)
  
setwd("~/Documents/Labis/Dados") 
  
#add rates in morphodataset ====
time_tree <- read.tree("3.trees/mimosoid_calibrated_clean_updated.tre")
substitution_tree <- read.tree("3.trees/mimosoid_branchesoptimized_clean_updated.tre")
  
all.equal.phylo(time_tree, substitution_tree, use.edge.length = FALSE) #checking topology 
  
time_lengths <- time_tree$edge.length[time_tree$edge[,2] <= length(time_tree$tip.label)]
substitution_lengths <- substitution_tree$edge.length[substitution_tree$edge[,2] <= length(substitution_tree$tip.label)]
  
rates <- substitution_lengths / time_lengths
rates_df <- data.frame(rates = rates)
  
molecular_evolutionary_rate <- data.frame(
    taxon = time_tree$tip.label,
    rate = rates
  )
molecular_evolutionary_rate$taxon <- gsub("_", " ", molecular_evolutionary_rate$taxon)
tree_table_filtered <- read.csv("4.outputs/morphoData_mimosa.csv")
tree_table_filtered <- merge(tree_table_filtered, molecular_evolutionary_rate, by = "taxon")

#0 are spike, 1 are raceme, 2 capitula and 3 umbel
#for the variable filament color, 0 is white, 1 is yellow, 2 is pink and 3 is red
#0 if they are diplo and 1 if they are haplo

set.seed(1234)
# データ読込
phylo_tree <- read.tree("3.trees/mimosoid_calibrated_clean_updated.tre")
phylo_tree$tip.label <- gsub("_", " ", phylo_tree$tip.label)
#because names in tree are Mimosa_myriadenia and in dataset are Mimosa myriadenia

#154 obs and 13 variables

#imputing data because are some NAs ====
tree_table <- tree_table_filtered %>%
  select(-nectary_presence, -CalShape)

#now we have 11 variables
get_moda <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# 1. integer
integer_cols <- names(tree_table)[sapply(tree_table, is.integer)]

tree_table[integer_cols] <- lapply(tree_table[integer_cols], function(x) {
  tidyr::replace_na(x, get_moda(x))
})

# 2. numéricas contínuas
numeric_cols <- names(tree_table)[sapply(tree_table, is.numeric)]
numeric_cols <- setdiff(numeric_cols, integer_cols)

tree_table[numeric_cols] <- lapply(tree_table[numeric_cols], function(x) {
  tidyr::replace_na(x, mean(x, na.rm = TRUE))
})

# 3. categóricas
categorical_cols <- names(tree_table)[sapply(tree_table, function(x) is.character(x) || is.factor(x))]

tree_table[categorical_cols] <- lapply(tree_table[categorical_cols], function(x) {
  tidyr::replace_na(x, get_moda(x))
})
#now all NA in dataset has a value 
# Na除去
tree_table <- tree_table[!is.na(tree_table$rate), ]
nrow(tree_table)
# 表を参照して系統樹をトリミング
phylo_tree_filtered <- drop.tip(
  phylo_tree,
  setdiff(phylo_tree$tip.label, tree_table$taxon)
)

# 相関行列 correlation matrix ====
vcv_matrix <- vcv(phylo_tree_filtered, corr = TRUE) 

# Aの並びをデータに合わせる
tip_order <- intersect(rownames(vcv_matrix), tree_table$taxon)
A <- vcv_matrix[tip_order, tip_order]
rownames(A) <- colnames(A) <- tip_order

cat("Is A symmetric? ", isSymmetric(A), "\n")
cat("Are all eigenvalues positive? ", all(eigen(A)$values > 0), "\n")

# データ側も Tip を A に一致させる
tree_table <- subset(tree_table, taxon %in% tip_order)
tree_table$taxon <- factor(tree_table$taxon, levels = tip_order)
# 分布確認 ====
p <- ggplot(tree_table, aes(x = rate)) +
  geom_histogram(bins = 40) +
  theme_minimal()
ggsave("rates_histogram.pdf", plot = p, width = 8, height = 6)

# brmsモデル（切片なし） without intercept ====
#using all traits
options(mc.cores = 16)
brms_model_all <- brm(
  formula = rate ~  stamenNumber + filament_color + InfType + FlowMerosity + inflorescence_peduncle_length_mean +
    inflorescence_length_mean + CalLengthMean + ColLengthMean +
    FilLengthMean + (1 | gr(taxon, cov = A)),  #bayesian PGLMM
  data   = tree_table,
  family = Gamma(link = "log"),
  data2  = list(A = A),
  prior  = c(
    prior(normal(0, 5), class = "b"),
    prior(exponential(1), class = "sd")
  ),
  control = list(adapt_delta = 0.95),
  iter = 5000, warmup = 2500, chains = 4, cores = 16, seed = 1234
)

#there were 52 divergent transition

sink("brms_mimosa_all.txt"); print(summary(brms_model_all)); sink()

# サマリー・事後サンプルの保存 ====
# summary and posterior predictive 
out <- list(
  fixef   = fixef(brms_model_all),
  ranef   = tryCatch(ranef(brms_model_all), error = function(e) NULL),
  draws   = posterior::as_draws_df(brms_model_all),
  formula = formula(brms_model_all),
  family  = family(brms_model_all)
)
saveRDS(out, "brmsmimosa_postproc_bundle.rds")

fixef_df <- fixef(brms_model_all) %>% as.data.frame() %>% rownames_to_column("Parameter")
write.csv(fixef_df, "brms_fixef_summary_testing_allmimosa.csv", row.names = FALSE)
saveRDS(fixef_df, "brms_fixef_summary_testing_allmimosa.rds")

plot(brms_model_all)

# ランダム効果 ==== 
pdf("random_effects_mimosa.pdf", width = 10, height = 6)
try({
  plot(ranef(brms_model_all)$taxon[, , "Intercept"],
       main = "Random Effects (taxon)")
})
dev.off()

# 残差 vs フィット ====
fitted_vals    <- fitted(brms_model_all)[, "Estimate"]
residuals_vals <- residuals(brms_model_all)[, "Estimate"]
pdf("residuals_vs_fitted_mimosa.pdf", width = 8, height = 6)
plot(fitted_vals, residuals_vals,
     xlab = "Fitted values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red")
dev.off()

#=== 事後予測チェック====
pdf("posterior_predictive_check_mimosa.pdf", width = 8, height = 6)
pp_check(brms_model_all)
dev.off()

#=== 固定効果の図示（95% CI）====
library(bayesplot)
color_scheme_set("blue")

pdf("fixed_effects_forest_mimosall.pdf", width = 8, height = 6)
mcmc_plot(brms_model_all, 
          variable = NULL, 
          prob = 0.95, 
          type = "areas") +
  ggtitle("Posterior distributions with 95% CI")
dev.off()

##testing other plot type 
library(tidyverse)

fixef_df <- fixef(brms_model_all) |>
  as.data.frame() |>
  rownames_to_column("Parameter") |>
  filter(Parameter != "Intercept") |>
  mutate(
    Parameter = recode(Parameter,
                       stamenNumber = "Stamen number",
                       filament_color = "Filament color",
                       InfType = "Inflorescence type",
                       FlowMerosity = "Flower merosity",
                       inflorescence_peduncle_length_mean = "Peduncle length",
                       inflorescence_length_mean = "Inflorescence length",
                       CalLengthMean = "Calyx length",
                       ColLengthMean = "Corolla length",
                       FilLengthMean = "Filament length"
    )
  )

##make a plot 
ggplot(fixef_df,
       aes(y = reorder(Parameter, Estimate))) +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             colour = "grey50") +
  geom_segment(aes(x = Q2.5,
                   xend = Q97.5,
                   yend = reorder(Parameter, Estimate)),
               linewidth = 0.8) +
  geom_point(aes(x = Estimate),
             size = 2.8) +
  labs(
    x = "Posterior estimate",
    y = NULL,
    title = "Effects of floral traits on molecular evolutionary rate"
  ) +
  theme_classic(base_size = 13)


