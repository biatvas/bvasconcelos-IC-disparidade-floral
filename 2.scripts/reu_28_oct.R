# Phylomorphospace com hill smith
# Plot da contribuicao das variaveis

# Plotar os estados coloridos com hill smith

# Fazer de disparidade (sv e sr) usando hill smith para os clados
# e para os biomas

# violin plots de sv para biomas e para os clados

# gerar tabelas com os resultados (obs e median dos boots)

# gerar tabelas com os resultados de p-valor do Wilcoxon e 
# bhattacharyya coef para comparacao pareada dos clados

# gerar tabelas com os resultados de p-valor do Wilcoxon e 
# bhattacharyya coef para comparacao pareada dos biomas

# exemplos do violin pra mimosa e adenanthera
sr_clade$disparity$`Adenanthera clade`[[3]] #dados rarefeiros para usar no violin plot
sr_clade$disparity$`Cojoba clade`[[2]] #dados observados pois n é igual ao numero de rarefacao

sum <- as.vector(sv_clade$disparity$`Adenanthera clade`[[3]])
C1 <- c(rep("Adenanthera", bootstrap))

sv <- data.frame(C1, sum)
colnames(sv) <- c("Group", "Sum_of_variance")

sum_2 <- as.vector(sv_clade$disparity$`Mimosa clade`[[3]])
C2 <- c(rep("Mimosa", bootstrap))

sv_2 <- data.frame(C2, sum_2)
colnames(sv_2) <- c("Group", "Sum_of_variance")

clades_to_plot<-rbind(sv,sv_2)
ylab_list<-c(rep("Sum of variances",2))

ci <- function(x)
{
  exit <- quantile(x, probs = c(0.025,0.975))
  names(exit) <- c("ymin","ymax")
  return(exit) 
}

median <- function(x)
{
  exit <- quantile(x, probs = c(0.5))
  names(exit) <- c("y")
  return(exit) 
}

clades_to_plot$Group<-factor(clades_to_plot$Group,levels=c("Adenanthera","Mimosa"))
# fams_plot_disp[[i]]<-
ggplot(clades_to_plot, aes(x=Group, y=Sum_of_variance, 
                           colour=Group, fill=Group)) +
  geom_violin(alpha=0.5, show.legend=T) + labs(x="", y="") +
  xlab("") +
  ylab(ylab_list) +
  scale_fill_manual(values=c("dodgerblue3","gold2")) +
  scale_colour_manual(values=c("dodgerblue3","gold2")) +
  stat_summary(fun=ci, geom="line", size=0.6) +
  stat_summary(fun=median, geom="point", size=2.5) +
  theme_bw()  +
  theme(legend.title=element_blank(), legend.position="bottom", 
        legend.text=element_text(size=12), plot.title=element_text(size=14,face= "bold",hjust=0.5),
        axis.text.x=element_blank(), axis.text.y=element_text(size=10), axis.ticks.x=element_blank(), axis.title.x=element_text(size=12),axis.title.y=element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 



## só testando com o will smith ignorar
library(dispRity)
sub_clade <- custom.subsets(gower, groups_clade, dist.data = TRUE) #dataset for mpd 
sub_clade_hill <- custom.subsets(hs_full$li, groups_clade) #dataset for sv and sr

sub_biomes_hill <- custom.subsets(hs_full$li, groups_biome)
sub_biomes <- custom.subsets(gower_biomes, groups_biome, dist.data = TRUE)

bootstrapped_clades <- boot.matrix(sub_clade_hill, 
                                   bootstraps = bootstrap, 
                                   rarefaction = rarefaction)
sv_clade <- dispRity(bootstrapped_clades, metric = c(sum, variances))
sv_clade$disparity$`Adenanthera clade` #valores de sv do bootstrap
sv_clade$disparity$`Albizia clade` #valores de sv do bootstrap
clade_sv_smry <- summary(sv_clade)

clade_sv_wilcox <- test.dispRity(sv_clade, 
                                 wilcox.test,
                                 correction = "holm")

clade_sv_bhatt <- test.dispRity(sv_clade, bhatt.coeff)

bootstrapped_biomes <- boot.matrix(sub_biomes_hill, bootstraps = bootstrap, 
                                   rarefaction = rarefaction)



