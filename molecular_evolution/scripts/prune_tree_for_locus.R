library(ape)
library(phytools)

#calibrated tree
tree <- read.tree("~/Documents/Labis/BEPE/2.trees/mimosoid_branchesoptimized_clean_updated.tre")

#check mimosa tips
mimosa_tips <- tree$tip.label[grep("Mimosa_", tree$tip.label)]
cat("Número de tips Mimosa encontrados:", length(mimosa_tips), "\n")
print(mimosa_tips)

# Extrair o subtree
mimosa_tree <- keep.tip(tree, mimosa_tips)
# Checar se ficou correto
plot(mimosa_tree, cex = 0.4)
# Salvar
write.tree(mimosa_tree, "~/Documents/Labis/BEPE/2.trees/mimosa_subtree.nwk")

#Now, we need to prune tree for each locus 
library(ape)
library(stringr)

mimosa_tree <- read.tree("2.trees/mimosa_subtree.nwk")

# Função para extrair "Genero_epiteto" do nome do fasta
# "Mimosoid_Mimosa_pudica_KIB_P025_WE02" → "Mimosa_pudica"
extrair_nome <- function(nome_fasta) {
  partes <- str_split(nome_fasta, "_")[[1]]
  # Encontra a posição de "Mimosa" (ou outro gênero)
  idx <- which(partes == "Mimosa")
  if (length(idx) == 0) return(NA)
  paste(partes[idx], partes[idx + 1], sep = "_")
}

#listar todos os fastas
fastas <- list.files("data_peelab/Mimosa-analysis/5.mimosa_trimal/", 
                     pattern = "*.fasta", 
                     full.names = TRUE)

for (fasta_path in fastas) {
  # Nome do locus para salvar os arquivos
  nome_locus <- tools::file_path_sans_ext(basename(fasta_path))
  
  # Ler os nomes das sequências do fasta
  linhas <- readLines(fasta_path)
  headers <- linhas[grep("^>", linhas)]
  nomes_fasta <- gsub("^>", "", headers)
  
  # Converter para formato da árvore
  nomes_arvore <- sapply(nomes_fasta, extrair_nome)
  nomes_arvore <- nomes_arvore[!is.na(nomes_arvore)]
  
  # Checar quais nomes existem na árvore
  nomes_validos <- nomes_arvore[nomes_arvore %in% mimosa_tree$tip.label]
  
  # Avisar se algum nome não foi encontrado
  ausentes <- nomes_arvore[!nomes_arvore %in% mimosa_tree$tip.label]
  if (length(ausentes) > 0) {
    cat("Locus", nome_locus, "— tips não encontrados na árvore:\n")
    print(ausentes)
  }
  
  # Fazer o prune se tiver pelo menos 4 tips
  if (length(nomes_validos) >= 4) {
    tree_podada <- keep.tip(mimosa_tree, nomes_validos)
    write.tree(tree_podada, 
               paste0("arvores_por_locus_trimmed/", nome_locus, "_tree.nwk"))
    cat("Locus", nome_locus, "— árvore salva com", 
        length(nomes_validos), "tips\n")
  } else {
    cat("Locus", nome_locus, "— tips insuficientes, pulando\n")
  }
}

locus9tree <-read.tree("~/Documents/Labis/BEPE/arvores_por_locus_trimmed/Locus9.mimosa.codon_aln.trimmed_tree.nwk")