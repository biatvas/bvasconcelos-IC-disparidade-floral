setwd("~/Documents/Labis/BEPE")
library(ape)
library(phytools)
library(stringr)

#calibrated tree
tree <- read.tree("~/Documents/Labis/BEPE/2.trees/mimosoid_branchesoptimized_clean_updated.tre")

# Função para extrair "Genero_epiteto" do nome do fasta
# "Mimosoid_Mimosa_pudica_KIB_P025_WE02" → "Mimosa_pudica"
extrair_nome <- function(nome_fasta) {
  partes <- str_split(nome_fasta, "_")[[1]]
  if (length(partes) < 3) return(NA)
  paste(partes[2], partes[3], sep = "_")  # genero = 2a parte, epiteto = 3a
}

#listar todos os fastas
fastas <- list.files("data_peelab/5. trimmed/", 
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
  nomes_validos <- nomes_arvore[nomes_arvore %in% tree$tip.label]
  
  # Avisar se algum nome não foi encontrado
  ausentes <- nomes_arvore[!nomes_arvore %in% tree$tip.label]
  if (length(ausentes) > 0) {
    cat("Locus", nome_locus, "— tips não encontrados na árvore:\n")
    print(ausentes)
  }
  
  # Fazer o prune se tiver pelo menos 4 tips
  if (length(nomes_validos) >= 4) {
    tree_podada <- keep.tip(tree, nomes_validos)
    write.tree(tree_podada, 
               paste0("arvores_locus_mimoseae/", nome_locus, "_tree.nwk"))
    cat("Locus", nome_locus, "— árvore salva com", 
        length(nomes_validos), "tips\n")
  } else {
    cat("Locus", nome_locus, "— tips insuficientes, pulando\n")
  }
}

locus1tree <- read.tree("~/Documents/Labis/BEPE/arvores_locus_mimoseae/Locus1.codon_aln.trimmed_tree.nwk")
locus9tree <- read.tree("~/Documents/Labis/BEPE/arvores_locus_mimoseae/Locus9.codon_aln.trimmed_tree.nwk")
locus94tree <- read.tree("~/Documents/Labis/BEPE/arvores_locus_mimoseae/Locus94.codon_aln.trimmed_tree.nwk")
#o numero de especies na arvore e no alinhamento nao esta consistente. 
#posso inicialmente limpar o nome das especies no alinhamento, pra ficar somente genero + epiteto especifico e em seguida, 
#filtrar as sequencias somente pro que tiver na arvore
