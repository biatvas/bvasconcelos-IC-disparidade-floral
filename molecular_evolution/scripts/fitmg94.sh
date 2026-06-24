cd /Users/biatvas/Documents/Labis/BEPE/

python scripts/rename_and_filter.py \
        --fasta  data_peelab/5.trimmed \
        --trees  arvores_locus_mimoseae \
        --out6   data_peelab/6.codon_name_clean \
        --out7   data_peelab/7.match_tree_codon

# Para um locus de teste
hyphy /Users/biatvas/Documents/Labis/BEPE/data_peelab/Mimosa-analysis/hyphy-analyses/FitMG94/FitMG94.bf \
  --alignment data_peelab/7.match_tree_codon/Locus1.codon_aln.trimmed.fasta \
  --tree arvores_locus_mimoseae/Locus1.codon_aln.trimmed_tree.nwk\
  --output data_peelab/8.FitMGresults/mimoseae_locus1_resultado.json \
  --type local 

#para rodar em loop
for locus in 1 2 3; do
  hyphy /Users/biatvas/Documents/Labis/BEPE/data_peelab/Mimosa-analysis/hyphy-analyses/FitMG94/FitMG94.bf \
    --alignment data_peelab/7.match_tree_codon/Locus${locus}.codon_aln.trimmed.fasta \
    --tree arvores_locus_mimoseae/Locus${locus}.codon_aln.trimmed_tree.nwk \
    --output data_peelab/8.FitMGresults/mimoseae_locus${locus}_resultado.json \
    --type local \
done



