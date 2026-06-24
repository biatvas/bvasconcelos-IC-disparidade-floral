#!/usr/bin/env python3
"""
Pipeline de limpeza de nomes e filtragem de alinhamentos por árvore.

Etapa 1 (pasta 6): Renomeia sequências para Genus_species
Etapa 2 (pasta 7): Remove sequências ausentes na árvore correspondente

Script de exemplo:
    python rename_and_filter.py \
        --fasta  data_peelab/5.trimmed \
        --trees  data_peelab/trees \
        --out6   data_peelab/6.codon_name_clean \
        --out7   data_peelab/7.match_tree_codon
"""

import argparse
import re
import sys
from pathlib import Path

from Bio import SeqIO

# ---------------------------------------------------------------------------
# Utilitários
# ---------------------------------------------------------------------------

def shorten_name(seq_id: str) -> str:
    """
    Extrai Genus_species de um ID no formato:
        Clade_Genus_species_VOUCHER_...
    ou qualquer variante com underscores.

    Estratégia: pega os dois primeiros tokens separados por '_' que
    começam com letra maiúscula (gênero) e minúscula (epíteto).
    Se não achar esse padrão, devolve os dois primeiros tokens simples.
    """
    parts = seq_id.split("_")

    #Genus (maiúscula) + species (minúscula)
    for i in range(len(parts) - 1):
        if parts[i][0].isupper() and parts[i + 1][0].islower():
            return f"{parts[i]}_{parts[i + 1]}"

    # Fallback: dois primeiros tokens
    if len(parts) >= 2:
        return f"{parts[0]}_{parts[1]}"
    return parts[0]


def taxa_in_tree(tree_file: Path) -> set:
    """
    Lê os nomes de táxons de um arquivo Newick.
    """
    text = tree_file.read_text()

    # Remove comentários entre colchetes
    text = re.sub(r"\[.*?\]", "", text, flags=re.DOTALL)

    # Extrai labels: qualquer token alfanumérico/_/. antes de : ou , ou ) ou (
    taxa = set(re.findall(r"([A-Za-z][A-Za-z0-9_.]+)(?=\s*[,):;])", text))
    return taxa


# ---------------------------------------------------------------------------
# Etapa 1 – renomear
# ---------------------------------------------------------------------------

def step1_rename(fasta_dir: Path, out_dir: Path):
    out_dir.mkdir(parents=True, exist_ok=True)
    fastas = sorted(fasta_dir.glob("*.fasta")) + sorted(fasta_dir.glob("*.fa"))

    if not fastas:
        print(f"[AVISO] Nenhum .fasta/.fa encontrado em {fasta_dir}", file=sys.stderr)
        return

    for fasta in fastas:
        records = []
        seen = {}  # para detectar duplicatas após renomeação

        for rec in SeqIO.parse(fasta, "fasta"):
            short = shorten_name(rec.id)

            # Se o nome curto já existir, adiciona sufixo numérico
            if short in seen:
                seen[short] += 1
                short = f"{short}_{seen[short]}"
            else:
                seen[short] = 0

            rec.id = short
            rec.name = short
            rec.description = ""
            records.append(rec)

        out_file = out_dir / fasta.name
        SeqIO.write(records, out_file, "fasta")
        print(f"[1] {fasta.name}: {len(records)} seqs → {out_file.name}")


# ---------------------------------------------------------------------------
# Etapa 2 – filtrar por árvore
# ---------------------------------------------------------------------------

def find_tree(gene: str, tree_dir: Path) -> Path | None:
    """
    Compara arquivo de alinhamento com a arvore 
    Aceita extensões .tre, .tree, .nwk, .nex.
    """
    exts = [".tre", ".tree", ".nwk", ".nex", ".treefile"]
    for ext in exts:
        candidate = tree_dir / (gene + ext)
        if candidate.exists():
            return candidate

    # Busca parcial: árvore cujo nome contém o gene
    for tree_file in tree_dir.iterdir():
        if tree_file.suffix in exts and gene in tree_file.stem:
            return tree_file

    return None


def step2_filter(fasta_dir: Path, tree_dir: Path, out_dir: Path):
    out_dir.mkdir(parents=True, exist_ok=True)
    fastas = sorted(fasta_dir.glob("*.fasta")) + sorted(fasta_dir.glob("*.fa"))

    if not fastas:
        print(f"[AVISO] Nenhum .fasta/.fa encontrado em {fasta_dir}", file=sys.stderr)
        return

    for fasta in fastas:
        gene = fasta.stem  # nome sem extensão usado para buscar a árvore
        tree_file = find_tree(gene, tree_dir)

        if tree_file is None:
            print(f"[AVISO] Árvore não encontrada para {fasta.name} — pulando.", file=sys.stderr)
            continue

        tree_taxa = taxa_in_tree(tree_file)
        kept, dropped = [], []

        for rec in SeqIO.parse(fasta, "fasta"):
            if rec.id in tree_taxa:
                kept.append(rec)
            else:
                dropped.append(rec.id)

        out_file = out_dir / fasta.name
        SeqIO.write(kept, out_file, "fasta")

        status = f"{len(kept)} mantidas"
        if dropped:
            status += f", {len(dropped)} removidas: {', '.join(dropped)}"
        print(f"[2] {fasta.name} (árvore: {tree_file.name}): {status}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fasta",  required=True, help="Pasta de entrada com os .fasta (etapa 1)")
    p.add_argument("--trees",  required=True, help="Pasta com as árvores podadas")
    p.add_argument("--out6",   required=True, help="Saída etapa 1: nomes limpos")
    p.add_argument("--out7",   required=True, help="Saída etapa 2: filtrado por árvore")
    p.add_argument("--skip1",  action="store_true",
                   help="Pula etapa 1 (usa --out6 como entrada da etapa 2)")
    return p.parse_args()


def main():
    args = parse_args()

    fasta_dir = Path(args.fasta)
    tree_dir  = Path(args.trees)
    out6      = Path(args.out6)
    out7      = Path(args.out7)

    if not args.skip1:
        print("\n=== Etapa 1: renomeando sequências ===")
        step1_rename(fasta_dir, out6)
    else:
        print("[INFO] Etapa 1 pulada. Usando", out6, "como entrada.")

    print("\n=== Etapa 2: filtrando por árvore ===")
    step2_filter(out6, tree_dir, out7)

    print("\nConcluído.")


if __name__ == "__main__":
    main()