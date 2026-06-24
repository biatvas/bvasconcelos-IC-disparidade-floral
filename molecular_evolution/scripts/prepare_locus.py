#!/usr/bin/env python3
# prepare_locus.py
# Usage: python3 prepare_locus.py <input.fasta> <full_tree.nwk> <out.fasta> <out.nwk>

import sys
from Bio import SeqIO
from ete3 import Tree

def shorten_name(full_name):
    """
    'Mimosoid_Acacia_sabulosa_FMN_P038_WB04' → 'Acacia_sabulosa'
    Keeps fields at index [1] and [2] (genus + species epithet).
    """
    parts = full_name.lstrip(">").split("_")
    return "_".join(full_name.lstrip(">").split("_")[1:3])

input_fasta, full_tree, out_fasta, out_tree = sys.argv[1:]

# ── 1. Read alignment and rename sequences ──────────────────────────────
records = []
short_names = []

for rec in SeqIO.parse(input_fasta, "fasta"):
    short = shorten_name(rec.id)
    rec.id = short
    rec.name = short
    rec.description = ""
    records.append(rec)
    short_names.append(short)

SeqIO.write(records, out_fasta, "fasta")

# ── 2. Prune tree to taxa present in this alignment ─────────────────────
t = Tree(full_tree)

tree_leaves = {leaf.name for leaf in t.get_leaves()}
to_keep = [n for n in short_names if n in tree_leaves]

missing_in_tree = set(short_names) - tree_leaves
if missing_in_tree:
    print(f"  [WARN] Not in tree, will be skipped: {missing_in_tree}", file=sys.stderr)

if len(to_keep) < 4:
    print(f"  [FAIL] Too few taxa after pruning ({len(to_keep)}), skipping.", file=sys.stderr)
    sys.exit(1)

t.prune(to_keep, preserve_branch_length=True)
t.write(outfile=out_tree, format=1)

print(f"  [OK] {len(to_keep)} taxa kept", file=sys.stderr)


