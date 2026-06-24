#test hyphy for one locus 
cd Documents/Labis/Dados/PEE_data/analysis

python3 prepare_locus.py 1.filtered/Locus1.fasta trees/mimosoid_branchesoptimized_clean_updated.tre Locus1clean.fasta pruned.nwk 
mkdir prunedtrees
mkdir -p cleanfastas

for f in 1.filtered/*.fasta; do
    locus=$(basename "$f" .fasta)
    python3 prepare_locus.py "$f" trees/mimosoid_branchesoptimized_clean_updated.tre "cleanfastas/${locus}clean.fasta" "prunedtrees/${locus}_pruned.nwk"
done

nano run_hyphy.sh
#!/bin/bash
# ── paths ──────────────────────────────────────────
CLEAN_DIR="cleanfastas"        # output from prepare_locus.py
TREE_DIR="prunedtrees"         # pruned trees from prepare_locus.py
OUT_DIR="hyphy_results"
HYPHY="hyphy"
# ────────────────────────────────────────────────────────

mkdir -p "$OUT_DIR"

for FASTA in "$CLEAN_DIR"/*clean.fasta; do

    LOCUS=$(basename "$FASTA" clean.fasta)
    PRUNED_TREE="$TREE_DIR/${LOCUS}_pruned.nwk"
    OUTFILE="$OUT_DIR/${LOCUS}.json"

    # check pruned tree exists for this locus
    if [[ ! -f "$PRUNED_TREE" ]]; then
        echo "[SKIP] $LOCUS — no pruned tree found"
        continue
    fi

    if [[ -f "$OUTFILE" ]]; then
        echo "[SKIP] $LOCUS — already done"
        continue
    fi

    echo "[RUN] $LOCUS"

    "$HYPHY" aBSREL \
        --alignment "$FASTA" \
        --tree "$PRUNED_TREE" \
        --type local \
        --output "$OUTFILE" \
        2> "$OUT_DIR/${LOCUS}.log"

    if [[ $? -eq 0 ]]; then
        echo "[OK]  $LOCUS"
    else
        echo "[FAIL] $LOCUS — check $OUT_DIR/${LOCUS}.log"
    fi

done