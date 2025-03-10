#!/bin/bash
#SBATCH --job-name=metabuli_db
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --time=10-23:00:00
#SBATCH --partition=gpu
#SBATCH --nodelist=devbox001
#SBATCH --output=/fast/lunajang/metabuli/exclusion_test/new_metabuli/logs/%j_output.log
#SBATCH --error=/fast/lunajang/metabuli/exclusion_test/new_metabuli/logs/%j_error.log

set -e

source ~/.bashrc
source activate base
conda activate metabuli

DEPTH=0.25
OUTPUT_DIR="/fast/lunajang/metabuli/exclusion_test/new_metabuli/fasta"
MASON2_SIMULATOR="/home/lunajang/src/mason2-2.0.9-Linux-x86_64_sse4/bin/mason_simulator"

mkdir -p "$OUTPUT_DIR/reads/single/depth_$DEPTH"
mkdir -p "$OUTPUT_DIR/reads/single/depth_$DEPTH/mason_result"
mkdir -p "$OUTPUT_DIR/temp_single_$DEPTH"

echo " Run Mason2 simulator"

awk -F '/' '{print $0, $NF }' "$OUTPUT_DIR/query.list" | while read -r FNA_FILE ACCESSION; do
    ACCESSION=$(basename "$ACCESSION" .fna)
    echo "Mason2 for $ACCESSION, FNA_FILE: $FNA_FILE"

    TEMP_FASTA="$OUTPUT_DIR/temp_single_$DEPTH/temp_$ACCESSION.fasta"
    seqkit seq -m 500 -g "$FNA_FILE" > "$TEMP_FASTA"

    if [ ! -s "$TEMP_FASTA" ]; then
        echo "Skipping $ACCESSION (no sequences ≥ 500bp)"
        rm -f "$TEMP_FASTA"
        continue
    fi

    GENOME_SIZE=$(wc -c < "$TEMP_FASTA")
    N=$(awk "BEGIN {print int((1/150) * $DEPTH * $GENOME_SIZE + 0.5)}")

    $MASON2_SIMULATOR \
        -q \
        --force-single-end \
        --illumina-read-length 150 \
        --illumina-prob-mismatch 0.0011 \
        --illumina-prob-mismatch-begin 0.00055 \
        --illumina-prob-mismatch-end 0.0022 \
        --fragment-mean-size 500 \
        --read-name-prefix "${ACCESSION}_" \
        -ir "$TEMP_FASTA" \
        -n "$N" \
        -o "$OUTPUT_DIR/reads/single/depth_$DEPTH/mason_result/$ACCESSION.fasta" \

    rm -f "$TEMP_FASTA"
    rm -f "$TEMP_FASTA.fai"
done

echo "Make query file and shuffle it"
cat $OUTPUT_DIR/reads/single/depth_$DEPTH/mason_result/*.fasta > "$OUTPUT_DIR/reads/single/depth_$DEPTH/query.fasta"
python "/fast/lunajang/metabuli/exclusion_test/Metabuli_DB/shuffle_fasta.py" \
    --query "$OUTPUT_DIR/reads/single/depth_$DEPTH/query.fasta" \
    --shuffled_query "$OUTPUT_DIR/reads/single/depth_$DEPTH/shuffled_query.fasta" 

echo ""

rm -rf "$OUTPUT_DIR/temp_single_$DEPTH"