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

DEPTH=8  
OUTPUT_DIR="/fast/lunajang/metabuli/exclusion_test/new_metabuli/fasta"
MASON2_SIMULATOR="/home/lunajang/src/mason2-2.0.9-Linux-x86_64_sse4/bin/mason_simulator"

mkdir -p "$OUTPUT_DIR/reads/paired/depth_$DEPTH"
mkdir -p "$OUTPUT_DIR/reads/paired/depth_$DEPTH/mason_result"
mkdir -p "$OUTPUT_DIR/temp_paired_$DEPTH"

echo " Run Mason2 simulator"

awk -F '/' '{print $0, $NF }' "$OUTPUT_DIR/query.list" | while read -r FNA_FILE ACCESSION; do
    ACCESSION=$(basename "$ACCESSION" .fna)
    echo "Mason2 for $ACCESSION, FNA_FILE: $FNA_FILE"

    TEMP_FASTA="$OUTPUT_DIR/temp_paired_$DEPTH/temp_$ACCESSION.fasta"
    seqkit seq -m 500 -g "$FNA_FILE" > "$TEMP_FASTA"

    if [ ! -s "$TEMP_FASTA" ]; then
        echo "Skipping $ACCESSION (no sequences ≥ 500bp)"
        rm -f "$TEMP_FASTA"
        continue
    fi

    GENOME_SIZE=$(wc -c < "$TEMP_FASTA")
    N=$(awk "BEGIN {print int((1/300) * $DEPTH * $GENOME_SIZE + 0.5)}")

    $MASON2_SIMULATOR \
        -q \
        --illumina-read-length 150 \
        --illumina-prob-mismatch 0.0011 \
        --illumina-prob-mismatch-begin 0.00055 \
        --illumina-prob-mismatch-end 0.0022 \
        --fragment-mean-size 500 \
        --read-name-prefix "${ACCESSION}_" \
        -ir "$TEMP_FASTA" \
        -n "$N" \
        -o "$OUTPUT_DIR/reads/paired/depth_$DEPTH/mason_result/$ACCESSION_1.fasta" \
        -or "$OUTPUT_DIR/reads/paired/depth_$DEPTH/mason_result/$ACCESSION_2.fasta"

    rm -f "$TEMP_FASTA"
    rm -f "$TEMP_FASTA.fai"
done

echo "Make query file and shuffle it"
cat $OUTPUT_DIR/reads/paired/depth_$DEPTH/mason_result/*_1.fasta > "$OUTPUT_DIR/reads/paired/depth_$DEPTH/query_1.fasta"
cat $OUTPUT_DIR/reads/paired/depth_$DEPTH/mason_result/*_2.fasta > "$OUTPUT_DIR/reads/paired/depth_$DEPTH/query_2.fasta"
python "/fast/lunajang/metabuli/exclusion_test/Metabuli_DB/shuffle_paired_fasta.py" \
    --query_l "$OUTPUT_DIR/reads/paired/depth_$DEPTH/query_1.fasta" \
    --query_r "$OUTPUT_DIR/reads/paired/depth_$DEPTH/query_2.fasta" \
    --shuffled_query_l "$OUTPUT_DIR/reads/paired/depth_$DEPTH/shuffled_query_1.fasta" \
    --shuffled_query_r "$OUTPUT_DIR/reads/paired/depth_$DEPTH/shuffled_query_2.fasta"

echo ""

rm -rf "$OUTPUT_DIR/temp_paired_$DEPTH"