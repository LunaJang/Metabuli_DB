#!/bin/bash
#SBATCH --job-name=metabuli_db
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 2
#SBATCH --time=10-23:00:00
#SBATCH --partition=gpu
#SBATCH --nodelist=devbox002
#SBATCH --output=/fast/lunajang/metabuli/exclusion_test/new_metabuli/logs/%j_output.log
#SBATCH --error=/fast/lunajang/metabuli/exclusion_test/new_metabuli/logs/%j_error.log

set -e

source ~/.bashrc
source activate base
conda activate metabuli

OUTPUT_DIR="/fast/lunajang/metabuli/exclusion_test/new_metabuli"
NUM_GENUS=700
MIN_SPECIES_PER_GENUS=3
QUERY_FRACTION=0.3                   
ACCESSION2TAXID_SCRIPT="/home/lunajang/src/Metabuli/util/prepare_gtdb_taxonomy.sh"
GTDB_TAXDUMP="/fast/lunajang/metabuli/exclusion_test/new_metabuli/gtdb-taxdump/R220"
METABULI="/home/lunajang/src/Metabuli/build_devbox/bin/metabuli"  
MASON2_SIMULATOR="/home/lunajang/src/mason2-2.0.9-Linux-x86_64_sse4/bin/mason_simulator"
DEPTH=0.25        

# mkdir -p "$OUTPUT_DIR/fasta"

echo "1. Download FASTA files"
echo "1-1. Download GTDB Metadata"
# wget https://data.gtdb.ecogenomic.org/releases/release220/220.0/bac120_metadata_r220.tsv.gz -O "$OUTPUT_DIR/bac120_metadata_r220.tsv.gz"
# gunzip -c "$OUTPUT_DIR/bac120_metadata_r220.tsv.gz" > "$OUTPUT_DIR/bac120_metadata_r220.tsv"
# rm "$OUTPUT_DIR/bac120_metadata_r220.tsv.gz"

echo "1-2. Download assembly summary"
# wget "ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt" -O "$OUTPUT_DIR/assembly_summary.txt"

echo "1-3. Generate accession list"
python "/fast/lunajang/metabuli/exclusion_test/Metabuli_DB/generate_accession_list.py" \
    --metadata "$OUTPUT_DIR/bac120_metadata_r220.tsv" \
    --assembly_summary "$OUTPUT_DIR/assembly_summary.txt" \
    --gtdb_taxid "$GTDB_TAXDUMP/taxid.map" \
    --output "$OUTPUT_DIR/fasta" \
    --num_genus "$NUM_GENUS"\
    --min_species_per_genus "$MIN_SPECIES_PER_GENUS"

python "/fast/lunajang/metabuli/exclusion_test/Metabuli_DB/multi_genus_sequence_check.py" \
    --mapping_file "$OUTPUT_DIR/fasta/genus_fasta_mapping.txt" 

echo "1-4. Download FASTA files"
ncbi-genome-download --formats fasta -o "$OUTPUT_DIR/fasta" -A "$OUTPUT_DIR/fasta/accessions.txt" bacteria --section refseq --retries 10 --parallel 1

echo ""

echo "2. Separate query and reference"
python "/fast/lunajang/metabuli/exclusion_test/Metabuli_DB/get_reference_query_fa.py" \
    --mapping_file "$OUTPUT_DIR/fasta/genus_fasta_mapping.txt" \
    --fasta_dir "$OUTPUT_DIR/fasta/refseq/bacteria" \
    --output "$OUTPUT_DIR/fasta" \
    --query_fraction "$QUERY_FRACTION"

python "/fast/lunajang/metabuli/exclusion_test/Metabuli_DB/taxid_duplicates_check.py" \
    --mapping_file "$GTDB_TAXDUMP/taxid.map" \
    --query_list "$OUTPUT_DIR/fasta/query.list" \
    --reference_list "$OUTPUT_DIR/fasta/reference.list"

echo ""

echo "3. Prepare GTDB taxonomy and accession2taxid"
echo "3-1. Download NCBI-style taxonomy dump"
# latest_release=$(curl -s https://api.github.com/repos/shenwei356/gtdb-taxdump/releases/latest | jq -r '.assets[] | select(.name | contains("gtdb-taxdump.tar.gz")) | .browser_download_url')
# wget gtdb-taxdump.tar.gz "$latest_release" -O "$OUTPUT_DIR/gtdb-taxdump.tar.gz"
# tar -xvzf "$OUTPUT_DIR/gtdb-taxdump.tar.gz" -C "$OUTPUT_DIR"
# rm "$OUTPUT_DIR/gtdb-taxdump.tar.gz"

echo "3-2. Prepare accession2taxid"
$METABULI editNames "$GTDB_TAXDUMP/names.dmp" "$GTDB_TAXDUMP/taxid.map"
$METABULI accession2taxid  "$OUTPUT_DIR/fasta/reference.list" "$GTDB_TAXDUMP/taxid.map"

echo ""

echo "4. Make query reads"
echo "4-1. Run Mason2 simulator"
mkdir "$OUTPUT_DIR/fasta/reads"
mkdir "$OUTPUT_DIR/fasta/reads/query"
mkdir "$OUTPUT_DIR/fasta/temp"
awk -F '/' '{print $0, $NF }' "$OUTPUT_DIR/fasta/query.list" | while read -r FNA_FILE ACCESSION; do
    ACCESSION=$(basename "$ACCESSION" .fna)
    echo "Mason2 for $ACCESSION, FNA_FILE: $FNA_FILE"

    TEMP_FASTA="$OUTPUT_DIR/fasta/temp/temp_$ACCESSION.fasta"
    seqkit seq -m 500 -g "$FNA_FILE" > "$TEMP_FASTA"

    if [ ! -s "$TEMP_FASTA" ]; then
        echo "Skipping $ACCESSION (no sequences ≥ 500bp)"
        rm -f "$TEMP_FASTA"
        continue
    fi
    
    TEMP_FASTA="$OUTPUT_DIR/fasta/temp/temp_$ACCESSION.fasta"
    seqkit seq -m 500 -g "$FNA_FILE" > "$TEMP_FASTA"

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
        -o "$OUTPUT_DIR/fasta/reads/depth_$DEPTH/mason_result/$ACCESSION.fasta"
        #--force-single-end \

    rm -f "$TEMP_FASTA"
    rm -f "$TEMP_FASTA.fai"

done

echo "4-2. Make query file and shuffle it"
cat $OUTPUT_DIR/fasta/reads/depth_$DEPTH/mason_result/*.fasta > "$OUTPUT_DIR/fasta/reads/depth_$DEPTH/query.fasta"
python "/fast/lunajang/metabuli/exclusion_test/Metabuli_DB/shuffle_fasta.py" \
    --query "$OUTPUT_DIR/fasta/reads/depth_$DEPTH/query.fasta" \
    --shuffled_query "$OUTPUT_DIR/fasta/reads/depth_$DEPTH/shuffled_query.fasta"

echo ""

echo "5. Create Metabuli Custom Database"
$METABULI build "$OUTPUT_DIR/db" "$OUTPUT_DIR/fasta/reference.list" "$GTDB_TAXDUMP/taxid.accession2taxid" --taxonomy-path "$GTDB_TAXDUMP" 

echo "Metabuli Database generated: $OUTPUT_DIR/db"
