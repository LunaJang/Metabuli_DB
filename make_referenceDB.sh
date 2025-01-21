#!/bin/bash
#SBATCH --job-name=metabuli_db
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 1
#SBATCH --time=10-23
#SBATCH --partition=compute
#SBATCH --nodelist=super003
#SBATCH --output=/fast/lunajang/metabuli/exclusion_test/new_metabuli/logs/%j_output.log
#SBATCH --error=/fast/lunajang/metabuli/exclusion_test/new_metabuli/logs/%j_error.log

set -e

source ~/.bashrc
source activate base
conda activate metabuli

OUTPUT_DIR="/fast/lunajang/metabuli/exclusion_test/new_metabuli"
NUM_GENUS=500
MIN_SPECIES_PER_GENUS=3
QUERY_FRACTION=0.2                   
ACCESSION2TAXID_SCRIPT="/home/lunajang/src/Metabuli/util/prepare_gtdb_taxonomy.sh"
GTDB_TAXDUMP="/fast/lunajang/metabuli/exclusion_test/new_metabuli/gtdb-taxdump/R202"
METABULI="/home/lunajang/src/Metabuli/build/bin/metabuli"  

# mkdir -p "$OUTPUT_DIR/taxonomy" "$OUTPUT_DIR/fasta"  "$OUTPUT_DIR/metabuli_db" 

echo "1. Download FASTA files"
echo "1-1. Download GTDB Metadata"
# wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release202/202.0/bac120_metadata_r202.tar.gz -O "$OUTPUT_DIR/bac120_metadata_r202.tar.gz"
# tar -xvzf "$OUTPUT_DIR/bac120_metadata_r202.tar.gz" -C "$OUTPUT_DIR"
# rm "$OUTPUT_DIR/bac120_metadata_r202.tar.gz"

echo "1-2. Download assembly summary"
# wget "ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt" -O "$OUTPUT_DIR/assembly_summary.txt"

echo "1-3. Generate accession list"
# python "/fast/lunajang/metabuli/exclusion_test/new_metabuli/code/generate_accession_list.py" \
#     --metadata "$OUTPUT_DIR/bac120_metadata_r202.tsv" \
#     --assembly_summary "$OUTPUT_DIR/assembly_summary.txt" \
#     --output "$OUTPUT_DIR/fasta" \
#     --num_genus "$NUM_GENUS"\
#     --min_species_per_genus "$MIN_SPECIES_PER_GENUS"

echo "1-4. Download FASTA files"
# ncbi-genome-download --formats fasta -o "$OUTPUT_DIR/fasta" -A "$OUTPUT_DIR/fasta/accessions.txt" bacteria --section refseq


echo "2. Separate query and reference"
python "/fast/lunajang/metabuli/exclusion_test/new_metabuli/code/get_reference_query_fa.py" \
    --mapping_file "$OUTPUT_DIR/fasta/genus_fasta_mapping.txt" \
    --fasta_dir "$OUTPUT_DIR/refseq/bacteria" \
    --output "$OUTPUT_DIR/fasta/contigs"


echo "3. Run Mason2"
# mason2 illumina -n 1000000 -N 1000000 -i -sq -o "$OUTPUT_DIR/fasta/reads/" "$OUTPUT_DIR/fasta/contigs/query.fa"


echo "4. Prepare GTDB taxonomy and accession2taxid"
echo "4-1. Download NCBI-style taxonomy dump"
# latest_release=$(curl -s https://api.github.com/repos/shenwei356/gtdb-taxdump/releases/latest | jq -r '.assets[] | select(.name | contains("gtdb-taxdump.tar.gz")) | .browser_download_url')
# wget gtdb-taxdump.tar.gz "$latest_release" -O "$OUTPUT_DIR/gtdb-taxdump.tar.gz"
# tar -xvzf "$OUTPUT_DIR/gtdb-taxdump.tar.gz" -C "$OUTPUT_DIR"
# rm "$OUTPUT_DIR/gtdb-taxdump.tar.gz"

echo "4-2. Prepare accession2taxid"
# $METABULI editNames "$GTDB_TAXDUMP/names.dmp" "$GTDB_TAXDUMP/taxid.map"
# echo "$OUTPUT_DIR/fasta/reference.fa" > "$OUTPUT_DIR/metabuli_db/reference_list.txt"
# $METABULI accession2taxid  <FASTA_LIST> <GTDB_TAXDUMP/taxid.map>

echo "5. Create Metabuli Custom Database"
# metabuli build <DBDIR> <FASTA_LIST> <GTDB_TAXDUMP/taxid.accession2taxid> --taxonomy-path GTDB_TAXDUMP [options]

# echo "Metabuli Database generated: $OUTPUT_DIR/metabuli_db"
