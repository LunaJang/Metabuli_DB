import os
import gzip
import glob
import random
import argparse


def parse_mapping_file(mapping_file):
    """
    Parse the genus_fasta_mapping.txt file.
    """
    genus_to_accessions = {}
    with open(mapping_file, 'r') as f:
        for line in f:
            genus, accessions = line.strip().split(": ")
            genus_to_accessions[genus] = accessions.split(", ")
    return genus_to_accessions


def determine_query_and_reference(accessions, query_fraction, num_query=None):
    """
    Split accession IDs into query and reference sets.
    """
    num_species = len(accessions)
    num_query = max(1, int(query_fraction * num_species))
    query_accessions = set(random.sample(accessions, num_query))
    reference_accessions = set(accessions) - query_accessions
    return query_accessions, reference_accessions


def process_accession_file(accession, fasta_dir):
    """
    Read and return the content of a compressed FASTA file.
    """
    fasta_files = glob.glob(os.path.join(fasta_dir, accession, "*.fna.gz"))
    if fasta_files:
        with gzip.open(fasta_files[0], 'rt') as fasta_f:
            return fasta_f.read()
    else:
        print(f"Warning: No .fna.gz file found for {accession}.")
        return None


def split_and_merge_fasta(mapping_file, fasta_dir, output_dir, query_fraction):
    """
    Split sequences into reference and query for each genus and save to FASTA files.
    """
    genus_to_accessions = parse_mapping_file(mapping_file)

    query_file = os.path.join(output_dir, "query.list")
    reference_file = os.path.join(output_dir, "reference.list")
    open(query_file, 'w').close()
    open(reference_file, 'w').close()
    
    num_query = 0
    num_reference = 0
    
    for genus, accessions in genus_to_accessions.items():
        print(f"Processing genus: {genus}")
        
        # Split accessions into query and reference
        query_accessions, reference_accessions = determine_query_and_reference(accessions, query_fraction)
        
        # Collect FASTA contents
        num_genus_query = 0
        num_genus_reference = 0
        for accession in accessions:
            content = process_accession_file(accession, fasta_dir)
            if content:
                accession_file = os.path.join(fasta_dir, accession, f"{accession}.fna")
                open(accession_file, 'w').close()        
                with open(accession_file, 'a+') as f:
                    f.write(content)
                    
                if accession in query_accessions:
                    num_genus_query += 1
                    with open(query_file, "a") as f:
                        f.write(f"{accession_file}\n")
                else:
                    num_genus_reference += 1
                    with open(reference_file, "a") as f:
                        f.write(f"{accession_file}\n")
                        
        num_query += num_genus_query
        num_reference += num_genus_reference

    print(f"Query / reference splitting completed. Query saved to {query_file} ({num_query}), Reference saved to {reference_file} ({num_reference}).")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare Accession and Genus Mapping Files")
    parser.add_argument("--mapping_file", type=str, required=True, help="genus_fasta_mapping.txt file path")
    parser.add_argument("--fasta_dir", type=str, required=True, help="input directory for FASTA files")    
    parser.add_argument("--output", type=str, required=True, help="output directory for query and reference lists")
    parser.add_argument("--query_fraction", type=float, required=True, help="ratio of query")
    
    args = parser.parse_args()
    

    split_and_merge_fasta(args.mapping_file, args.fasta_dir, args.output, args.query_fraction)
