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


def determine_query_and_reference(accessions, num_query=None):
    """
    Split accession IDs into query and reference sets.
    """
    num_species = len(accessions)
    if num_query is None:
        if num_species > 4:
            num_query = int(0.2 * num_species)
        else:
            num_query = 1
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


def write_fasta_file(output_file, content_list):
    """
    Write the combined content to a FASTA file.
    """
    with open(output_file, 'a+') as f:
        for content in content_list:
            f.write(content)


def split_and_merge_fasta(mapping_file, fasta_dir, output_dir):
    """
    Split sequences into reference and query for each genus and save to FASTA files.
    """
    os.makedirs(output_dir, exist_ok=True)
    query_file = os.path.join(output_dir, "query.fasta")
    reference_file = os.path.join(output_dir, "reference.fasta")
    open(query_file, 'w').close()
    open(reference_file, 'w').close()
    
    genus_to_accessions = parse_mapping_file(mapping_file)

    for genus, accessions in genus_to_accessions.items():
        print(f"Processing genus: {genus}")
        genus_file = os.path.join(output_dir, f"{genus}.fasta")
        open(genus_file, 'w').close()
        
        # Split accessions into query and reference
        query_accessions, reference_accessions = determine_query_and_reference(accessions)
        
        # Collect FASTA contents
        query_content = []
        reference_content = []
        for accession in accessions:
            content = process_accession_file(accession, fasta_dir)
            if content:
                if accession in query_accessions:
                    query_content.append(content)
                else:
                    reference_content.append(content)

        # Write to files
        write_fasta_file(query_file, query_content)
        write_fasta_file(reference_file, reference_content)
        write_fasta_file(genus_file, query_content + reference_content)
        
        print(f"  Query from genus {genus} added to: {query_file}")
        print(f"  Reference from genus {genus} added to: {reference_file}")
        print(f"  Genus sequences saved to: {genus_file}")

    print(f"FASTA splitting completed. Query saved to {query_file}, Reference saved to {reference_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare Accession and Genus Mapping Files")
    parser.add_argument("--mapping_file", type=str, required=True, help="GTDB metadata file path")
    parser.add_argument("--fasta_dir", type=str, required=True, help="refseq assembly summary file path")
    parser.add_argument("--output", type=str, required=True, help="Output directory for accession and mapping files")
    
    args = parser.parse_args()
    

    split_and_merge_fasta(args.mapping_file, args.fasta_dir, args.output)
