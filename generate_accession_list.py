import os
import pandas as pd
import subprocess
from pathlib import Path
import argparse


def generate_target_genus(metadata_file, num_genus=None, genus_filter=None, min_species_per_genus=3):
    """
    Select target genera and filter metadata for selected genera.
    """
    # Load metadata
    metadata = pd.read_csv(metadata_file, sep="\t", low_memory=False)
    metadata["genus"] = metadata["gtdb_taxonomy"].str.split(";").str[-2].str.strip()

    # Filter metadata by quality and assembly level
    filtered_metadata = metadata[
        (metadata["checkm_completeness"] > 90)
        & (metadata["checkm_contamination"] < 5)
        & (metadata["ncbi_assembly_level"].isin(["Complete Genome", "Chromosome"]))
    ]

    # Filter by specific genus if provided
    if genus_filter:
        filtered_metadata = filtered_metadata[filtered_metadata["genus"].isin(genus_filter)]

    # Filter by number of species per genus
    genus_counts = filtered_metadata.groupby("genus").size()
    valid_genera = genus_counts[genus_counts >= min_species_per_genus].index
    filtered_metadata = filtered_metadata[filtered_metadata["genus"].isin(valid_genera)]
    print(f"Filtered to {len(valid_genera)} genera with at least {min_species_per_genus} species/subspecies.")    
    
    # Randomly sample genera if num_genus is specified
    unique_genera = filtered_metadata["genus"].dropna().unique()
    if num_genus is not None:
        if num_genus > len(unique_genera):
            print(f"Requested number of genera ({num_genus}) exceeds available genera ({len(unique_genera)}). Using all available genera.")
            num_genus = len(unique_genera)

        selected_genera = pd.Series(unique_genera).sample(n=num_genus, random_state=42)
        filtered_metadata = filtered_metadata[filtered_metadata["genus"].isin(selected_genera)]  

    print(f"Selected {filtered_metadata['genus'].nunique()} unique genera.")
    return filtered_metadata[["genus", "ncbi_species_taxid"]]
    # return filtered_metadata[["genus", "ncbi_taxid"]]


def filter_redundant_taxids(filtered_assemblies):
    """
    Randomly select one representative assembly per taxid.
    """
    # taxid별로 랜덤하게 하나의 assembly 선택
    unique_assemblies = filtered_assemblies.groupby("taxid").sample(n=1, random_state=42)

    print(f"Filtered to {len(unique_assemblies)} unique taxids (random selection).")
    return unique_assemblies


def make_accession_file(filtered_metadata, assembly_summary_file, output_dir):
    """
    Generate accession list from filtered metadata and assembly summary file.
    """
    assembly_summary = pd.read_csv(assembly_summary_file, sep="\t", skiprows=1)
    assembly_summary.columns = assembly_summary.columns.str.lstrip("#")
    os.makedirs(output_dir, exist_ok=True)
    
    # strain_taxids = filtered_metadata['ncbi_taxid'].dropna().unique()
    strain_taxids = filtered_metadata['ncbi_species_taxid'].dropna().unique()
    filtered_assemblies = assembly_summary[assembly_summary['taxid'].isin(strain_taxids)]


    # Create a temporary file to store accession IDs
    accession_all_file = os.path.join(output_dir, "accessions_all.txt")
    accession_all_list = filtered_assemblies["assembly_accession"].unique()
    with open(accession_all_file, "w") as f:
        for accession in accession_all_list:
            f.write(accession + "\n")
            
    # Create a temporary file to store accession IDs
    accession_file = os.path.join(output_dir, "accessions.txt")
    filtered_assemblies_unique = filter_redundant_taxids(filtered_assemblies)
    accession_list = filtered_assemblies_unique["assembly_accession"].unique()
    with open(accession_file, "w") as f:
        for accession in accession_list:
            f.write(accession + "\n")
    
    print(f"Accession file saved to: {accession_file}")    
    return filtered_assemblies[["assembly_accession", "taxid"]]
            
            
def save_taxid_fasta_mapping(filtered_assemblies, output_dir):
    """
    Save a mapping of taxid to expected FASTA file names.
    Args:
        filtered_assemblies (pd.DataFrame): Filtered assemblies with taxid and assembly_accession.
        output_dir (str): Directory to save the mapping file.
    """
    os.makedirs(output_dir, exist_ok=True)

    taxid_mapping_file = os.path.join(output_dir, "taxid_fasta_mapping.txt")

    grouped_assemblies = filtered_assemblies.groupby("taxid")
    with open(taxid_mapping_file, "w") as f:
        for taxid, group in grouped_assemblies:
            # Convert taxid to string and generate accession list
            accession_ids = ", ".join(map(str, group["assembly_accession"].tolist()))
            f.write(f"{taxid}: {accession_ids}\n")

    print(f"Taxid to FASTA mapping file saved to: {taxid_mapping_file}")


def save_genus_fasta_mapping(filtered_assemblies, filtered_metadata, output_dir):
    """
    Save a mapping of genus to expected FASTA file names.
    """
    os.makedirs(output_dir, exist_ok=True)

    # Merge filtered_metadata with filtered_assemblies to include genus information
    merged_data = pd.merge(
        filtered_assemblies,
        filtered_metadata,
        left_on="taxid",
        right_on="ncbi_species_taxid",
        # right_on="ncbi_taxid",
        how="inner"
    )
    
    merged_data = merged_data.drop_duplicates(subset=["genus", "assembly_accession"])

    genus_mapping_file = os.path.join(output_dir, "genus_fasta_mapping.txt")

    # Group by genus and write the mapping to file
    grouped_assemblies = merged_data.groupby("genus")
    with open(genus_mapping_file, "w") as f:
        for genus, group in grouped_assemblies:
            accession_ids = ", ".join(group["assembly_accession"].tolist())
            f.write(f"{genus}: {accession_ids}\n")

    print(f"Genus to FASTA mapping file saved to: {genus_mapping_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare Accession and Genus Mapping Files")
    parser.add_argument("--metadata", type=str, required=True, help="GTDB metadata file path")
    parser.add_argument("--assembly_summary", type=str, required=True, help="refseq assembly summary file path")
    parser.add_argument("--output", type=str, required=True, help="Output directory for accession and mapping files")
    parser.add_argument("--num_genus", type=int, default=None, help="Number of genera to select")
    parser.add_argument("--genus_filter", type=str, nargs="*", default=None, help="Specific genera to include")
    parser.add_argument("--min_species_per_genus", type=int, default=3, help="Minimum number of species per genus")

    args = parser.parse_args()

    if args.genus_filter is None and args.num_genus is None:
        raise ValueError("You must specify either --genus_filter or --num_genus.")

    # Generate target genera and filter metadata
    print("Filtering target genus from metadata...")
    filtered_metadata = generate_target_genus(args.metadata, args.num_genus, args.genus_filter, args.min_species_per_genus)

    # Download genomes
    print("Generating accession file...")
    filtered_assemblies = make_accession_file(filtered_metadata, args.assembly_summary, args.output)

    # Save genus-specific mapping files
    print("Generating genus-specific mapping files...")
    save_taxid_fasta_mapping(filtered_assemblies, args.output)
    save_genus_fasta_mapping(filtered_assemblies, filtered_metadata, args.output)

    print("Preparation completed. You can now run ncbi-genome-download manually.")