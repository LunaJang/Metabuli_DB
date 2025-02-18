import os
import pandas as pd
import argparse
from collections import defaultdict

def load_files(metadata_file, taxid_map_file, assembly_summary_file):
    """Load metadata, taxid map, and assembly summary files with only required columns."""
    
    metadata = pd.read_csv(metadata_file, sep="\t", low_memory=False, usecols=["gtdb_taxonomy", "ncbi_species_taxid"])
    metadata["genus"] = metadata["gtdb_taxonomy"].str.split(";").str[-2].str.strip()
    metadata["species"] = metadata["gtdb_taxonomy"].str.split(";").str[-1].str.strip()
    metadata = metadata[["genus", "species", "ncbi_species_taxid"]].astype({"ncbi_species_taxid": "int32"})
    
    taxid_map = pd.read_csv(taxid_map_file, sep="\t", header=None, names=["accession", "gtdb_taxid"], usecols=[0, 1])
    
    assembly_summary = pd.read_csv(assembly_summary_file, sep="\t", skiprows=1)
    assembly_summary.columns = assembly_summary.columns.str.lstrip("#")
    assembly_summary = assembly_summary[["assembly_accession", "taxid"]].astype({"taxid": "int32"})

    return metadata, taxid_map, assembly_summary


def merge_datasets(metadata, taxid_map, assembly_summary):
    merged_data = metadata.merge(assembly_summary, left_on="ncbi_species_taxid", right_on="taxid", how="inner")
    merged_data = merged_data.merge(taxid_map, left_on="assembly_accession", right_on="accession", how="inner")

    print(f"\t\tAfter merging, {len(merged_data)} records remain.")
    return merged_data


def filter_duplicated_seq(merged_data):
    genus_counts = merged_data.groupby("assembly_accession")["genus"].nunique()
    duplicated_accessions = genus_counts[genus_counts > 1].index
    filtered_data = merged_data[~merged_data["assembly_accession"].isin(duplicated_accessions)].copy()
    
    print(f"\t\tFiltered out {len(duplicated_accessions)} sequences that appeared in multiple genera.")
    return filtered_data

def filter_by_species(merged_data):
    filtered_by_species = (merged_data.groupby("species", group_keys=False, as_index=False)
                            .apply(lambda x: x.sample(n=1, random_state=42), include_groups=False))
    filtered_data = (filtered_by_species.groupby("gtdb_taxid", group_keys=False, as_index=False)
                     .apply(lambda x: x.sample(n=1, random_state=42), include_groups=False))

    print(f"\t\tFiltered down to {len(filtered_data)} unique species and taxid representatives.")
    return filtered_data

def filter_by_genus(merged_data, num_genus=None, min_species_per_genus=3,
                    max_species_per_genus=20):
    """Filter merged data by genus count."""
    # Count number of species per genus
    genus_counts = merged_data.groupby("genus").size()
    valid_genera = genus_counts[genus_counts >= min_species_per_genus].index
    filtered_data = merged_data[merged_data["genus"].isin(valid_genera)]

    # Randomly sample genera if num_genus is specified
    if num_genus is not None and num_genus < len(valid_genera):
        selected_genera = pd.Series(valid_genera).sample(n=num_genus, random_state=42)
        filtered_data = filtered_data[filtered_data["genus"].isin(selected_genera)]

    # 각 속에서 최대 max_species_per_genus 개수만큼 무작위 샘플링
    filtered_data = (
        filtered_data.groupby("genus", group_keys=False)
        .apply(lambda group: group.sample(n=min(len(group), max_species_per_genus), random_state=42))
        .reset_index(drop=True)
    )

    # 'genus' 컬럼이 존재하는지 확인
    if "genus" not in filtered_data.columns:
        raise ValueError("Error: 'genus' column is missing after applying groupby and sample.")

    print(f"\t\tAfter filtering, {filtered_data['genus'].nunique()} genera remain.")
    print(f"\t\tAfter filtering, {len(filtered_data.drop_duplicates())} unique data points remain.")

    return filtered_data


def save_accession_files(filtered_data, output_dir):
    """Save final accession list."""
    os.makedirs(output_dir, exist_ok=True)

    accession_file = os.path.join(output_dir, "accessions.txt")
    filtered_data["assembly_accession"].drop_duplicates().to_csv(accession_file, index=False, header=False)

    print(f"\t\tAccession file saved to: {accession_file}")


def save_taxid_fasta_mapping(filtered_data, output_dir):
    """Save a mapping of taxid to expected FASTA file names."""
    os.makedirs(output_dir, exist_ok=True)

    taxid_mapping_file = os.path.join(output_dir, "assacc_to_taxid.tsv")
    
    # 중복 제거 후 데이터프레임 유지
    grouped_assemblies = filtered_data.drop_duplicates(subset=["gtdb_taxid", "assembly_accession"])
    
    # 파일 저장
    with open(taxid_mapping_file, "w") as f:
        for _, row in grouped_assemblies.iterrows():  # iterrows 사용 가능
            f.write(f"{row['assembly_accession']}\t{row['gtdb_taxid']}\n")

    print(f"\t\tTaxid to FASTA mapping file saved to: {taxid_mapping_file}")
    

def save_genus_fasta_mapping(filtered_data, output_dir):
    """Save a mapping of genus to expected FASTA file names."""
    os.makedirs(output_dir, exist_ok=True)

    genus_mapping_file = os.path.join(output_dir, "genus_fasta_mapping.txt")
    grouped_assemblies = filtered_data.drop_duplicates(subset=["genus", "assembly_accession"])
    grouped_assemblies = grouped_assemblies.groupby("genus")

    with open(genus_mapping_file, "w") as f:
        for genus, group in grouped_assemblies:
            accession_ids = ", ".join(group["assembly_accession"].tolist())
            f.write(f"{genus}: {accession_ids}\n")

    print(f"\t\tGenus to FASTA mapping file saved to: {genus_mapping_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare Accession and Genus Mapping Files")
    parser.add_argument("--metadata", type=str, required=True, help="GTDB metadata file path")
    parser.add_argument("--assembly_summary", type=str, required=True, help="RefSeq assembly summary file path")
    parser.add_argument("--gtdb_taxid", type=str, required=True, help="GTDB taxdump taxid map file path")
    parser.add_argument("--output", type=str, required=True, help="Output directory for accession files")
    parser.add_argument("--num_genus", type=int, default=None, help="Number of genera to select")
    parser.add_argument("--min_species_per_genus", type=int, default=5, help="Minimum number of species per genus")
    parser.add_argument("--max_species_per_genus", type=int, default=20, help="Maximum number of species per genus")
    args = parser.parse_args()

    print("\tLoading files...")
    metadata, taxid_map, assembly_summary = load_files(args.metadata, args.gtdb_taxid, args.assembly_summary)

    print("\tMerging datasets...")
    merged_data = merge_datasets(metadata, taxid_map, assembly_summary)

    print("\tFiltering...")
    filtered_data = filter_duplicated_seq(merged_data)
    filtered_data = filter_by_species(filtered_data)
    filtered_data = filter_by_genus(merged_data, args.num_genus, args.min_species_per_genus, args.max_species_per_genus)

    print("\tSaving final accession list...")
    save_accession_files(filtered_data, args.output)

    print("\tSaving taxid mapping file...")
    save_taxid_fasta_mapping(filtered_data, args.output)

    print("\tSaving genus mapping file...")
    save_genus_fasta_mapping(filtered_data, args.output)

    print("\tPreparation completed. You can now run ncbi-genome-download manually.")
