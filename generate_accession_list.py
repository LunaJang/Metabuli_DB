import os
import pandas as pd
import argparse


def load_files(metadata_file, taxid_map_file, assembly_summary_file):
    """Load metadata, taxid map, and assembly summary files with only required columns."""
    
    metadata = pd.read_csv(metadata_file, sep="\t", low_memory=False, usecols=["gtdb_taxonomy", "ncbi_species_taxid"])
    metadata["genus"] = metadata["gtdb_taxonomy"].str.split(";").str[-2].str.strip()
    metadata = metadata[["genus", "ncbi_species_taxid"]]
    
    taxid_map = pd.read_csv(taxid_map_file, sep="\t", header=None, names=["accession", "taxid"], usecols=[0, 1])
    
    assembly_summary = pd.read_csv(assembly_summary_file, sep="\t", skiprows=1)
    assembly_summary.columns = assembly_summary.columns.str.lstrip("#")
    assembly_summary = assembly_summary[["assembly_accession", "taxid"]]

    return metadata, taxid_map, assembly_summary


def merge_datasets(metadata, taxid_map, assembly_summary):
    """Merge GTDB metadata, RefSeq assembly summary, and GTDB taxdump to retain only valid accessions."""
    
    # Step 1: Merge GTDB metadata and RefSeq assembly summary using species taxid
    merged_data = pd.merge(metadata, assembly_summary, left_on="ncbi_species_taxid", right_on="taxid", how="inner")
    
    # Step 2: Merge the above with GTDB taxdump using accession
    merged_data = pd.merge(merged_data, taxid_map, left_on="assembly_accession", right_on="accession", how="inner")
    
    # taxid_x, taxid_y가 존재하는 경우
    if "taxid_y" in merged_data.columns:
        merged_data["taxid"] = merged_data["taxid_y"]  # 최종적으로 사용할 taxid 선택
        merged_data = merged_data.drop(columns=["taxid_x", "taxid_y"])  # 필요 없는 열 삭제


    print(f"After merging, {len(merged_data)} records remain.")
    return merged_data



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

    print(f"After filtering, {filtered_data['genus'].nunique()} genera remain.")
    print(f"After filtering, {len(filtered_data.drop_duplicates())} unique data points remain.")

    return filtered_data


def  save_accession_files(filtered_data, output_dir):
    """Save final accession list."""
    os.makedirs(output_dir, exist_ok=True)

    accession_file = os.path.join(output_dir, "accessions.txt")
    filtered_data["assembly_accession"].drop_duplicates().to_csv(accession_file, index=False, header=False)

    print(f"Accession file saved to: {accession_file}")


def save_taxid_fasta_mapping(filtered_data, output_dir):
    """Save a mapping of taxid to expected FASTA file names."""
    os.makedirs(output_dir, exist_ok=True)

    taxid_mapping_file = os.path.join(output_dir, "assacc_to_taxid.tsv")
    grouped_assemblies = filtered_data.drop_duplicates(subset=["taxid", "assembly_accession"])
    grouped_assemblies = grouped_assemblies.groupby("taxid")

    with open(taxid_mapping_file, "w") as f:
        for _, row in grouped_assemblies.iterrows():
            f.write(f"{row['assembly_accession']}\t{row['taxid']}\n")

    print(f"Taxid to FASTA mapping file saved to: {taxid_mapping_file}")


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

    print(f"Genus to FASTA mapping file saved to: {genus_mapping_file}")


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

    print("Loading files...")
    metadata, taxid_map, assembly_summary = load_files(args.metadata, args.gtdb_taxid, args.assembly_summary)

    print("Merging datasets...")
    merged_data = merge_datasets(metadata, taxid_map, assembly_summary)

    print("Filtering by genus...")
    filtered_data = filter_by_genus(merged_data, args.num_genus, args.min_species_per_genus, args.max_species_per_genus)

    print("Saving final accession list...")
    save_accession_files(filtered_data, args.output)

    print("Saving taxid mapping file...")
    save_taxid_fasta_mapping(filtered_data, args.output)

    print("Saving genus mapping file...")
    save_genus_fasta_mapping(filtered_data, args.output)

    print("Preparation completed. You can now run ncbi-genome-download manually.")
