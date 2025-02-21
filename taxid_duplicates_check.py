from collections import defaultdict
import argparse

def check_ref_query_speceis_level(mapping_file, query_list, reference_list):
    accession_to_species = {}

    with open(mapping_file) as f:
        for line in f:
            acc, taxid = line.strip().split("\t")
            accession_to_species[acc] = taxid 

    query_accessions = set(open(query_list).read().splitlines())
    reference_accessions = set(open(reference_list).read().splitlines())

    query_species = {accession_to_species[acc] for acc in query_accessions if acc in accession_to_species}
    reference_species = {accession_to_species[acc] for acc in reference_accessions if acc in accession_to_species}

    overlapping_species = query_species & reference_species

    if overlapping_species:
        print("❌ Query와 Reference가 같은 taxid를 가지는 경우가 존재합니다.")
        print(f"중복된 species 개수: {len(overlapping_species)}")
        print(f"중복된 species 예시: {list(overlapping_species)[:10]}")
    else:
        print("✅ Query와 Reference는 모두 다른 taxid를 가지고 있습니다.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare Accession and Genus Mapping Files")
    parser.add_argument("--mapping_file", type=str, required=True, help="GTDB taxdump taxid.map file path")
    parser.add_argument("--query_list", type=str, required=True, help="query.list file path")
    parser.add_argument("--reference_list", type=str, required=True, help="reference.list file path")
    args = parser.parse_args()    

    check_ref_query_speceis_level(args.mapping_file, args.query_list, args.reference_list)