from collections import defaultdict
import argparse

def check_multi_genus_seq(mapping_file):
    accession_to_genus = defaultdict(set)

    with open(mapping_file, 'r') as f:
        for line in f:
            genus, accessions = line.strip().split(": ")
            for acc in accessions.split(", "):
                accession_to_genus[acc].add(genus)

    # 여러 genus에 속한 accession 찾기
    multi_genus_accessions = {acc: genus for acc, genus in accession_to_genus.items() if len(genus) > 1}

    if multi_genus_accessions:
        print("❌ 하나의 sequence가 여러 genus에 속해 있습니다.")
        print(f"중복된 accession 개수: {len(multi_genus_accessions)}")
        print("예시:", list(multi_genus_accessions.items())[:10])
    else:
        print("✅ 모든 sequence는 하나의 genus에만 속해 있습니다.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare Accession and Genus Mapping Files")
    parser.add_argument("--mapping_file", type=str, required=True, help="genus_fasta_mapping.txt file path")
    args = parser.parse_args()    

    check_multi_genus_seq(args.mapping_file)