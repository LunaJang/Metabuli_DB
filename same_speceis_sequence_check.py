from collections import defaultdict

taxid_mapping_file = "/fast/lunajang/metabuli/exclusion_test/new_metabuli/fasta/assacc_to_taxid.tsv"

# Accession-to-TaxID 매핑 불러오기
accession_to_species = {}

with open(taxid_mapping_file) as f:
    for line in f:
        acc, taxid = line.strip().split("\t")
        accession_to_species[acc] = taxid  # TaxID가 species level인지 확인 필요

# Query 및 Reference 목록 불러오기
query_list_file = "/fast/lunajang/metabuli/exclusion_test/new_metabuli/fasta/query.list"
reference_list_file = "/fast/lunajang/metabuli/exclusion_test/new_metabuli/fasta/reference.list"

query_accessions = set(open(query_list_file).read().splitlines())
reference_accessions = set(open(reference_list_file).read().splitlines())

# Query 및 Reference의 Species 목록 생성
query_species = {accession_to_species[acc] for acc in query_accessions if acc in accession_to_species}
reference_species = {accession_to_species[acc] for acc in reference_accessions if acc in accession_to_species}

# Query와 Reference에서 중복된 Species 찾기
overlapping_species = query_species & reference_species

if overlapping_species:
    print("❌ Query와 Reference가 같은 species를 포함하고 있습니다.")
    print(f"중복된 species 개수: {len(overlapping_species)}")
    print(f"중복된 species 예시: {list(overlapping_species)[:10]}")
else:
    print("✅ Query와 Reference는 species level에서 완전히 분리되었습니다.")
