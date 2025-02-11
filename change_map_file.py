import sys

# 입력 및 출력 파일 지정
input_file = "/fast/lunajang/metabuli/exclusion_test/new_metabuli/fasta/taxid_fasta_mapping.txt"   # 기존 파일 (예: taxid: GCF_xxxxxx.x)
output_file = "/fast/lunajang/metabuli/exclusion_test/new_metabuli/fasta/assacc_to_taxid.tsv"  # 변환된 파일 (GCA_xxxxxx.x \t taxid)

with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    for line in infile:
        taxid, accession = line.strip().split(": ")
        outfile.write(f"{accession}\t{taxid}\n")

print(f"Converted file saved to {output_file}")
