import random

def read_fasta(filepath):
    with open(filepath, 'r') as f:
        lines = f.read().splitlines()
        return [(lines[i], lines[i + 1]) for i in range(0, len(lines), 2)]

def write_fasta(pairs, filepath):
    with open(filepath, 'w') as f:
        for header, seq in pairs:
            f.write(f"{header}\n{seq}\n")

# Read paired-end FASTA files
queries = read_fasta("/fast/lunajang/metabuli/exclusion_test/new_metabuli/fasta/reads/query/query.fasta")

# Shuffle the combined pairs
random.shuffle(queries)

# Write the shuffled pairs back to new FASTA files
write_fasta(queries, "/fast/lunajang/metabuli/exclusion_test/new_metabuli/fasta/reads/query/shuffled_query.fasta")

print("Shuffling complete!")
