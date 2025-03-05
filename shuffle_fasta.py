import random
import argparse

def read_fasta(filepath):
    with open(filepath, 'r') as f:
        lines = f.read().splitlines()
        return [(lines[i], lines[i + 1]) for i in range(0, len(lines), 2)]

def write_fasta(pairs, filepath):
    with open(filepath, 'w') as f:
        for header, seq in pairs:
            f.write(f"{header}\n{seq}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare Accession and Genus Mapping Files")
    parser.add_argument("--query", type=str, required=True, help="query.fasta file path")
    parser.add_argument("--shuffled_query", type=str, required=True, help="shuffled_query.fasta file path")
    args = parser.parse_args()   
    
    queries = read_fasta(args.query)
    # Shuffle the combined pairs
    random.shuffle(queries)
    # Write the shuffled pairs back to new FASTA files
    write_fasta(queries, args.shuffled_query)

    print("Shuffling complete!")
