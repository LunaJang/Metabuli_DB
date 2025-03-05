import random

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
    parser.add_argument("--query_l", type=str, required=True, help="query.fasta file path")
    parser.add_argument("--query_r", type=str, required=True, help="query.fasta file path")
    parser.add_argument("--shuffled_query_r", type=str, required=True, help="shuffled_query.fasta file path")
    parser.add_argument("--shuffled_query_l", type=str, required=True, help="shuffled_query.fasta file path")
    args = parser.parse_args()   
    
    pairs_1 = read_fasta(args.query_l)
    pairs_2 = read_fasta(args.query_r)
    
    assert len(pairs_1) == len(pairs_2), "The number of reads in query_1 and query_2 must match."
    
    combined_pairs = list(zip(pairs_1, pairs_2))
    random.shuffle(combined_pairs)
    shuffled_pairs_1, shuffled_pairs_2 = zip(*combined_pairs)    
    
    write_fasta(shuffled_pairs_1, args.shuffled_query_r)
    write_fasta(shuffled_pairs_2, args.shuffled_query_l)
    print("Shuffling complete!")

