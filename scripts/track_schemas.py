import argparse
import pandas as pd
import numpy as np
import itertools
from collections import defaultdict
import os
import csv

def load_data(filepath):
    data = []
    with open(filepath, 'r') as f:
        header = f.readline().strip().split(',')
        for line in f:
            parts = line.strip().split(',')
            if not parts or len(parts) < 4:
                continue
            cycle = int(parts[0])
            fitness = float(parts[-1])
            sequence = parts[-2]
            data.append({
                'cycle': cycle,
                'sequence': sequence,
                'fitness': fitness
            })
    return pd.DataFrame(data)

def format_schema(schema_tuple, length):
    parts = [f"{pos}{aa}" for pos, aa in schema_tuple]
    return "_".join(parts)

def main():
    parser = argparse.ArgumentParser(description="Track schemas across GA generations")
    parser.add_argument('--csv', required=True, help="Path to all_sequences.csv")
    parser.add_argument('--k', type=int, default=3, help="Order of the schema")
    parser.add_argument('--out_dir', default='.', help="Directory to save the results")
    parser.add_argument('--contiguous', action='store_true', help="Only track contiguous k-mers")
    parser.add_argument('--min_count', type=int, default=2, help="Minimum count for a schema to be tracked")
    args = parser.parse_args()

    print(f"Loading data from {args.csv}...", flush=True)
    df = load_data(args.csv)

    cycles = sorted(df['cycle'].unique())

    print(f"Tracking Schemas of order k={args.k}", flush=True)
    if args.contiguous:
        print("Using CONTIGUOUS schemas only.", flush=True)
    else:
        print("Using NON-CONTIGUOUS schemas.", flush=True)

    os.makedirs(args.out_dir, exist_ok=True)
    out_path = os.path.join(args.out_dir, 'schema_statistics.csv')
    
    # Initialize CSV writer
    with open(out_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['cycle', 'schema', 'count_overall', 'freq_overall', 'count_top10', 'freq_top10', 'mean_fitness', 'enrichment_score'])

        for cycle in cycles:
            cycle_df = df[df['cycle'] == cycle]
            pop_size = len(cycle_df)
            
            top10_cutoff = np.percentile(cycle_df['fitness'], 90)
            
            schema_data = defaultdict(lambda: {'count': 0, 'count_top10': 0, 'sum_fitness': 0.0})

            seqs = cycle_df['sequence'].tolist()
            fits = cycle_df['fitness'].tolist()
            seq_len = len(seqs[0]) if seqs else 0
            
            for seq, fit in zip(seqs, fits):
                is_top10 = (fit >= top10_cutoff)

                if args.contiguous:
                    if seq_len >= args.k:
                        for i in range(seq_len - args.k + 1):
                            schema = tuple((i + j, seq[i + j]) for j in range(args.k))
                            schema_data[schema]['count'] += 1
                            schema_data[schema]['sum_fitness'] += fit
                            if is_top10:
                                schema_data[schema]['count_top10'] += 1
                else:
                    for schema in itertools.combinations(enumerate(seq), args.k):
                        schema_data[schema]['count'] += 1
                        schema_data[schema]['sum_fitness'] += fit
                        if is_top10:
                            schema_data[schema]['count_top10'] += 1

            top10_size = sum(f >= top10_cutoff for f in fits)
            if top10_size == 0: top10_size = 1
            
            # Write out passing schemas directly to save RAM
            saved_count = 0
            for schema, data in schema_data.items():
                if data['count'] < args.min_count:
                    continue
                    
                saved_count += 1
                freq_overall = data['count'] / pop_size
                freq_top10 = data['count_top10'] / top10_size
                mean_fitness = data['sum_fitness'] / data['count']
                
                eps = 1e-4
                enrichment = np.log2((freq_top10 + eps) / (freq_overall + eps))

                schema_str = format_schema(schema, seq_len)
                writer.writerow([
                    cycle, schema_str, data['count'], freq_overall, 
                    data['count_top10'], freq_top10, mean_fitness, enrichment
                ])
                
            print(f"Processed cycle {cycle}: {len(schema_data)} unique schemas, {saved_count} saved (>= {args.min_count} count).", flush=True)

    print(f"Saved statistics to {out_path}", flush=True)

if __name__ == "__main__":
    main()
