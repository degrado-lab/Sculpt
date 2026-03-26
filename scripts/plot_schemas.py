import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np

def main():
    parser = argparse.ArgumentParser(description="Plot schema tracking statistics across generations")
    parser.add_argument('--csv', required=True, help="Path to schema_statistics.csv")
    parser.add_argument('--out_dir', default='.', help="Directory to save the plots")
    parser.add_argument('--top_n', type=int, default=5, help="Number of top schemas to plot")
    args = parser.parse_args()

    df = pd.read_csv(args.csv)
    max_cycle = df['cycle'].max()

    # Find top schemas by Frequency in the last cycle
    final_cycle_df = df[df['cycle'] == max_cycle]
    
    top_freq_schemas = final_cycle_df.sort_values(by='freq_overall', ascending=False).head(args.top_n)['schema'].tolist()
    top_enrich_schemas = final_cycle_df.sort_values(by='enrichment_score', ascending=False).head(args.top_n)['schema'].tolist()

    os.makedirs(args.out_dir, exist_ok=True)

    # Prepare complete cycles index so missing cycles are handled gracefully
    cycles = sorted(df['cycle'].unique())

    # --- Plot Frequency over time ---
    plt.figure(figsize=(10, 6))
    freq_pivot = df.pivot(index='cycle', columns='schema', values='freq_overall').reindex(cycles).fillna(0)
    
    for schema in top_freq_schemas:
        if schema in freq_pivot.columns:
            plt.plot(freq_pivot.index, freq_pivot[schema], marker='o', label=schema)
            
    plt.title(f'Frequency of Top {args.top_n} Schemas (ranked by Final Generation Frequency)')
    plt.xlabel('Generation (Cycle)')
    plt.ylabel('Frequency Overall')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    out_freq = os.path.join(args.out_dir, 'schema_frequency.png')
    plt.savefig(out_freq, dpi=300)
    plt.close()

    # --- Plot Fitness over time ---
    plt.figure(figsize=(10, 6))
    fit_pivot = df.pivot(index='cycle', columns='schema', values='mean_fitness').reindex(cycles)
    
    for schema in top_freq_schemas:
        if schema in fit_pivot.columns:
            plt.plot(fit_pivot.index, fit_pivot[schema], marker='s', label=schema)
            
    plt.title(f'Mean Fitness of Top {args.top_n} Schemas (ranked by Final Gen Freq)')
    plt.xlabel('Generation (Cycle)')
    plt.ylabel('Mean Fitness')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    out_fit = os.path.join(args.out_dir, 'schema_fitness.png')
    plt.savefig(out_fit, dpi=300)
    plt.close()

    # --- Plot Enrichment Score of Top Enriched over time ---
    plt.figure(figsize=(10, 6))
    enrich_pivot = df.pivot(index='cycle', columns='schema', values='enrichment_score').reindex(cycles)
    
    for schema in top_enrich_schemas:
        if schema in enrich_pivot.columns:
            plt.plot(enrich_pivot.index, enrich_pivot[schema], marker='^', label=schema)
            
    plt.title(f'Enrichment Score of Top {args.top_n} Schemas (ranked by Final Gen Enrichment)')
    plt.xlabel('Generation (Cycle)')
    plt.ylabel('Log2(Freq in Top 10% / Freq Overall)')
    plt.axhline(0, color='gray', linestyle='--')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    out_enrich = os.path.join(args.out_dir, 'schema_enrichment.png')
    plt.savefig(out_enrich, dpi=300)
    plt.close()

    print(f"Saved plots to {args.out_dir}")

if __name__ == "__main__":
    main()
