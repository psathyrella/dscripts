import pandas as pd
import os
import glob
import re

path_to_super = '/fh/fast/matsen_e/shared/bcr-mut-sel/supertest'

loris_output_path = '/fh/fast/matsen_e/shared/bcr-mut-sel/working/data/loris/vanwinkle/'
os.makedirs(loris_output_path, exist_ok=True)

for locus in ['igh', 'igk', 'igl']:
    print(f"\nProcessing locus: {locus}")
    super_file = f"{path_to_super}/vanwinkle-170_{locus}_super.csv.gz"
    if not os.path.exists(super_file):
        print(f"Warning: Super file for locus {locus} does not exist: {super_file}")
    all_locus_data = super_file

    # First, get unique sample IDs from the file
    sample_ids = set()
    chunk_size = 50000
    for chunk in pd.read_csv(all_locus_data, chunksize=chunk_size, dtype={'sample_id': str}, low_memory=False):
        sample_ids.update(chunk['sample_id'].unique())
        
    sample_ids = sorted(list(sample_ids))
    print(f"Found {len(sample_ids)} unique sample IDs: {sample_ids}")

    # Read the large file in chunks to manage memory efficiently
    donor_dataframes = {sample_id: [] for sample_id in sample_ids}

    print("Reading and splitting by sample ID...")
    chunk_count = 0
    for chunk in pd.read_csv(all_locus_data, chunksize=chunk_size, dtype={'sample_id': str}, low_memory=False):
        chunk_count += 1
        if chunk_count % 10 == 0:
            print(f"  Processed {chunk_count} chunks...")

        # Split chunk by sample_id
        for sample_id in sample_ids:
            sample_chunk = chunk[chunk['sample_id'] == sample_id].copy()
            if len(sample_chunk) > 0:
                donor_dataframes[sample_id].append(sample_chunk)

    print(f"Finished reading {chunk_count} chunks")

    # Combine chunks for each sample and computationally sorting
    for sample_id in sample_ids:
        print(f"\nProcessing Sample {sample_id}...")
        
        if not donor_dataframes[sample_id]:
            print(f"  No data found for Sample {sample_id}")
            continue
        
        # Combine all chunks for this sample
        donor_df = pd.concat(donor_dataframes[sample_id], ignore_index=True)
        print(f"  All sequences: {len(donor_df)}")

        # Trim sequence columns to length that is a multiple of 3
        if locus == 'igh':
            sequence_columns = [
                'sequence_heavy', 'indel_reversed_seqs_heavy', 'naive_seq_heavy'
            ]
        else:
            sequence_columns = [
                'sequence_light', 'indel_reversed_seqs_light', 'naive_seq_light'
            ]
        
        for col in sequence_columns:
            if col in donor_df.columns:
                # Trim to make length a multiple of 3
                donor_df[col] = donor_df[col].astype(str).apply(lambda seq: seq[:len(seq) - len(seq) % 3])
                print(f"  Trimmed {col} to multiple of 3 length")
            else:
                print(f"  Warning: column {col} not found")
        
        # Computationally sort
        if locus == 'igh':
            naive_df = donor_df[donor_df['n_mutations_heavy'] == 0].copy()
        else:
            naive_df = donor_df[donor_df['n_mutations_light'] == 0].copy()
        print(f"  Naive sequences: {len(naive_df)}")


        # Save naive dataframe if it has data
        if len(naive_df) > 0:
            naive_output_file = os.path.join(loris_output_path, f'vanwinkle-170_{sample_id}_{locus}_naive_comp_sort.csv.gz')
            naive_df.to_csv(naive_output_file, index=False, compression='gzip')

        else:
            print(f"  No naive sequences found for Sample {sample_id}")
        
        # Clear memory by deleting the sample dataframe chunks
        del donor_dataframes[sample_id]

print("\nProcessing complete!")

# Print summary
print("\nSummary:")
for sample_id in sample_ids:
    print(f"Sample {sample_id}:")

    for locus in ['igh', 'igk', 'igl']:
        nomut_file = os.path.join(loris_output_path, f'vanwinkle-170_{sample_id}_{locus}_naive_comp_sort.csv.gz')
    
        nomut_exists = os.path.exists(nomut_file)

        print(f"  No mutations file: {'✓' if nomut_exists else '✗'} {nomut_file}")