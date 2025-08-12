#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import os

def create_padded_alignment(ref_fasta_file, pattern_fasta_file, map_file, output_fasta_file):
    """
    Creates a new FASTA alignment by padding patterns according to a map file
    (PatternID\tStartPosition\tStrand) and combining them with reference sequences.

    Args:
        ref_fasta_file (str): Path to the reference alignment FASTA (file1.fas).
        pattern_fasta_file (str): Path to the patterns FASTA (file2.fas).
        map_file (str): Path to the mapping table (file3.txt / alingment_index.tsv).
                        Expected format: PatternID\tStartPosition\tStrand (+ or -).
        output_fasta_file (str): Path for the combined output FASTA file.
    """
    print(f"Starting alignment creation process...", file=sys.stderr)

    # --- 1. Read Reference Sequences and Find Max Length ---
    ref_seqs = {}
    max_len = 0
    print(f"Reading reference sequences from {ref_fasta_file}...", file=sys.stderr)
    try:
        for record in SeqIO.parse(ref_fasta_file, "fasta"):
            ref_seqs[record.id] = record.seq
            if len(record.seq) > max_len:
                max_len = len(record.seq)

        if not ref_seqs:
            print(f"Error: No sequences found in reference file: {ref_fasta_file}", file=sys.stderr)
            sys.exit(1)
        print(f"Found {len(ref_seqs)} reference sequences. Determined maximum alignment length: {max_len}", file=sys.stderr)
    except FileNotFoundError:
        print(f"Error: Reference FASTA file not found: {ref_fasta_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading reference FASTA file '{ref_fasta_file}': {e}", file=sys.stderr)
        sys.exit(1)

    # --- 2. Read Pattern Sequences ---
    pattern_seqs = {}
    print(f"Reading pattern sequences from {pattern_fasta_file}...", file=sys.stderr)
    try:
        for record in SeqIO.parse(pattern_fasta_file, "fasta"):
            if record.id not in pattern_seqs:
                 pattern_seqs[record.id] = record.seq # Store Seq object directly
            else:
                 if str(pattern_seqs[record.id]) != str(record.seq):
                     print(f"Warning: Duplicate pattern ID '{record.id}' found in {pattern_fasta_file} with different sequences. Using the first one encountered.", file=sys.stderr)

        if not pattern_seqs:
            print(f"Warning: No sequences found in pattern file: {pattern_fasta_file}", file=sys.stderr)
    except FileNotFoundError:
        print(f"Error: Pattern FASTA file not found: {pattern_fasta_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading pattern FASTA file '{pattern_fasta_file}': {e}", file=sys.stderr)
        sys.exit(1)
    print(f"Found {len(pattern_seqs)} unique pattern sequences.", file=sys.stderr)

    # --- 3. Prepare Output Sequences and Pattern Counter ---
    output_records = []
    # Add original reference sequences first, padded to max_len if needed
    for ref_id, ref_seq in ref_seqs.items():
        padding_needed = max_len - len(ref_seq)
        if padding_needed < 0:
             print(f"Warning: Reference sequence {ref_id} is longer ({len(ref_seq)}) than determined max_len ({max_len}). Truncating reference sequence.", file=sys.stderr)
             padded_seq = ref_seq[:max_len]
        elif padding_needed > 0:
             padded_seq = ref_seq + Seq("-" * padding_needed)
        else:
             padded_seq = ref_seq
        output_records.append(SeqRecord(padded_seq, id=ref_id, description=""))

    pattern_counts = {}
    padded_pattern_records = []

    # --- 4. Process Map File (alingment_index.tsv) with Strand ---
    print(f"Processing map file {map_file} (expecting 3 columns)...", file=sys.stderr)
    line_num = 0
    patterns_added = 0
    try:
        with open(map_file, 'r') as f_map:
            for line in f_map:
                line_num += 1
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                parts = line.split('\t')
                # *** ADJUSTED: Expect 3 columns ***
                if len(parts) != 3:
                    print(f"Warning: Skipping line {line_num} in {map_file}. Expected 3 tab-separated columns (PatternID StartPosition Strand), got {len(parts)}. Line: '{line}'", file=sys.stderr)
                    continue

                # *** ADJUSTED: Assign parsed parts ***
                pattern_id, start_pos_str, strand_str = parts
                strand = strand_str.strip()

                # Validate Pattern ID
                if pattern_id not in pattern_seqs:
                    print(f"Warning: Skipping line {line_num}. Pattern ID '{pattern_id}' not found in {pattern_fasta_file}.", file=sys.stderr)
                    continue

                # Validate Strand
                if strand not in ['+', '-']:
                    print(f"Warning: Skipping line {line_num}. Invalid strand value '{strand}'. Expected '+' or '-'.", file=sys.stderr)
                    continue

                # Validate Start Position
                try:
                    start_pos = int(start_pos_str)
                    if start_pos < 1:
                         raise ValueError("Start position must be 1 or greater.")
                    num_leading_gaps = start_pos - 1
                except ValueError as e:
                    print(f"Warning: Skipping line {line_num}. Invalid start position '{start_pos_str}': {e}.", file=sys.stderr)
                    continue

                # Check if pattern starts beyond max_len
                if num_leading_gaps >= max_len:
                     print(f"Warning: Skipping line {line_num}. Start position {start_pos} (index {num_leading_gaps}) is at or beyond max alignment length {max_len} for pattern '{pattern_id}'.", file=sys.stderr)
                     continue

                # Retrieve original pattern sequence
                original_pattern_seq = pattern_seqs[pattern_id]

                # *** ADJUSTMENT: Apply reverse complement if strand is '-' ***
                if strand == '-':
                    try:
                        processed_pattern_seq = original_pattern_seq.reverse_complement()
                        # print(f"  Applied reverse complement to {pattern_id} for line {line_num}", file=sys.stderr) # Optional debug
                    except Exception as e:
                         print(f"Warning: Skipping line {line_num}. Could not reverse complement pattern '{pattern_id}': {e}", file=sys.stderr)
                         continue
                else: # strand == '+'
                    processed_pattern_seq = original_pattern_seq

                # Use processed_pattern_seq for subsequent steps
                pattern_len = len(processed_pattern_seq)

                # Create padded sequence string (ensuring uppercase)
                padded_seq_part = "-" * num_leading_gaps + str(processed_pattern_seq).upper()

                # Truncate if pattern extends beyond max_len
                final_padded_len = len(padded_seq_part)
                if final_padded_len > max_len:
                    print(f"Warning: Pattern '{pattern_id}' (Strand: {strand}) starting at {start_pos} extends beyond max length {max_len}. Truncating pattern.", file=sys.stderr)
                    padded_seq_part = padded_seq_part[:max_len]
                    final_padded_len = max_len # Update length after truncation

                # Add trailing gaps
                num_trailing_gaps = max_len - final_padded_len
                padded_seq_final = Seq(padded_seq_part + "-" * num_trailing_gaps) # Create Seq object

                # Generate unique output ID
                count = pattern_counts.get(pattern_id, -1) + 1
                pattern_counts[pattern_id] = count
                output_id = f"{pattern_id}({count})"

                # Add to list of new records, including strand in description
                padded_pattern_records.append(SeqRecord(padded_seq_final, id=output_id, description=f"OriginalID={pattern_id} Start={start_pos} Strand={strand}"))
                patterns_added += 1

    except FileNotFoundError:
        print(f"Error: Map file not found: {map_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error processing map file {map_file} at line ~{line_num}: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Processed {line_num} lines from map file. Added {patterns_added} padded pattern sequences.", file=sys.stderr)

    # Combine reference and new pattern records
    output_records.extend(padded_pattern_records)

    # --- 5. Write Output FASTA ---
    print(f"Writing combined alignment ({len(output_records)} sequences) to {output_fasta_file}...", file=sys.stderr)
    try:
        with open(output_fasta_file, "w") as f_out:
            SeqIO.write(output_records, f_out, "fasta")
    except IOError as e:
        print(f"Error: Could not write to output file {output_fasta_file}: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"\nProcess completed successfully.", file=sys.stderr)
    print(f"Combined alignment saved to: {output_fasta_file}", file=sys.stderr)


# --- Main execution block ---
if __name__ == "__main__":
    # *** ADJUSTED: Update description and help text for map file format ***
    parser = argparse.ArgumentParser(description="Combine a reference FASTA alignment with pattern sequences, padding the patterns based on a 3-column mapping file (PatternID\\tStartPosition\\tStrand). Handles reverse complementation.")

    parser.add_argument("-r", "--reference",
                        required=True,
                        help="Path to the reference alignment FASTA file.")

    parser.add_argument("-p", "--patterns",
                        required=True,
                        help="Path to the pattern FASTA file.")

    parser.add_argument("-m", "--map",
                        required=True,
                        help="Path to the map file. Format: PatternID\\tStartPosition\\tStrand (+ or -). Tab-separated, 1-based start.")

    parser.add_argument("-o", "--output",
                        required=True,
                        help="Path for the combined output FASTA file.")

    args = parser.parse_args()

    # Basic input file existence check
    for file_path in [args.reference, args.patterns, args.map]:
        if not os.path.isfile(file_path):
            print(f"Error: Input file not found: {file_path}", file=sys.stderr)
            sys.exit(1)

    create_padded_alignment(args.reference, args.patterns, args.map, args.output)
