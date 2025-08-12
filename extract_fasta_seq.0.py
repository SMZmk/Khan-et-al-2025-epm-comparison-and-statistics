import sys

def extract_fasta_sequences(fasta_file, selection_list_file, output_file):
    """
    Extracts sequences from a multi-FASTA file based on a list of headers.

    Args:
        fasta_file (str): Path to the input multi-FASTA file.
        selection_list_file (str): Path to the file containing the list of headers.
        output_file (str): Path to the output FASTA file.
    """
    try:
        with open(selection_list_file, 'r') as f:
            selected_headers = [line.strip() for line in f if line.strip()]
    except FileNotFoundError:
        print(f"Error: Selection list file not found at {selection_list_file}", file=sys.stderr)
        return
    except Exception as e:
        print(f"Error reading selection list file: {e}", file=sys.stderr)
        return

    fasta_sequences = {}
    current_header = None
    current_sequence = []

    try:
        with open(fasta_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if current_header:
                        fasta_sequences[current_header] = "".join(current_sequence)
                    current_header = line[1:].split()[0] # Get header, strip '>' and potential description
                    current_sequence = []
                else:
                    current_sequence.append(line)
            # Add the last sequence
            if current_header:
                fasta_sequences[current_header] = "".join(current_sequence)
    except FileNotFoundError:
        print(f"Error: FASTA file not found at {fasta_file}", file=sys.stderr)
        return
    except Exception as e:
        print(f"Error reading FASTA file: {e}", file=sys.stderr)
        return


    try:
        with open(output_file, 'w') as f:
            for header in selected_headers:
                if header in fasta_sequences:
                    f.write(f">{header}\n")
                    # Write sequence in lines of 60 characters (standard FASTA format)
                    sequence = fasta_sequences[header]
                    for i in range(0, len(sequence), 60):
                        f.write(sequence[i:i+60] + '\n')
                else:
                    print(f"Warning: Header '{header}' not found in FASTA file.", file=sys.stderr)
    except Exception as e:
        print(f"Error writing output file: {e}", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python extract_fasta.py <fasta_file> <selection_list_file> <output_file>")
        sys.exit(1)

    fasta_input = sys.argv[1]
    selection_list_input = sys.argv[2]
    output_file_name = sys.argv[3]

    extract_fasta_sequences(fasta_input, selection_list_input, output_file_name)
