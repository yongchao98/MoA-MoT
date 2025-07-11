def find_protein_sequence():
    """
    Finds and prints the amino acid sequence of the first protein translated
    from a given DNA sequence.
    """
    # The provided 5' to 3' forward DNA sequence
    dna_sequence = "AGTTGCTTCCCTACCAATAGCTTCTTTTGCCATGGGG"

    # Standard genetic code mapping (RNA codons to single-letter amino acid code)
    # '*' represents a stop codon.
    genetic_code = {
        'AUA': 'I', 'AUC': 'I', 'AUU': 'I', 'AUG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
        'AAC': 'N', 'AAU': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGU': 'S', 'AGA': 'R', 'AGG': 'R',
        'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
        'CAC': 'H', 'CAU': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',
        'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
        'GAC': 'D', 'GAU': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
        'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S',
        'UUC': 'F', 'UUU': 'F', 'UUA': 'L', 'UUG': 'L',
        'UAC': 'Y', 'UAU': 'Y', 'UAA': '*', 'UAG': '*',
        'UGC': 'C', 'UGU': 'C', 'UGA': '*', 'UGG': 'W',
    }

    # Step 1: Transcribe the DNA coding strand to mRNA by replacing T with U.
    mrna_sequence = dna_sequence.replace('T', 'U')

    # Step 2: Find the index of the first start codon 'AUG'.
    start_index = mrna_sequence.find('AUG')

    protein_sequence = ""
    # Proceed only if a start codon is found.
    if start_index != -1:
        # Step 3 & 4: Translate the mRNA sequence starting from the start codon.
        # Iterate through the sequence in steps of 3 (for each codon).
        for i in range(start_index, len(mrna_sequence), 3):
            # Check if there is a complete codon of 3 bases left.
            if i + 3 <= len(mrna_sequence):
                codon = mrna_sequence[i:i+3]
                amino_acid = genetic_code.get(codon, '') # Look up the codon.

                # Step 5: Stop translation if a stop codon is found.
                if amino_acid == '*':
                    break
                
                protein_sequence += amino_acid
            else:
                # Break the loop if an incomplete codon is at the end.
                break

    # Print the resulting protein sequence.
    if protein_sequence:
        print(protein_sequence)
    else:
        print("No protein could be translated (no start codon found).")

# Execute the function to find and print the protein sequence.
find_protein_sequence()