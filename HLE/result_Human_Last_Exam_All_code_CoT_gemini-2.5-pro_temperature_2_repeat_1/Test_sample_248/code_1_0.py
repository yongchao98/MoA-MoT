def translate_first_protein(dna_sequence):
    """
    Translates the first open reading frame from a given 5' to 3' DNA sequence.

    Args:
        dna_sequence (str): The forward DNA sequence in the 5' to 3' direction.

    Returns:
        str: The amino acid sequence in single-letter code.
    """
    # Standard RNA codon table mapping codons to single-letter amino acids
    # '_' represents a stop codon.
    codon_table = {
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
        'UAC': 'Y', 'UAU': 'Y', 'UAA': '_', 'UAG': '_', 'UGA': '_'
    }

    # Step 1: Transcribe DNA to mRNA
    mrna_sequence = dna_sequence.upper().replace('T', 'U')

    # Step 2: Find the first start codon 'AUG'
    start_codon_index = mrna_sequence.find('AUG')

    protein_sequence = ""
    # If no start codon is found, return an empty string
    if start_codon_index == -1:
        print("No start codon (AUG) found.")
        return

    # Step 3, 4 & 5: Translate from the start codon
    # Iterate through the mRNA in 3-nucleotide steps (codons)
    for i in range(start_codon_index, len(mrna_sequence), 3):
        codon = mrna_sequence[i:i+3]

        # Check if a full codon is available
        if len(codon) < 3:
            break

        # Look up the amino acid for the codon
        amino_acid = codon_table.get(codon, '') # Default to empty string if not found

        # If a stop codon is found, end translation
        if amino_acid == '_':
            break
        
        protein_sequence += amino_acid
    
    # Print the resulting amino acid sequence
    if protein_sequence:
        print(f"The resulting amino acid sequence is:")
        print(protein_sequence)
    else:
        print("Translation started but no full codons could be read.")


# The provided forward DNA sequence
dna = "AGTTGCTTCCCTACCAATAGCTTCTTTTGCCATGGGG"
translate_first_protein(dna)
