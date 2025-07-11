def translate_dna(dna_sequence):
    """
    Translates a DNA sequence into an amino acid sequence of the first protein.
    """
    # Standard codon to single-letter amino acid table
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

    protein_sequence = ""
    is_translating = False
    
    # Step 2 & 3: Find start codon and translate
    # Find the index of the first start codon 'AUG'
    try:
        start_index = mrna_sequence.index('AUG')
    except ValueError:
        print("No start codon (AUG) found in the sequence.")
        return

    # Iterate through the mRNA sequence by codons from the start codon
    for i in range(start_index, len(mrna_sequence), 3):
        # Ensure we have a full 3-base codon
        if i + 3 <= len(mrna_sequence):
            codon = mrna_sequence[i:i+3]
            amino_acid = codon_table.get(codon, '?')  # '?' for unknown codons

            # Step 4: Stop translation at a stop codon
            if amino_acid == '_':
                break
            
            # Step 5: Append the amino acid to the protein sequence
            protein_sequence += amino_acid
        else:
            # Reached the end of the mRNA with an incomplete codon
            break

    if protein_sequence:
        print(f"The amino acid sequence for the first translated protein is:")
        print(protein_sequence)
    else:
        print("Translation started but no full codons could be translated.")


# The provided forward DNA sequence
dna = "AGTTGCTTCCCTACCAATAGCTTCTTTTGCCATGGGG"
translate_dna(dna)
