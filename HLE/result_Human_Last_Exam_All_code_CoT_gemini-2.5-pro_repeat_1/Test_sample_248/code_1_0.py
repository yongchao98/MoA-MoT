def translate_dna_to_protein(dna_sequence):
    """
    Translates a given 5' to 3' DNA sequence into an amino acid sequence.
    """
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
        'UAC': 'Y', 'UAU': 'Y', 'UAA': '_', 'UAG': '_', 'UGA': '_'
    }

    # Step 1: Transcribe DNA to mRNA
    rna_sequence = dna_sequence.replace('T', 'U')

    # Step 2: Find the first start codon 'AUG'
    start_index = rna_sequence.find('AUG')

    protein_sequence = ""
    if start_index != -1:
        # Step 3 & 4: Translate from start codon until a stop codon is found
        # or the sequence ends.
        for i in range(start_index, len(rna_sequence), 3):
            codon = rna_sequence[i:i+3]
            if len(codon) < 3:
                break  # Incomplete codon at the end

            amino_acid = genetic_code.get(codon, '?')
            if amino_acid == '_':
                break  # Stop codon reached
            protein_sequence += amino_acid
    
    print(f"The final amino acid sequence is: {protein_sequence}")


# Provided DNA sequence
dna = "AGTTGCTTCCCTACCAATAGCTTCTTTTGCCATGGGG"
translate_dna_to_protein(dna)