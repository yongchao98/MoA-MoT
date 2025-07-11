def translate_dna(dna_sequence):
    """
    Translates a DNA sequence to an amino acid sequence.

    This function simulates transcription and translation to find the first protein
    sequence from a given 5' to 3' forward DNA strand.
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
        'UAC': 'Y', 'UAU': 'Y', 'UAA': '_', 'UAG': '_',
        'UGC': 'C', 'UGU': 'C', 'UGA': '_', 'UGG': 'W',
    }
    
    # Step 1: Transcription (DNA to mRNA)
    mrna = dna_sequence.replace('T', 'U')
    
    protein = ""
    start_codon_found = False
    
    # Step 2: Find the start codon 'AUG'
    try:
        start_index = mrna.index('AUG')
    except ValueError:
        print("No start codon (AUG) found in the sequence.")
        return
        
    # Step 3 & 4: Translate from start codon until stop codon
    for i in range(start_index, len(mrna), 3):
        # Ensure we have a full codon
        if i + 3 <= len(mrna):
            codon = mrna[i:i+3]
            amino_acid = genetic_code.get(codon, '?') # '?' for unknown codons
            
            # Stop translation at a stop codon
            if amino_acid == '_':
                break
                
            protein += amino_acid
    
    # Step 5: Output the result
    if protein:
        print("Translated Amino Acid Sequence (single-letter code):")
        # Output each amino acid in the final sequence
        for amino_acid_code in protein:
            print(amino_acid_code, end="")
        print() # for a newline at the end
    else:
        print("No protein was translated (translation did not complete after start codon).")

# Provided forward DNA sequence
dna = "AGTTGCTTCCCTACCAATAGCTTCTTTTGCCATGGGG"
translate_dna(dna)