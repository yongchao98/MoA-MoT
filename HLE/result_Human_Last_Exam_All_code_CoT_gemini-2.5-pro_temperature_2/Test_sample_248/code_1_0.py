def translate_dna_to_protein(dna_sequence):
    """
    Transcribes a DNA coding strand to mRNA and translates it to an amino acid sequence.
    The translation starts from the first start codon (AUG).

    Args:
        dna_sequence (str): The forward DNA sequence in the 5' to 3' direction.

    Returns:
        str: The translated amino acid sequence in single-letter code.
             Returns an empty string if no start codon is found.
    """
    # Standard genetic code: mRNA codon to single-letter amino acid
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
        'UAC': 'Y', 'UAU': 'Y', 'UAA': 'Stop', 'UAG': 'Stop',
        'UGC': 'C', 'UGU': 'C', 'UGA': 'Stop', 'UGG': 'W',
    }

    # Step 1: Transcribe DNA to mRNA by replacing T with U
    mrna_sequence = dna_sequence.replace('T', 'U')

    # Step 2: Find the first start codon 'AUG'
    start_index = mrna_sequence.find('AUG')

    protein_sequence = []
    
    if start_index != -1:
        # Step 3: Translate from the start codon
        # Iterate through the mRNA sequence in steps of 3 (for each codon)
        for i in range(start_index, len(mrna_sequence), 3):
            codon = mrna_sequence[i:i+3]

            # Ensure we have a full codon of 3 nucleotides
            if len(codon) == 3:
                amino_acid = codon_table.get(codon, '') # Look up amino acid
                
                # Step 4: Check for a stop codon
                if amino_acid == 'Stop':
                    break # Stop translation
                
                protein_sequence.append(amino_acid)

    # Step 5: Join the list of amino acids into the final sequence string
    final_sequence = "".join(protein_sequence)
    
    # Print the resulting amino acid sequence
    if final_sequence:
        print("The amino acid sequence is:")
        print(final_sequence)
    else:
        print("No protein could be translated (no start codon found).")

# Provided DNA sequence
dna = "AGTTGCTTCCCTACCAATAGCTTCTTTTGCCATGGGG"
translate_dna_to_protein(dna)
<<<MG>>>