def translate_dna_to_protein(dna_sequence):
    """
    Finds the first protein sequence by transcribing and translating a given DNA sequence.
    """
    # The standard genetic code (RNA codon to amino acid)
    genetic_code = {
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L', 'UCU': 'S', 'UCC': 'S',
        'UCA': 'S', 'UCG': 'S', 'UAU': 'Y', 'UAC': 'Y', 'UAA': 'Stop', 'UAG': 'Stop',
        'UGU': 'C', 'UGC': 'C', 'UGA': 'Stop', 'UGG': 'W', 'CUU': 'L', 'CUC': 'L',
        'CUA': 'L', 'CUG': 'L', 'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGU': 'R', 'CGC': 'R',
        'CGA': 'R', 'CGG': 'R', 'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAU': 'N', 'AAC': 'N',
        'AAA': 'K', 'AAG': 'K', 'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V', 'GCU': 'A', 'GCC': 'A',
        'GCA': 'A', 'GCG': 'A', 'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }

    # Step 1: Transcribe DNA to mRNA
    mRNA_sequence = dna_sequence.replace('T', 'U')

    # Step 2: Find the start codon 'AUG'
    start_codon_index = -1
    for i in range(len(mRNA_sequence)):
        if mRNA_sequence[i:i+3] == 'AUG':
            start_codon_index = i
            break

    protein_sequence = ""
    # Step 3 & 4: Translate mRNA to protein if start codon is found
    if start_codon_index != -1:
        # Start translation from the found start codon
        for i in range(start_codon_index, len(mRNA_sequence), 3):
            codon = mRNA_sequence[i:i+3]
            # Ensure the codon is 3 characters long
            if len(codon) < 3:
                break
            
            amino_acid = genetic_code.get(codon, '')
            # Step 5: Stop translation at a stop codon
            if amino_acid == 'Stop':
                break
            
            protein_sequence += amino_acid

    # Step 6: Print the final amino acid sequence
    if protein_sequence:
        print(protein_sequence)
    else:
        print("No protein could be translated (no start codon found).")

# Provided DNA sequence
dna = "AGTTGCTTCCCTACCAATAGCTTCTTTTGCCATGGGG"
translate_dna_to_protein(dna)