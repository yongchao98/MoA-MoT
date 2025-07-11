def translate_dna_to_protein(dna_sequence):
    """
    Translates a given 5' to 3' DNA sequence into its first potential protein sequence.
    """
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
        'UAC': 'Y', 'UAU': 'Y', 'UAA': '*', 'UAG': '*',
        'UGC': 'C', 'UGU': 'C', 'UGA': '*', 'UGG': 'W',
    }

    # Step 1: Transcribe DNA to mRNA
    mRNA_sequence = dna_sequence.replace('T', 'U')

    # Step 2: Find the start codon 'AUG'
    start_index = mRNA_sequence.find('AUG')
    
    protein_sequence = ""
    
    if start_index != -1:
        # Step 3 & 4: Translate the sequence from the start codon
        for i in range(start_index, len(mRNA_sequence), 3):
            codon = mRNA_sequence[i:i+3]
            
            # Check if there's a full codon left
            if len(codon) < 3:
                break
            
            amino_acid = codon_table.get(codon, '')
            
            # Stop translation at a stop codon
            if amino_acid == '*':
                break
            
            protein_sequence += amino_acid

    # Step 5: Print the final amino acid sequence
    print(protein_sequence)

# The provided forward DNA sequence
dna = "AGTTGCTTCCCTACCAATAGCTTCTTTTGCCATGGGG"
translate_dna_to_protein(dna)