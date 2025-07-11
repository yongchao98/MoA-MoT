def translate_dna(dna_sequence):
    """
    Translates a given 5' to 3' DNA sequence into an amino acid sequence.
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
        'UAC': 'Y', 'UAU': 'Y', 'UAA': 'Stop', 'UAG': 'Stop',
        'UGC': 'C', 'UGU': 'C', 'UGA': 'Stop', 'UGG': 'W',
    }
    
    # Step 1: Transcribe DNA to mRNA
    mRNA = dna_sequence.replace('T', 'U')
    
    # Step 2: Find the first start codon 'AUG'
    try:
        start_index = mRNA.find('AUG')
    except ValueError:
        print("No start codon (AUG) found.")
        return

    if start_index == -1:
        print("No start codon (AUG) found.")
        return

    protein_sequence = []
    
    # Step 3 & 4: Translate the mRNA from the start codon
    for i in range(start_index, len(mRNA), 3):
        codon = mRNA[i:i+3]
        
        # Ensure we have a full codon
        if len(codon) < 3:
            break
            
        amino_acid = codon_table.get(codon, '?') # '?' for unknown codons
        
        # Step 5: Stop at a stop codon
        if amino_acid == 'Stop':
            break
        
        protein_sequence.append(amino_acid)
        
    print("".join(protein_sequence))

# The provided DNA sequence
dna = "AGTTGCTTCCCTACCAATAGCTTCTTTTGCCATGGGG"
translate_dna(dna)