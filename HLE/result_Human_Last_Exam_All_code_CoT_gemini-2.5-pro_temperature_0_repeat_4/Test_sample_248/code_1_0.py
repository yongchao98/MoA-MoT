def translate_dna_to_protein(dna_sequence):
    """
    Translates a given DNA sequence to an amino acid sequence for the first protein.

    Args:
        dna_sequence (str): The forward DNA sequence in 5' to 3' direction.

    Returns:
        str: The amino acid sequence in single-letter format.
    """
    # Standard genetic code (codon to single-letter amino acid)
    codon_table = {
        'AUA':'I', 'AUC':'I', 'AUU':'I', 'AUG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACU':'T',
        'AAC':'N', 'AAU':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGU':'S', 'AGA':'R', 'AGG':'R',
        'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P',
        'CAC':'H', 'CAU':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGU':'R',
        'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A',
        'GAC':'D', 'GAU':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G',
        'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S',
        'UUC':'F', 'UUU':'F', 'UUA':'L', 'UUG':'L',
        'UAC':'Y', 'UAU':'Y', 'UAA':'Stop', 'UAG':'Stop',
        'UGC':'C', 'UGU':'C', 'UGA':'Stop', 'UGG':'W',
    }

    # Step 1: Transcribe DNA to mRNA
    rna_sequence = dna_sequence.replace('T', 'U')

    # Step 2: Find the start codon 'AUG'
    start_codon_index = rna_sequence.find('AUG')

    protein_sequence = ""
    if start_codon_index != -1:
        # Step 3 & 4: Translate the mRNA from the start codon
        # Iterate through the sequence in steps of 3 (codons)
        for i in range(start_codon_index, len(rna_sequence), 3):
            codon = rna_sequence[i:i+3]
            
            # Ensure we have a full codon
            if len(codon) < 3:
                break
            
            amino_acid = codon_table.get(codon, '?') # '?' for unknown codons
            
            # Step 5: Stop translation at a stop codon
            if amino_acid == 'Stop':
                break
            
            protein_sequence += amino_acid

    # Step 6: Print the final amino acid sequence
    print(protein_sequence)

# The provided forward DNA sequence
dna = "AGTTGCTTCCCTACCAATAGCTTCTTTTGCCATGGGG"
translate_dna_to_protein(dna)