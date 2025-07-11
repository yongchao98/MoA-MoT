def translate_dna_to_protein(dna_sequence):
    """
    Translates a given 5' to 3' DNA sequence into an amino acid sequence.

    Args:
        dna_sequence (str): The forward DNA sequence.

    Returns:
        str: The resulting amino acid sequence in single-letter format.
    """
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
        'UAC':'Y', 'UAU':'Y', 'UAA':'_', 'UAG':'_',
        'UGC':'C', 'UGU':'C', 'UGA':'_', 'UGG':'W',
    }
    stop_codons = ['UAA', 'UAG', 'UGA']

    # Step 1: Transcribe DNA to mRNA
    mrna_sequence = dna_sequence.replace('T', 'U')

    # Step 2: Find the start codon 'AUG'
    try:
        start_index = mrna_sequence.index('AUG')
    except ValueError:
        print("No start codon (AUG) found.")
        return ""

    # Step 3 & 4: Translate the mRNA from the start codon
    protein_sequence = []
    # Create the reading frame starting from the start codon
    reading_frame = mrna_sequence[start_index:]
    
    # Iterate through the reading frame by codons (3 bases)
    for i in range(0, len(reading_frame), 3):
        codon = reading_frame[i:i+3]
        
        # Ensure the codon is complete (3 bases long)
        if len(codon) < 3:
            break
            
        # Stop translation if a stop codon is encountered
        if codon in stop_codons:
            break
        
        # Translate codon to amino acid
        amino_acid = codon_table.get(codon, '?') # '?' for unknown codons
        protein_sequence.append(amino_acid)

    # Step 5: Output the final sequence
    final_sequence = "".join(protein_sequence)
    print(f"DNA Sequence: {dna_sequence}")
    print(f"mRNA Sequence: {mrna_sequence}")
    print(f"Start Codon found at position: {start_index}")
    print(f"Reading Frame: {' '.join([reading_frame[i:i+3] for i in range(0, len(reading_frame) - len(reading_frame)%3, 3)])}")
    print(f"Amino Acid Sequence: {final_sequence}")

# Provided DNA sequence
dna = "AGTTGCTTCCCTACCAATAGCTTCTTTTGCCATGGGG"
translate_dna_to_protein(dna)