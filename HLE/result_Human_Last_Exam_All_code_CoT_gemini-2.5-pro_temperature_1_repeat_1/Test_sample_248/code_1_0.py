def translate_first_protein(dna_sequence):
    """
    Finds the amino acid sequence of the first protein translated from a given DNA coding strand.

    Args:
        dna_sequence (str): The 5' to 3' DNA coding strand.

    Returns:
        str: The amino acid sequence in single-letter format. Returns an error message
             if no start codon is found.
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
        'UAC':'Y', 'UAU':'Y', 'UAA':'*', 'UAG':'*',
        'UGC':'C', 'UGU':'C', 'UGA':'*', 'UGG':'W',
    }

    # Step 1: Transcribe DNA to mRNA
    mrna_sequence = dna_sequence.replace('T', 'U')

    # Step 2: Find the first start codon 'AUG'
    start_index = mrna_sequence.find('AUG')

    if start_index == -1:
        return "No start codon (AUG) found in the sequence."

    protein_sequence = ""
    # Step 3 & 4: Translate the reading frame from the start codon
    for i in range(start_index, len(mrna_sequence), 3):
        codon = mrna_sequence[i:i+3]
        
        # Stop if the codon is incomplete
        if len(codon) < 3:
            break
            
        amino_acid = codon_table.get(codon)
        
        # Stop if a stop codon is found
        if amino_acid == '*':
            break
            
        protein_sequence += amino_acid

    return protein_sequence

# The provided forward DNA sequence
dna = "AGTTGCTTCCCTACCAATAGCTTCTTTTGCCATGGGG"

# Get the amino acid sequence
amino_acid_sequence = translate_first_protein(dna)

# Print the final result
print(amino_acid_sequence)