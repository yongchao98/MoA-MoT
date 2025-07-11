def translate_dna(dna_sequence):
    """
    Finds the first protein sequence translated from a given DNA coding strand.
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
        'UAC':'Y', 'UAU':'Y', 'UAA':'Stop', 'UAG':'Stop',
        'UGC':'C', 'UGU':'C', 'UGA':'Stop', 'UGG':'W',
    }

    # Step 1: Transcribe DNA to mRNA
    mrna_sequence = dna_sequence.replace('T', 'U')

    # Step 2: Find the start codon 'AUG'
    start_codon_index = mrna_sequence.find('AUG')

    protein_sequence = ""
    if start_codon_index != -1:
        # Step 3 & 4: Translate the sequence starting from the start codon
        for i in range(start_codon_index, len(mrna_sequence), 3):
            codon = mrna_sequence[i:i+3]

            # If there isn't a full codon at the end, stop
            if len(codon) < 3:
                break
            
            amino_acid = codon_table.get(codon, '?')

            # Step 5: Stop translation if a stop codon is found
            if amino_acid == 'Stop':
                break
            
            protein_sequence += amino_acid
            
    # Step 6: Print the result
    if protein_sequence:
        print(f"The amino acid sequence is: {protein_sequence}")
    else:
        print("No start codon ('AUG') found in the provided DNA sequence, so no protein was translated.")

# The provided forward DNA sequence
dna = "AGTTGCTTCCCTACCAATAGCTTCTTTTGCCATGGGG"
translate_dna(dna)