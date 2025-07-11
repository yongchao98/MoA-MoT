def translate_first_protein(dna_sequence):
    """
    Finds the amino acid sequence of the first protein translated from a given DNA sequence.
    """
    genetic_code = {
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
        'UAC':'Y', 'UAU':'Y', 'UAA':'_STOP_', 'UAG':'_STOP_',
        'UGC':'C', 'UGU':'C', 'UGA':'_STOP_', 'UGG':'W',
    }

    # Step 1: Transcribe DNA to mRNA
    mrna_sequence = dna_sequence.replace('T', 'U')

    # Step 2: Find the first start codon 'AUG'
    start_codon_index = mrna_sequence.find('AUG')

    protein_sequence = ""
    # If a start codon is found
    if start_codon_index != -1:
        # Step 3: Start translating from the start codon
        # Create the reading frame
        coding_sequence = mrna_sequence[start_codon_index:]
        
        # Step 4 & 5: Translate codon by codon until a stop codon or the end of the sequence
        for i in range(0, len(coding_sequence), 3):
            codon = coding_sequence[i:i+3]
            
            # Check if we have a full codon
            if len(codon) == 3:
                amino_acid = genetic_code.get(codon, '?') # Use '?' for unknown codons
                
                # Check for stop codon
                if amino_acid == '_STOP_':
                    break
                protein_sequence += amino_acid
            else:
                # Reached the end of the sequence with an incomplete codon
                break
    else:
        # No protein is translated if no start codon is found
        protein_sequence = "No start codon found."

    # Step 6: Print the final amino acid sequence
    print(protein_sequence)

# The provided forward DNA sequence
dna = "AGTTGCTTCCCTACCAATAGCTTCTTTTGCCATGGGG"
translate_first_protein(dna)