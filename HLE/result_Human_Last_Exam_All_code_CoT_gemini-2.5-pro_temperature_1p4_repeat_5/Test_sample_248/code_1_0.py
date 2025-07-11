def translate_dna_to_protein(dna_sequence):
    """
    Translates a given 5' to 3' forward DNA sequence to its corresponding
    amino acid sequence for the first open reading frame.
    """
    # Codon to single-letter amino acid code mapping
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
    
    # 1. Transcription: Convert DNA to mRNA
    mRNA = dna_sequence.replace('T', 'U')
    
    # 2. Find the first start codon 'AUG'
    start_codon_index = mRNA.find('AUG')
    
    protein_sequence = ""
    
    if start_codon_index != -1:
        # 3. Translate from the start codon
        for i in range(start_codon_index, len(mRNA), 3):
            codon = mRNA[i:i+3]
            # Ensure the codon is complete (3 bases long)
            if len(codon) == 3:
                amino_acid = codon_table.get(codon, '?') # Use '?' for unknown codons
                # 4. Stop translation at a stop codon
                if amino_acid == '_':
                    break
                protein_sequence += amino_acid
            else:
                break # Stop if we reach the end of the sequence with an incomplete codon
                
    print(f"The DNA sequence is: {dna_sequence}")
    print(f"The transcribed mRNA sequence is: {mRNA}")
    
    if protein_sequence:
        print(f"The first start codon 'AUG' was found.")
        print(f"The resulting amino acid sequence is: {protein_sequence}")
    else:
        print("No start codon 'AUG' was found, so no protein was translated.")

# Provided DNA sequence
dna = "AGTTGCTTCCCTACCAATAGCTTCTTTTGCCATGGGG"
translate_dna_to_protein(dna)

<<<MG>>>