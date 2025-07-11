def translate_dna_to_protein():
    """
    Finds the first protein sequence translated from a given DNA sequence.
    """
    # The forward DNA sequence in the 5' to 3' direction
    dna_sequence = "AGTTGCTTCCCTACCAATAGCTTCTTTTGCCATGGGG"

    # Standard genetic code dictionary (RNA codon -> single-letter amino acid)
    # '_' represents a stop codon.
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
        'UAC':'Y', 'UAU':'Y', 'UAA':'_', 'UAG':'_', 'UGA':'_'
    }
    stop_codons = ['UAA', 'UAG', 'UGA']

    print(f"Original 5' to 3' DNA Sequence:\n{dna_sequence}\n")

    # Step 1: Transcribe DNA to mRNA by replacing 'T' with 'U'
    mrna_sequence = dna_sequence.replace('T', 'U')
    print(f"Transcribed 5' to 3' mRNA Sequence:\n{mrna_sequence}\n")

    # Step 2: Find the start codon 'AUG'
    start_codon_index = mrna_sequence.find('AUG')

    protein_sequence = ""
    if start_codon_index != -1:
        print(f"Start codon 'AUG' found. Translation begins from this point.\n")
        
        # Step 3 & 4: Translate the mRNA sequence from the start codon
        codons_translated = []
        for i in range(start_codon_index, len(mrna_sequence), 3):
            codon = mrna_sequence[i:i+3]
            
            # Ensure we have a full 3-base codon
            if len(codon) == 3:
                codons_translated.append(codon)
                if codon in stop_codons:
                    break  # Stop translation at a stop codon
                
                amino_acid = genetic_code.get(codon, '?') # Use '?' for unknown codons
                protein_sequence += amino_acid
        
        print(f"Codons translated: {' -> '.join(codons_translated)}")
        print(f"Amino acids:       {' -> '.join(list(protein_sequence))}\n")

    else:
        protein_sequence = "No start codon 'AUG' found, so no protein is translated."

    # Step 5: Print the final amino acid sequence
    print("Final Amino Acid Sequence (Single-Letter Code):")
    print(protein_sequence)

# Execute the function
translate_dna_to_protein()