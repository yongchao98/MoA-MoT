def translate_first_protein():
    """
    Finds the first Open Reading Frame (ORF) in a DNA sequence and translates it
    into a protein sequence using the single-letter amino acid code.
    """
    # The forward DNA sequence in 5' to 3' direction
    dna_seq = "AGTTGCTTCCCTACCAATAGCTTCTTTTGCCATGGGG"

    # Standard DNA codon table for translation
    codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_STOP_', 'TAG':'_STOP_',
        'TGC':'C', 'TGT':'C', 'TGA':'_STOP_', 'TGG':'W',
    }

    protein_seq = ""
    start_codon = "ATG"

    # Find the starting position of the first start codon
    start_index = dna_seq.find(start_codon)

    if start_index == -1:
        print("No start codon (ATG) found in the sequence.")
        return

    print("Found start codon at the end of the DNA sequence.")
    print("Translating the following codons:")
    
    # Iterate through the sequence by codons (3 bases at a time) from the start codon
    for i in range(start_index, len(dna_seq), 3):
        codon = dna_seq[i:i+3]
        
        # Ensure we have a full codon
        if len(codon) == 3:
            amino_acid = codon_table.get(codon, "X") # 'X' for unknown codons
            
            # Stop translation if a stop codon is found
            if amino_acid == "_STOP_":
                print(f"{codon} -> STOP")
                break
            
            # Add the amino acid to the protein sequence
            protein_seq += amino_acid
            print(f"{codon} -> {amino_acid}")
        else:
            print("\nReached the end of the DNA sequence.")

    print("\nFinal Amino Acid Sequence:")
    print(protein_seq)

# Execute the function
translate_first_protein()
<<<MG>>>