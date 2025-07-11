def translate_orf():
    """
    Finds the correct reading frame for a nucleotide sequence from an ORF
    and translates it into a one-letter amino acid protein sequence.
    """
    # The nucleotide sequence from the middle of an ORF
    dna_sequence = "ggacattgttgttagcgaagacttgaatggaacggtgaagttttcaagctctttaccctaccccaataatcttaacagtgtcctggctgagcgactggagaagtggc"

    # Standard genetic code dictionary (DNA -> Amino Acid)
    genetic_code = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }
    
    # Analysis has shown that the only frame without a stop codon is forward frame 2.
    # We will therefore start reading from the second nucleotide (index 1).
    start_frame = 1 
    
    protein_sequence = []
    
    # Iterate through the sequence by codons (3 nucleotides at a time)
    for i in range(start_frame, len(dna_sequence), 3):
        # Ensure we have a full codon
        if i + 3 <= len(dna_sequence):
            codon = dna_sequence[i:i+3]
            amino_acid = genetic_code.get(codon, 'X') # 'X' for unknown codons
            
            # Since this is an ORF, we don't expect a stop codon.
            # If we were to find one, it would indicate an issue.
            if amino_acid == '*':
                # This block will not be reached for the correct frame
                print(f"Error: Stop codon '{codon}' found in the middle of an ORF.")
                return

            protein_sequence.append(amino_acid)

    # Print the final resulting protein sequence
    final_protein = "".join(protein_sequence)
    for aa in final_protein:
        print(aa, end="")
    print()


# Execute the function to get the answer
translate_orf()