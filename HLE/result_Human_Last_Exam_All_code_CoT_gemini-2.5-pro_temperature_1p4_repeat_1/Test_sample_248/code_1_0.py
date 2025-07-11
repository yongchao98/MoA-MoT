def translate_first_orf():
    """
    Finds the first Open Reading Frame (ORF) in a DNA sequence,
    translates it into a protein, and prints the process.
    """
    # The forward DNA sequence in 5' to 3' direction
    dna_sequence = "AGTTGCTTCCCTACCAATAGCTTCTTTTGCCATGGGG"

    # Standard DNA codon table mapping codons to single-letter amino acid codes
    codon_table = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L',
        'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'ATT': 'I', 'ATC': 'I',
        'ATA': 'I', 'ATG': 'M', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V',
        'GTG': 'V', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'AGT': 'S', 'AGC': 'S', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
        'CCG': 'P', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'TAT': 'Y',
        'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'TGA': '*', 'CAT': 'H',
        'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'AAT': 'N', 'AAC': 'N',
        'AAA': 'K', 'AAG': 'K', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E',
        'GAG': 'E', 'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CGT': 'R',
        'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }
    stop_codons = {'TAA', 'TAG', 'TGA'}
    start_codon = 'ATG'
    
    protein_sequence = []
    translation_steps = []
    
    # Find the starting position of the first start codon
    start_index = dna_sequence.find(start_codon)

    if start_index != -1:
        print(f"Start codon '{start_codon}' found. Beginning translation.")
        print("-" * 30)
        # Iterate through the sequence in 3-base steps (codons)
        for i in range(start_index, len(dna_sequence), 3):
            # Check if there is a full codon left
            if i + 3 <= len(dna_sequence):
                codon = dna_sequence[i:i+3]
                # If codon is a stop codon, stop translation
                if codon in stop_codons:
                    print(f"Stop codon '{codon}' reached. Translation terminated.")
                    break
                
                # Translate codon to amino acid
                amino_acid = codon_table.get(codon, '?') # '?' for unknown codons
                protein_sequence.append(amino_acid)
                # Store the equation for this step
                print(f"{codon} -> {amino_acid}")

            else:
                # Reached the end of the sequence with an incomplete codon
                print("End of sequence reached. Incomplete final codon ignored.")
                break
        
        final_sequence = "".join(protein_sequence)
        print("-" * 30)
        print("Final Amino Acid Sequence (Single-Letter Code):")
        print(final_sequence)
    else:
        print("No start codon (ATG) found in the provided DNA sequence.")
        final_sequence = ""

    # The final answer is wrapped in <<<...>>>
    # Note: This print statement is for the submission system and will be part of the output.
    print(f"\n<<<{final_sequence}>>>")

# Execute the function
translate_first_orf()