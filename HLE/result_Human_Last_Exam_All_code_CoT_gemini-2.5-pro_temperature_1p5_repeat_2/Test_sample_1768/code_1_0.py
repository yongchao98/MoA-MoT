def translate_orf():
    """
    Finds the correct reading frame for a nucleotide sequence from the middle of an ORF
    and translates it into a protein sequence.
    """
    # The nucleotide sequence provided by the user
    dna_seq = "ggacattgttgttagcgaagacttgaatggaacggtgaagttttcaagctctttaccctaccccaataatcttaacagtgtcctggctgagcgactggagaagtggc"
    dna_seq = dna_seq.upper()

    # Standard DNA to one-letter amino acid code dictionary
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

    # Iterate through the three possible reading frames
    for frame in range(3):
        protein_seq = []
        codon_list = []
        has_stop_codon = False
        
        # Get the sequence for the current frame
        frame_seq = dna_seq[frame:]
        
        # Process the sequence in 3-base codons
        for i in range(0, len(frame_seq) - len(frame_seq) % 3, 3):
            codon = frame_seq[i:i+3]
            amino_acid = genetic_code.get(codon, 'X') # 'X' for unknown codon
            
            # Check for stop codon
            if amino_acid == '*':
                has_stop_codon = True
                break
            
            codon_list.append(codon)
            protein_seq.append(amino_acid)

        # If the entire frame was translated without a stop codon, it's the correct one
        if not has_stop_codon:
            final_protein = "".join(protein_seq)
            print(f"The correct reading frame is Frame {frame + 1}.")
            print(f"The protein sequence is: {final_protein}\n")
            print("Translation Details (DNA -> Protein):")
            
            # Print the final equation: codons and their corresponding amino acids
            print("DNA:    ", ' '.join(codon_list))
            print("Protein:", '  '.join(protein_seq))
            
            # Return the protein sequence as the final answer
            return final_protein

    # If no valid ORF is found
    print("No open reading frame without a stop codon was found in any of the three frames.")
    return None

# Execute the function and capture the result
final_protein_sequence = translate_orf()
# The final answer is the protein sequence from the correct frame
# We add the <<<>>> marker to the captured result
# However, the instructions say to use print. 
# The function above prints the answer. Here we just call it.
# To satisfy the format, we can print the answer marker here.
print(f"\n<<<{final_protein_sequence}>>>")