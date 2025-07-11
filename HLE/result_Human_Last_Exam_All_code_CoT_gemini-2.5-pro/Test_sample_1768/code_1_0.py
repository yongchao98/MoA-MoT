def find_and_translate_orf():
    """
    Finds the correct reading frame for a nucleotide sequence from the middle of an ORF
    and translates it into a protein sequence.
    """
    # The nucleotide sequence from the middle of an ORF
    dna_sequence = "ggacattgttgttagcgaagacttgaatggaacggtgaagttttcaagctctttaccctaccccaataatcttaacagtgtcctggctgagcgactggagaagtggc"

    # Standard DNA codon table mapping to one-letter amino acid codes
    # '*' represents a stop codon.
    codon_table = {
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
        'AAT': 'N', 'AAC': 'N',
        'GAT': 'D', 'GAC': 'D',
        'TGT': 'C', 'TGC': 'C',
        'CAA': 'Q', 'CAG': 'Q',
        'GAA': 'E', 'GAG': 'E',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
        'CAT': 'H', 'CAC': 'H',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
        'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'AAA': 'K', 'AAG': 'K',
        'ATG': 'M',
        'TTT': 'F', 'TTC': 'F',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'AGT': 'S', 'AGC': 'S',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'TGG': 'W',
        'TAT': 'Y', 'TAC': 'Y',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'TAA': '*', 'TAG': '*', 'TGA': '*'
    }

    print(f"Analyzing sequence: {dna_sequence}\n")
    correct_protein_sequence = ""

    # Iterate through the three possible reading frames (0, 1, 2)
    for frame in range(3):
        print(f"--- Frame {frame + 1} ---")
        protein_sequence = []
        has_stop_codon = False
        
        # Process the sequence in chunks of 3 (codons)
        for i in range(frame, len(dna_sequence), 3):
            codon = dna_sequence[i:i+3]
            
            # Ensure it's a full codon before translating
            if len(codon) == 3:
                amino_acid = codon_table.get(codon, '?') # '?' for unknown codons
                if amino_acid == '*':
                    has_stop_codon = True
                    protein_sequence.append(amino_acid)
                    break  # A stop codon means this is not the correct ORF
                protein_sequence.append(amino_acid)
        
        translated_protein = "".join(protein_sequence)
        print(f"Translation: {translated_protein}")
        
        if has_stop_codon:
            print("Result: Stop codon found. Incorrect frame.\n")
        else:
            print("Result: No stop codon found. This is the likely correct frame.\n")
            correct_protein_sequence = translated_protein

    # Print the final result
    print("="*40)
    print("Final Result:")
    if correct_protein_sequence:
        print("The translated protein sequence from the correct open reading frame is:")
        print(correct_protein_sequence)
    else:
        print("No open reading frame without a stop codon was found.")
    print("="*40)

# Execute the function
find_and_translate_orf()