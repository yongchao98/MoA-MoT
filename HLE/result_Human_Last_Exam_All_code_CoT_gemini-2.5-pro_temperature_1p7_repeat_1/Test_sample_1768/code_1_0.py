import textwrap

def find_and_translate_orf():
    """
    Finds the correct reading frame and translates a nucleotide sequence
    from the middle of an ORF.
    """
    # The given nucleotide sequence
    dna_seq = "ggacattgttgttagcgaagacttgaatggaacggtgaagttttcaagctctttaccctaccccaataatcttaacagtgtcctggctgagcgactggagaagtggc"

    # Standard DNA codon to one-letter amino acid code mapping
    codon_map = {
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
    
    print("Analyzing translations for all three reading frames:")
    print("-" * 50)
    
    results = {}
    
    # Iterate through the three possible reading frames (0, 1, 2)
    for frame_start in range(3):
        protein_seq = []
        # Get the sequence for the current frame
        frame_dna = dna_seq[frame_start:]
        
        # Translate codon by codon
        for i in range(0, len(frame_dna) - 2, 3):
            codon = frame_dna[i:i+3]
            amino_acid = codon_map.get(codon, 'X')  # 'X' for unknown codons
            protein_seq.append(amino_acid)
        
        # Join the amino acids into a string
        results[frame_start + 1] = "".join(protein_seq)

    # Print results for all frames
    for frame, protein in results.items():
        print(f"Frame {frame}: {protein}")
        
    # Identify the correct frame (the one without a stop codon '*')
    correct_frame = None
    correct_protein = None
    for frame, protein in results.items():
        if '*' not in protein:
            correct_frame = frame
            correct_protein = protein
            break
            
    print("-" * 50)
    print("\nExplanation:")
    print("The sequence is from the middle of an Open Reading Frame (ORF), which means it should not contain stop codons (*).")
    print(f"Based on this, Frame {correct_frame} is the correct reading frame.")
    
    print("\nFinal Answer:")
    print(f"Frame: {correct_frame}")
    print(f"Protein: {correct_protein}")

# Run the function
find_and_translate_orf()
<<<Frame: 2
Protein: DIVVSEDLNGTVKFSSSLPTPIILNSVLAERLEKW>>>