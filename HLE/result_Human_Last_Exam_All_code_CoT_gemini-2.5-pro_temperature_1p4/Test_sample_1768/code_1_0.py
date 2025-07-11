import textwrap

def translate_orf():
    """
    Finds the correct reading frame and translates a nucleotide sequence into protein.
    """
    # The nucleotide sequence from the middle of an ORF
    dna_sequence = "ggacattgttgttagcgaagacttgaatggaacggtgaagttttcaagctctttaccctaccccaataatcttaacagtgtcctggctgagcgactggagaagtggc"
    dna_upper = dna_sequence.upper()

    # Standard genetic code mapping DNA codons to one-letter amino acid codes
    genetic_code = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '_', 'TAG': '_',
        'TGT': 'C', 'TGC': 'C', 'TGA': '_', 'TGG': 'W',
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

    # Frame 2 is the correct reading frame, starting from the second nucleotide (index 1)
    frame = 2
    start_index = frame - 1
    
    codons = []
    protein_sequence = []

    # Iterate through the sequence in steps of 3 from the start index
    for i in range(start_index, len(dna_upper), 3):
        codon = dna_upper[i:i+3]
        if len(codon) == 3:
            amino_acid = genetic_code.get(codon, 'X') # 'X' for unknown codons
            # Check for stop codons
            if amino_acid == '_':
                print(f"Error: Stop codon {codon} found in frame {frame}. This is not the correct ORF.")
                return
            codons.append(codon)
            protein_sequence.append(amino_acid)

    final_protein = "".join(protein_sequence)

    print(f"The correct reading frame is Frame {frame}.")
    print("-" * 30)
    print("Original Sequence (Frame 2):")
    # Wrap text for better readability
    print(textwrap.fill(dna_upper[start_index:], 60))
    print("\nCodons and corresponding Amino Acids:")
    
    # Print codons and amino acids aligned in chunks
    chunk_size = 20
    for i in range(0, len(codons), chunk_size):
        codon_chunk = codons[i:i+chunk_size]
        protein_chunk = protein_sequence[i:i+chunk_size]
        print(" ".join(codon_chunk))
        print("  " + "   ".join(protein_chunk))
        print()
        
    print("-" * 30)
    print("Final Protein Sequence:")
    print(final_protein)


translate_orf()
<<<DIVVSEDLNGEVKFSSSLPYPNNLNSVLAERLEKW>>>