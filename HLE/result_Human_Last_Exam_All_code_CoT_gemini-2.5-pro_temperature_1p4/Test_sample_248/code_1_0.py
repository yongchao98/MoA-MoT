def translate_first_orf():
    """
    Finds the first Open Reading Frame (ORF) in a DNA sequence,
    translates it, and prints the resulting amino acid sequence.
    """
    # The forward DNA sequence (5' to 3')
    dna_sequence = "AGTTGCTTCCCTACCAATAGCTTCTTTTGCCATGGGG"

    # DNA codon to single-letter amino acid code table
    # '_' represents a stop codon.
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
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }

    protein_sequence = ""
    translated_codons = []

    # Find the start of the first protein, which is the first 'ATG' codon.
    start_index = dna_sequence.find('ATG')

    if start_index == -1:
        print("No start codon 'ATG' found in the sequence.")
        return

    # Start the translation from the identified start codon.
    for i in range(start_index, len(dna_sequence), 3):
        codon = dna_sequence[i:i+3]

        # Stop if there isn't a full codon left
        if len(codon) < 3:
            break

        # Translate the codon to an amino acid
        amino_acid = codon_table.get(codon, '?') # '?' for unknown codons

        # Stop translation if a stop codon is found
        if amino_acid == '_':
            break

        protein_sequence += amino_acid
        translated_codons.append(codon)

    print("The DNA sequence is:")
    print(dna_sequence + "\n")
    print(f"The first start codon 'ATG' was found at position {start_index+1}.")
    print("The following codons were translated:")

    # Format the final "equation" showing the codons and resulting amino acids
    codons_str = " ".join(translated_codons)
    amino_acids_str = " ".join(list(protein_sequence))
    print(f"{codons_str} -> {amino_acids_str}\n")
    print("Final amino acid sequence (single-letter code):")
    print(protein_sequence)

translate_first_orf()
<<<MG>>>