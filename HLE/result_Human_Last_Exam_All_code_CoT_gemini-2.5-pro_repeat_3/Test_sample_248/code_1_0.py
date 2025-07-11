def translate_dna(dna_sequence):
    """
    Transcribes and translates a DNA sequence to find the first protein.

    Args:
        dna_sequence (str): The 5' to 3' forward DNA sequence.

    Returns:
        str: The amino acid sequence in single-letter code.
    """
    codon_table = {
        'AUA':'I', 'AUC':'I', 'AUU':'I', 'AUG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACU':'T',
        'AAC':'N', 'AAU':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGU':'S', 'AGA':'R', 'AGG':'R',
        'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P',
        'CAC':'H', 'CAU':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGU':'R',
        'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A',
        'GAC':'D', 'GAU':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G',
        'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S',
        'UUC':'F', 'UUU':'F', 'UUA':'L', 'UUG':'L',
        'UAC':'Y', 'UAU':'Y', 'UAA':'_', 'UAG':'_',
        'UGC':'C', 'UGU':'C', 'UGA':'_', 'UGG':'W',
    }

    # Step 1: Transcribe DNA to mRNA
    mRNA = dna_sequence.replace('T', 'U')

    # Step 2: Find the first start codon 'AUG'
    try:
        start_index = mRNA.index('AUG')
    except ValueError:
        # No start codon found
        print("No start codon (AUG) found in the sequence.")
        return ""

    protein_sequence = ""
    # Step 3 & 4: Translate the mRNA from the start codon
    for i in range(start_index, len(mRNA), 3):
        # Ensure we have a full codon
        if i + 3 > len(mRNA):
            break

        codon = mRNA[i:i+3]
        amino_acid = codon_table.get(codon, '?') # Use '?' for unknown codons

        # Stop translation if a stop codon is found
        if amino_acid == '_':
            break

        protein_sequence += amino_acid

    print(f"The DNA sequence is: {dna_sequence}")
    print(f"The transcribed mRNA is: {mRNA}")
    print(f"Translation starts at the first 'AUG' codon.")
    print(f"The resulting amino acid sequence is:")
    print(protein_sequence)

# The provided DNA sequence
dna = "AGTTGCTTCCCTACCAATAGCTTCTTTTGCCATGGGG"
translate_dna(dna)