def find_protein_sequence(dna_sequence):
    """
    This function finds the first protein sequence translated from a given DNA sequence.
    It transcribes the DNA to mRNA, finds the first start codon (AUG), and translates
    the subsequent codons until a stop codon or the end of the sequence is reached.
    """
    genetic_code = {
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
        'UAC':'Y', 'UAU':'Y', 'UAA':'_', 'UAG':'_', 'UGA':'_',
        'UGC':'C', 'UGU':'C', 'UGG':'W',
    }
    stop_codons = ['UAA', 'UAG', 'UGA']
    protein_sequence = ""

    # Step 1: Transcribe the DNA sequence to mRNA
    mrna = dna_sequence.upper().replace('T', 'U')

    # Step 2: Find the position of the first start codon 'AUG'
    try:
        start_index = mrna.index('AUG')
    except ValueError:
        print("No start codon (AUG) found in the sequence.")
        return

    # Step 3 & 4: Translate the mRNA from the start codon
    for i in range(start_index, len(mrna), 3):
        codon = mrna[i:i+3]
        
        # Stop if the codon is incomplete
        if len(codon) < 3:
            break
        
        # Stop if a stop codon is found
        if codon in stop_codons:
            break
            
        # Translate the codon to an amino acid
        amino_acid = genetic_code.get(codon, '?') # '?' for unknown codons
        protein_sequence += amino_acid

    print(protein_sequence)

# Provided DNA sequence
dna = "AGTTGCTTCCCTACCAATAGCTTCTTTTGCCATGGGG"
find_protein_sequence(dna)