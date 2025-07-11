def translate_first_protein(dna_sequence):
    """
    Finds the amino acid sequence of the first protein from a given DNA coding strand.
    """
    # Genetic code table (codon to single-letter amino acid)
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
        'UAC':'Y', 'UAU':'Y', 'UAA':'_', 'UAG':'_', 'UGA':'_' # '_' for STOP codons
    }

    # Step 1: Transcribe DNA to mRNA by replacing T with U
    mrna_sequence = dna_sequence.replace('T', 'U')

    # Step 2: Find the first start codon 'AUG'
    start_index = mrna_sequence.find('AUG')

    protein_sequence = ""
    if start_index != -1:
        # Step 3 & 4: Translate the mRNA sequence from the start codon
        # Iterate through the mRNA in 3-nucleotide steps (codons)
        for i in range(start_index, len(mrna_sequence), 3):
            codon = mrna_sequence[i:i+3]
            
            # Ensure we have a full codon
            if len(codon) == 3:
                amino_acid = codon_table.get(codon, '') # Get amino acid from table
                # If we hit a stop codon, end translation
                if amino_acid == '_':
                    break
                protein_sequence += amino_acid
            else:
                # Incomplete codon at the end of the sequence
                break
    
    # Step 5: Print the final amino acid sequence
    if protein_sequence:
        print(f"The first translated protein sequence is: {protein_sequence}")
    else:
        print("No start codon ('AUG') found. No protein was translated.")

# The provided forward DNA sequence
dna = "AGTTGCTTCCCTACCAATAGCTTCTTTTGCCATGGGG"
translate_first_protein(dna)