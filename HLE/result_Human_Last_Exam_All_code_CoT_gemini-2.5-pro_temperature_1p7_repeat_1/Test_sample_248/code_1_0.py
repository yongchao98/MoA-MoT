def translate_first_protein(dna_seq):
    """
    Finds the first start codon in a DNA sequence and translates the
    subsequent open reading frame into an amino acid sequence.
    """
    # Standard DNA codon table mapping to single-letter amino acid codes
    # '_' is used to represent a STOP codon.
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

    protein_seq = []
    
    # Find the starting position of the first start codon 'ATG'
    start_pos = dna_seq.find('ATG')

    if start_pos == -1:
        print("No start codon 'ATG' found in the sequence.")
        return

    # Extract the sequence starting from the first start codon
    coding_region = dna_seq[start_pos:]
    
    # Iterate through the coding region in steps of 3 (one codon at a time)
    for i in range(0, len(coding_region), 3):
        codon = coding_region[i:i+3]
        
        # Ensure we have a complete codon
        if len(codon) < 3:
            break
            
        # Translate the codon to an amino acid
        amino_acid = codon_table.get(codon)
        
        # If the codon is a stop codon, end translation
        if amino_acid == '_':
            break
        
        # Append the amino acid to our protein sequence list
        protein_seq.append(amino_acid)

    print("The amino acid sequence of the first protein is:")
    print("".join(protein_seq))

# The provided DNA sequence
dna_sequence = "AGTTGCTTCCCTACCAATAGCTTCTTTTGCCATGGGG"
translate_first_protein(dna_sequence)