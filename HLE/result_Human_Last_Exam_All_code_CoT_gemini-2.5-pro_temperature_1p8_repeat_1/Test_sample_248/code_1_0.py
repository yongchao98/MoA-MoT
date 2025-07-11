def translate_first_protein(dna_sequence):
    """
    Finds the first start codon in a DNA sequence and translates the
    open reading frame into a single-letter amino acid sequence.
    """
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
    # Find the starting position of the first start codon 'ATG'
    start_index = dna_sequence.find('ATG')
    
    # Proceed only if a start codon is found
    if start_index != -1:
        # Create the readable frame starting from the start codon
        reading_frame = dna_sequence[start_index:]
        
        # Translate codon by codon (in steps of 3 bases)
        for i in range(0, len(reading_frame), 3):
            # Check if there is a full codon left
            if i + 3 <= len(reading_frame):
                codon = reading_frame[i:i+3]
                amino_acid = codon_table.get(codon)
                
                # If a stop codon is found, terminate translation
                if amino_acid == '_':
                    break
                
                # Append the amino acid to the protein sequence
                if amino_acid:
                    protein_sequence += amino_acid
            else:
                # Incomplete codon at the end, stop translation
                break

    return protein_sequence

# The provided DNA sequence in 5' to 3' direction
dna = "AGTTGCTTCCCTACCAATAGCTTCTTTTGCCATGGGG"

# Perform the translation
amino_acid_sequence = translate_first_protein(dna)

# Print the final result
print("Amino acid sequence of the first protein:")
print(amino_acid_sequence)