def translate_dna(dna_sequence):
    """
    Translates a DNA sequence to a protein sequence starting from the first start codon.
    """
    # Standard DNA codon table for translation to single-letter amino acids
    codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'AGA':'R', 'AGG':'R', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'TTA':'L', 'TTG':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'GAA':'E', 'GAG':'E', 'GAC':'D', 'GAT':'D',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GGC':'G', 'GGG':'G', 'GGT':'G', 'GGU':'G',
        'TTC':'F', 'TTT':'F', 'TAC':'Y', 'TAT':'Y',
        'TGC':'C', 'TGT':'C', 'TGG':'W',
        'TAA':'_stop_', 'TAG':'_stop_', 'TGA':'_stop_'
    }

    protein_sequence = ""
    # Find the starting position of the first start codon 'ATG'
    start_codon_pos = dna_sequence.find('ATG')

    if start_codon_pos != -1:
        # Start translation from the start codon
        for i in range(start_codon_pos, len(dna_sequence), 3):
            codon = dna_sequence[i:i+3]
            # Ensure we have a full codon
            if len(codon) == 3:
                amino_acid = codon_table.get(codon, '?') # Use '?' for unknown codons
                # Stop translation if a stop codon is found
                if amino_acid == '_stop_':
                    break
                protein_sequence += amino_acid
            else:
                # End of sequence before a full codon
                break
    
    if protein_sequence:
        print("DNA sequence: {}".format(dna_sequence))
        print("Reading frame starts at: {}".format(dna_sequence[start_codon_pos:]))
        print("Translated Amino Acid Sequence (single-letter code):")
        print(protein_sequence)
    else:
        print("No start codon 'ATG' found in the provided DNA sequence.")

# Provided forward DNA sequence
dna = "AGTTGCTTCCCTACCAATAGCTTCTTTTGCCATGGGG"
translate_dna(dna)