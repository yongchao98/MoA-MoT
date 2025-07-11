def translate_dna_to_protein(dna_sequence):
    """
    Translates a DNA coding strand into an amino acid sequence.

    Args:
        dna_sequence (str): The 5' to 3' DNA coding strand.

    Returns:
        str: The translated amino acid sequence in single-letter format.
             Returns an empty string if no start codon is found.
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
    stop_codons = {'TAA', 'TAG', 'TGA'}
    protein_sequence = ""
    start_codon_pos = dna_sequence.find('ATG')

    if start_codon_pos != -1:
        # Start translation from the found start codon
        coding_region = dna_sequence[start_codon_pos:]
        
        # Process the sequence in 3-base codons
        for i in range(0, len(coding_region), 3):
            codon = coding_region[i:i+3]
            
            # Ensure we have a full codon
            if len(codon) == 3:
                if codon in stop_codons:
                    break  # Stop translation at a stop codon
                
                # Translate codon to amino acid
                amino_acid = codon_table.get(codon, '?') # '?' for unknown codons
                protein_sequence += amino_acid

    print(f"The DNA sequence is: {dna_sequence}")
    if protein_sequence:
        print(f"The first protein sequence is: {protein_sequence}")
    else:
        print("No start codon (ATG) was found in the provided DNA sequence.")

# The provided forward DNA sequence
dna = "AGTTGCTTCCCTACCAATAGCTTCTTTTGCCATGGGG"
translate_dna_to_protein(dna)