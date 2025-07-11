def translate_dna(dna_sequence):
    """
    Transcribes a DNA sequence into mRNA and translates it into an amino acid sequence.
    Translation starts at the first 'AUG' start codon.
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
        'UAC':'Y', 'UAU':'Y', 'UAA':'_STOP_', 'UAG':'_STOP_',
        'UGC':'C', 'UGU':'C', 'UGA':'_STOP_', 'UGG':'W',
    }
    
    # Step 1: Transcription (DNA to mRNA)
    mrna_sequence = dna_sequence.replace('T', 'U')
    
    # Step 2: Find the start codon 'AUG'
    start_codon_index = mrna_sequence.find('AUG')
    
    protein_sequence = ""
    
    if start_codon_index != -1:
        # Step 3 & 4: Translation
        # Start reading codons from the start codon
        for i in range(start_codon_index, len(mrna_sequence), 3):
            codon = mrna_sequence[i:i+3]
            
            # Check if there is a full codon to read
            if len(codon) == 3:
                amino_acid = codon_table.get(codon, '')
                # Stop translation at a stop codon
                if amino_acid == '_STOP_':
                    break
                protein_sequence += amino_acid
    
    # Step 5: Print the final sequence
    if protein_sequence:
        print(f"The DNA sequence is: {dna_sequence}")
        print(f"The transcribed mRNA sequence is: {mrna_sequence}")
        print(f"The first start codon 'AUG' was found starting at position {start_codon_index + 1}.")
        print("Translating the open reading frame...")
        print(f"Final Amino Acid Sequence (single-letter code): {protein_sequence}")
    else:
        print("No start codon ('AUG') found in the potential reading frames of this sequence.")


# The forward DNA sequence provided by the user
dna_fwd_strand = "AGTTGCTTCCCTACCAATAGCTTCTTTTGCCATGGGG"
translate_dna(dna_fwd_strand)

<<<MG>>>