def translate_dna():
    """
    Translates a given DNA sequence to find the first protein sequence.
    """
    dna_seq = "AGTTGCTTCCCTACCAATAGCTTCTTTTGCCATGGGG"

    # The standard genetic code (RNA codon to amino acid)
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

    # Step 1: Transcribe DNA to mRNA
    mrna_seq = dna_seq.replace('T', 'U')

    # Step 2: Find the start codon 'AUG'
    start_codon_index = -1
    for i in range(len(mrna_seq)):
        if mrna_seq[i:i+3] == 'AUG':
            start_codon_index = i
            break

    protein_sequence = []
    if start_codon_index != -1:
        print(f"Start codon 'AUG' found. Translation begins.")
        # Step 3 & 4: Translate the mRNA from the start codon
        for i in range(start_codon_index, len(mrna_seq), 3):
            codon = mrna_seq[i:i+3]
            # Ensure codon is complete (3 bases)
            if len(codon) == 3:
                amino_acid = codon_table.get(codon, '?') # Use '?' for unknown codons
                if amino_acid == '_STOP_':
                    print(f"{codon} -> STOP")
                    break
                else:
                    protein_sequence.append(amino_acid)
                    print(f"{codon} -> {amino_acid}")
            else:
                break # Reached the end of the sequence
    else:
        print("No start codon 'AUG' found in the mRNA sequence.")
        return

    # Step 5: Print the final sequence
    final_protein = "".join(protein_sequence)
    print(f"\nFinal protein sequence: {final_protein}")
    print(f"<<<{final_protein}>>>")


translate_dna()