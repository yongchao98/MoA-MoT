def translate_first_protein():
    """
    Finds and translates the first open reading frame in a given DNA sequence.
    """
    # The forward DNA sequence in the 5' to 3' direction
    dna_sequence = "AGTTGCTTCCCTACCAATAGCTTCTTTTGCCATGGGG"

    # Genetic code dictionary (RNA codon -> Amino Acid single letter code)
    # '_' represents a stop codon.
    codon_table = {
        'AUA': 'I', 'AUC': 'I', 'AUU': 'I', 'AUG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
        'AAC': 'N', 'AAU': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGU': 'S', 'AGA': 'R', 'AGG': 'R',
        'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
        'CAC': 'H', 'CAU': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',
        'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
        'GAC': 'D', 'GAU': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
        'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S',
        'UUC': 'F', 'UUU': 'F', 'UUA': 'L', 'UUG': 'L',
        'UAC': 'Y', 'UAU': 'Y', 'UAA': '_', 'UAG': '_', 'UGA': '_'
    }

    # Step 1: Transcribe DNA to mRNA by replacing T with U
    mrna_sequence = dna_sequence.replace('T', 'U')

    # Step 2: Find the starting position of the first start codon 'AUG'
    start_codon = 'AUG'
    start_index = mrna_sequence.find(start_codon)

    protein_sequence = []
    translated_codons = []

    if start_index != -1:
        # Step 3: Translate the mRNA sequence starting from the start codon
        # Iterate through the sequence in steps of 3 (one codon at a time)
        for i in range(start_index, len(mrna_sequence), 3):
            codon = mrna_sequence[i:i + 3]

            # Ensure we have a full 3-nucleotide codon
            if len(codon) == 3:
                amino_acid = codon_table.get(codon, '?')  # '?' for unknown codons

                # Step 4: Stop translation if a stop codon is found
                if amino_acid == '_':
                    break
                
                protein_sequence.append(amino_acid)
                translated_codons.append(f"{codon} ({amino_acid})")
            else:
                # Reached the end of the sequence without a full codon
                break
    
    # Step 5: Print the final result
    if protein_sequence:
        final_amino_acid_string = "".join(protein_sequence)
        print("The first protein sequence is translated from the following codons:")
        print(" -> ".join(translated_codons))
        print("\nFinal Amino Acid Sequence (single letter code):")
        print(final_amino_acid_string)
        # This is the final answer for the system.
        print(f"\n<<<{final_amino_acid_string}>>>")
    else:
        print("No start codon 'AUG' found in the transcribed mRNA sequence.")
        print("\n<<<No protein translated>>>")

# Execute the function to get the answer
translate_first_protein()