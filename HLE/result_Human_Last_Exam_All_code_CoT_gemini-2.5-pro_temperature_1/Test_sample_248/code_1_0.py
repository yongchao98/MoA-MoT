def translate_first_protein(dna_sequence):
    """
    Finds the amino acid sequence of the first protein translated from a given DNA sequence.
    """
    # Codon table mapping RNA codons to single-letter amino acid codes
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
    stop_codons = {'UAA', 'UAG', 'UGA'}
    start_codon = 'AUG'

    print(f"Original 5' to 3' DNA sequence:\n{dna_sequence}\n")

    # Step 1: Transcribe DNA to mRNA
    mrna_sequence = dna_sequence.replace('T', 'U')
    print(f"Transcribed mRNA sequence:\n{mrna_sequence}\n")

    # Step 2: Find the start codon
    start_index = mrna_sequence.find(start_codon)

    protein_sequence = ""
    if start_index == -1:
        print("No start codon (AUG) found. No protein will be translated.")
    else:
        # The part of the mRNA that will be translated
        reading_frame = mrna_sequence[start_index:]
        print(f"Translation starts at the first AUG codon.\nReading frame:\n{reading_frame}\n")

        # Step 3 & 4: Translate codons to amino acids
        print("Translating codons:")
        codons = []
        amino_acids = []

        for i in range(0, len(reading_frame), 3):
            codon = reading_frame[i:i+3]
            # Ensure we have a full codon of 3 bases
            if len(codon) < 3:
                print(f"'{codon}' -> Incomplete codon at the end. Translation stops.")
                break

            amino_acid = codon_table.get(codon, '?')

            # Step 5: Check for stop codon
            if codon in stop_codons:
                print(f"'{codon}' -> STOP codon. Translation terminates.")
                break
            else:
                print(f"'{codon}' -> {amino_acid}")
                protein_sequence += amino_acid

    print("\n-------------------------------------------")
    print("Final Amino Acid Sequence (single-letter code):")
    print(protein_sequence)
    print("-------------------------------------------")

    return protein_sequence

# The provided DNA sequence
dna_seq = "AGTTGCTTCCCTACCAATAGCTTCTTTTGCCATGGGG"
final_protein = translate_first_protein(dna_seq)