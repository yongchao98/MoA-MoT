import textwrap

def analyze_degeneracy():
    """
    Analyzes RNA sequences from multiple-choice options to determine which best represents
    amino acid degeneracy and its underlying condition.
    """

    # The standard RNA codon table mapping codons to amino acids
    CODON_TABLE = {
        'AUA': 'Isoleucine', 'AUC': 'Isoleucine', 'AUU': 'Isoleucine', 'AUG': 'Methionine',
        'ACA': 'Threonine', 'ACC': 'Threonine', 'ACG': 'Threonine', 'ACU': 'Threonine',
        'AAC': 'Asparagine', 'AAU': 'Asparagine', 'AAA': 'Lysine', 'AAG': 'Lysine',
        'AGC': 'Serine', 'AGU': 'Serine', 'AGA': 'Arginine', 'AGG': 'Arginine',
        'CUA': 'Leucine', 'CUC': 'Leucine', 'CUG': 'Leucine', 'CUU': 'Leucine',
        'CCA': 'Proline', 'CCC': 'Proline', 'CCG': 'Proline', 'CCU': 'Proline',
        'CAC': 'Histidine', 'CAU': 'Histidine', 'CAA': 'Glutamine', 'CAG': 'Glutamine',
        'CGA': 'Arginine', 'CGC': 'Arginine', 'CGG': 'Arginine', 'CGU': 'Arginine',
        'GUA': 'Valine', 'GUC': 'Valine', 'GUG': 'Valine', 'GUU': 'Valine',
        'GCA': 'Alanine', 'GCC': 'Alanine', 'GCG': 'Alanine', 'GCU': 'Alanine',
        'GAC': 'Aspartic Acid', 'GAU': 'Aspartic Acid', 'GAA': 'Glutamic Acid', 'GAG': 'Glutamic Acid',
        'GGA': 'Glycine', 'GGC': 'Glycine', 'GGG': 'Glycine', 'GGU': 'Glycine',
        'UCA': 'Serine', 'UCC': 'Serine', 'UCG': 'Serine', 'UCU': 'Serine',
        'UUC': 'Phenylalanine', 'UUU': 'Phenylalanine', 'UUA': 'Leucine', 'UUG': 'Leucine',
        'UAC': 'Tyrosine', 'UAU': 'Tyrosine', 'UAA': 'Stop', 'UAG': 'Stop', 'UGA': 'Stop',
        'UGC': 'Cysteine', 'UGU': 'Cysteine', 'UGG': 'Tryptophan'
    }

    # The answer choices from the problem
    CHOICES = [
        ('A', "5'-GAUACGUACGAU-3'", 'third position wobble effect.'),
        ('B', "5'-GUUUCAGAUUC-3'", 'presence of inosines.'),
        ('C', "5'-ACGGUCAACGU-3'", 'second position pyrimidine.'),
        ('D', "5'-CUUAUUGAUGU-3'", 'AUA as methionine in mitochondria.'),
        ('E', "5'-AUCGCAGCUAGC-3'", 'use of alternative splicing.')
    ]

    print("Analyzing each sequence for amino acid degeneracy...\n")

    for choice_letter, sequence_raw, condition in CHOICES:
        # Clean up the sequence string to get just the nucleotides
        sequence = sequence_raw.replace("5'-", "").replace("-3'", "")
        
        # Parse the sequence into codons
        codons = textwrap.wrap(sequence, 3)
        
        # Translate codons to amino acids, ignoring incomplete codons
        translatable_codons = [c for c in codons if len(c) == 3]
        amino_acids = [CODON_TABLE.get(c, 'Unknown') for c in translatable_codons]
        
        # Analysis
        unique_amino_acids = sorted(list(set(amino_acids)))
        
        print(f"--- Option {choice_letter} ---")
        print(f"Sequence: {sequence_raw}")
        print(f"Condition: {condition}\n")
        
        print(f"Codons: {translatable_codons}")
        print(f"Resulting Amino Acids: {amino_acids}")
        
        if amino_acids:
            # Outputting the numbers for the "equation" as requested
            print(f"Analysis: This sequence of {len(amino_acids)} codons translates to {len(unique_amino_acids)} unique amino acid(s): {unique_amino_acids}.")
        else:
            print("Analysis: No full codons could be translated.")

        # Check for explicit demonstration of degeneracy
        aa_to_codons = {}
        for i, aa in enumerate(amino_acids):
            if aa not in aa_to_codons:
                aa_to_codons[aa] = set()
            aa_to_codons[aa].add(translatable_codons[i])

        for aa, c_set in aa_to_codons.items():
            if len(c_set) > 1:
                print(f"Degeneracy shown: Codons {sorted(list(c_set))} all code for {aa}.")
        
        print("\n" + "="*40 + "\n")

    print("--- Final Conclusion ---")
    conclusion = (
        "The best choice must have a scientifically correct 'condition' that explains degeneracy. "
        "The conditions for options B, C, D, and E are all incorrect or misapplied in this context. "
        "Option A provides the only correct explanation: the 'third position wobble effect,' which is the fundamental principle of genetic code degeneracy. "
        "Although sequence E (AUC GCA GCU AGC -> Ile-Ala-Ala-Ser) shows a clearer example by using different codons (GCA, GCU) for the same amino acid (Alanine), its associated condition ('alternative splicing') is incorrect. "
        "Therefore, Option A is the superior answer as it correctly pairs a sequence with the valid scientific reason for degeneracy."
    )
    print(textwrap.fill(conclusion, 80))

# Execute the analysis function
analyze_degeneracy()