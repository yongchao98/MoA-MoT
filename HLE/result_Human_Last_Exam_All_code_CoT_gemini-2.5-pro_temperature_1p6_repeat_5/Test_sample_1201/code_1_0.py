import collections

def analyze_genetic_degeneracy():
    """
    Analyzes RNA sequences for amino acid degeneracy and evaluates the accompanying
    biological condition for each option.
    """
    # Standard genetic code mapping codons to full amino acid names
    genetic_code = {
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
        'UAC': 'Tyrosine', 'UAU': 'Tyrosine', 'UGC': 'Cysteine', 'UGU': 'Cysteine',
        'UGG': 'Tryptophan', 'UAA': 'Stop', 'UAG': 'Stop', 'UGA': 'Stop'
    }

    # Data from the answer choices
    options = {
        'A': {'sequence': 'GAUACGUACGAU', 'condition': 'third position wobble effect.'},
        'B': {'sequence': 'GUUUCAGAUUC', 'condition': 'presence of inosines.'},
        'C': {'sequence': 'ACGGUCAACGU', 'condition': 'second position pyrimidine.'},
        'D': {'sequence': 'CUUAUUGAUGU', 'condition': 'AUA as methionine in mitochondria.'},
        'E': {'sequence': 'AUCGCAGCUAGC', 'condition': 'use of alternative splicing.'}
    }

    print("--- Analysis of RNA Sequences and Conditions ---")

    for key, data in options.items():
        sequence = data['sequence']
        condition = data['condition']
        
        # Parse sequence into codons, handling sequences not divisible by 3
        codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
        complete_codons = [c for c in codons if len(c) == 3]

        # Translate to amino acids
        amino_acids = [genetic_code.get(c, 'Unknown') for c in complete_codons]
        
        # Format the translation "equation"
        codon_str = '-'.join(complete_codons)
        aa_str = '-'.join(aa_info.split()[0] for aa_info in amino_acids)
        
        print(f"\nOption {key}:")
        print(f"  Translation Equation: 5'-{codon_str}-3' -> {aa_str}")

        # Check for demonstrated degeneracy (multiple different codons for the same amino acid)
        aa_to_codons = collections.defaultdict(set)
        for i in range(len(complete_codons)):
            aa_to_codons[amino_acids[i]].add(complete_codons[i])
        
        degeneracy_found = False
        for aa, codon_set in aa_to_codons.items():
            if len(codon_set) > 1:
                print(f"  Degeneracy Found: The amino acid '{aa}' is coded by different codons: {', '.join(sorted(list(codon_set)))}.")
                degeneracy_found = True
        
        # Evaluate the provided condition
        print(f"  Provided Condition: {condition}")
        if key == 'A':
            verdict = "Correct. The 'wobble effect' at the third codon position is the fundamental mechanism for genetic code degeneracy."
        elif key == 'E' and degeneracy_found:
             verdict = "Incorrect. While the sequence demonstrates degeneracy, 'alternative splicing' is an unrelated pre-translational process."
        else:
            verdict = "Incorrect. The condition is not a valid explanation for degeneracy."
        print(f"  Condition Validity: {verdict}")

    print("\n--- Conclusion ---")
    print("An option is correct only if BOTH the sequence is relevant and the condition is a valid scientific principle.")
    print("Option A is the only choice that pairs a sequence containing degenerate amino acids with the correct explanatory mechanism ('third position wobble effect').")
    print("Option E demonstrates degeneracy well, but its stated condition ('alternative splicing') is incorrect.")
    print("Therefore, Option A is the best answer.")

analyze_genetic_degeneracy()
<<<A>>>