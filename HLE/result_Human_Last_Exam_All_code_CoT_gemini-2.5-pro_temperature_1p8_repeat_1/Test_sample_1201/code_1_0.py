import collections

def analyze_degeneracy():
    """
    Analyzes the amino acid degeneracy of a given RNA sequence and explains the reasoning.
    """
    # Standard genetic code mapping codons to amino acids (full names)
    genetic_code = {
        'UUU': 'Phenylalanine', 'UUC': 'Phenylalanine', 'UUA': 'Leucine', 'UUG': 'Leucine',
        'UCU': 'Serine', 'UCC': 'Serine', 'UCA': 'Serine', 'UCG': 'Serine',
        'UAU': 'Tyrosine', 'UAC': 'Tyrosine', 'UAA': 'Stop', 'UAG': 'Stop',
        'UGU': 'Cysteine', 'UGC': 'Cysteine', 'UGA': 'Stop', 'UGG': 'Tryptophan',
        'CUU': 'Leucine', 'CUC': 'Leucine', 'CUA': 'Leucine', 'CUG': 'Leucine',
        'CCU': 'Proline', 'CCC': 'Proline', 'CCA': 'Proline', 'CCG': 'Proline',
        'CAU': 'Histidine', 'CAC': 'Histidine', 'CAA': 'Glutamine', 'CAG': 'Glutamine',
        'CGU': 'Arginine', 'CGC': 'Arginine', 'CGA': 'Arginine', 'CGG': 'Arginine',
        'AUU': 'Isoleucine', 'AUC': 'Isoleucine', 'AUA': 'Isoleucine', 'AUG': 'Methionine',
        'ACU': 'Threonine', 'ACC': 'Threonine', 'ACA': 'Threonine', 'ACG': 'Threonine',
        'AAU': 'Asparagine', 'AAC': 'Asparagine', 'AAA': 'Lysine', 'AAG': 'Lysine',
        'AGU': 'Serine', 'AGC': 'Serine', 'AGA': 'Arginine', 'AGG': 'Arginine',
        'GUU': 'Valine', 'GUC': 'Valine', 'GUA': 'Valine', 'GUG': 'Valine',
        'GCU': 'Alanine', 'GCC': 'Alanine', 'GCA': 'Alanine', 'GCG': 'Alanine',
        'GAU': 'Aspartic Acid', 'GAC': 'Aspartic Acid', 'GAA': 'Glutamic Acid', 'GAG': 'Glutamic Acid',
        'GGU': 'Glycine', 'GGC': 'Glycine', 'GGA': 'Glycine', 'GGG': 'Glycine',
    }

    # Calculate the degeneracy for each amino acid
    amino_acid_counts = collections.Counter(genetic_code.values())
    
    # Sequence from option A
    sequence = 'GAUACGUACGAU'
    
    print("Analyzing Option A: Sequence: 5'-GAUACGUACGAU-3'; Condition: third position wobble effect.\n")
    print("The 'wobble effect' correctly identifies the primary mechanism for genetic code degeneracy, where the third base of a codon can often be changed without altering the resulting amino acid. This sequence serves as a clear example of this principle.\n")

    # Split sequence into codons
    codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
    amino_acids = [genetic_code[c] for c in codons]
    degeneracy_values = [amino_acid_counts[aa] for aa in amino_acids]

    # Format for output
    amino_acid_names_str = " - ".join([f"{aa}({c})" for aa, c in zip(amino_acids, codons)])
    degeneracy_equation_str = " + ".join(map(str, degeneracy_values))
    
    print(f"Sequence codons: {'-'.join(codons)}")
    print(f"Translated amino acids: {amino_acid_names_str}")
    print("\nThe degeneracy of an amino acid is the number of codons that code for it:")
    for aa, deg in sorted(list(set(zip(amino_acids, degeneracy_values)))):
        print(f"- {aa} is {deg}-fold degenerate.")
        
    print("\nFinal Degeneracy Equation:")
    print(f"Deg({amino_acids[0]}) + Deg({amino_acids[1]}) + Deg({amino_acids[2]}) + Deg({amino_acids[3]}) = {degeneracy_equation_str}")


analyze_degeneracy()