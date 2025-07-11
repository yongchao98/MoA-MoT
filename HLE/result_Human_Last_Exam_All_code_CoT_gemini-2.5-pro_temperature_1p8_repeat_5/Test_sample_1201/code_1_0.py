import collections

def analyze_degeneracy():
    """
    Analyzes an RNA sequence to calculate its amino acid degeneracy score.
    The function breaks the sequence into codons, identifies the corresponding amino acids,
    and sums the number of codons that code for each of those amino acids.
    """
    # The standard RNA codon table
    genetic_code = {
        'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu',
        'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
        'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile', 'AUG': 'Met',
        'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
        'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
        'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
        'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
        'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
        'UAU': 'Tyr', 'UAC': 'Tyr', 'UAA': 'STOP', 'UAG': 'STOP',
        'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
        'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
        'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
        'UGU': 'Cys', 'UGC': 'Cys', 'UGA': 'STOP', 'UGG': 'Trp',
        'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
        'AGU': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
        'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'
    }

    # Calculate degeneracy counts for each amino acid
    amino_acid_counts = collections.Counter(genetic_code.values())
    
    # Sequence from option A
    rna_sequence = "GAUACGUACGAU"

    print(f"Analyzing RNA sequence: 5'-{rna_sequence}-3'\n")

    codons = [rna_sequence[i:i+3] for i in range(0, len(rna_sequence), 3)]
    
    total_degeneracy = 0
    degeneracy_values = []

    for codon in codons:
        amino_acid = genetic_code[codon]
        degeneracy = amino_acid_counts[amino_acid]
        total_degeneracy += degeneracy
        degeneracy_values.append(str(degeneracy))
        print(f"Codon {codon} -> {amino_acid} (Amino acid is coded by {degeneracy} different codons)")

    # Format the final equation as requested
    equation_str = " + ".join(degeneracy_values)
    
    print("\n-------------------------------------------")
    print("The principle supporting this degeneracy is the 'third position wobble effect'.")
    print("Total degeneracy score for this sequence is the sum of the degeneracy for each amino acid.")
    print(f"Calculation: {equation_str} = {total_degeneracy}")
    print("-------------------------------------------")

# Run the analysis
analyze_degeneracy()