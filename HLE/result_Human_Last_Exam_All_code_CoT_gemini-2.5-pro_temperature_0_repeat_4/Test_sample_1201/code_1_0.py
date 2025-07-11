def analyze_rna_sequence():
    """
    Analyzes an RNA sequence to demonstrate translation and degeneracy.
    """
    # Standard Genetic Code mapping codons to amino acids
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

    # Sequence from option A
    rna_sequence = "GAUACGUACGAU"

    # Break the sequence into codons
    codons = [rna_sequence[i:i+3] for i in range(0, len(rna_sequence), 3)]

    # Translate codons into amino acids
    amino_acids = [genetic_code.get(codon, 'Unknown') for codon in codons]

    # Print the analysis
    print(f"Analysis of Option A:")
    print(f"RNA Sequence: 5'-{rna_sequence}-3'")
    
    # Format the codons and amino acids for clear output
    codon_str = " ".join(codons)
    amino_acid_str = "-".join(amino_acids)
    
    print(f"Codons: {codon_str}")
    print(f"Translated Amino Acid Sequence: {amino_acid_str}")
    print("\nExplanation:")
    print("This sequence contains 4 codons that translate into 3 unique amino acids, demonstrating degeneracy.")
    print("The amino acid Aspartic Acid (Asp), coded by GAU, is an example of an amino acid subject to degeneracy.")
    print("The codons for Aspartic Acid are GAU and GAC. This variation in the third position is explained by the 'third position wobble effect'.")
    print("Therefore, this option correctly pairs a sequence with the biological principle that explains its degeneracy.")

analyze_rna_sequence()