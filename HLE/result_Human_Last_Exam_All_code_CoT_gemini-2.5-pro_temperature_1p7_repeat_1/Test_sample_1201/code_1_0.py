def analyze_genetic_code_degeneracy():
    """
    Analyzes an RNA sequence to demonstrate genetic code degeneracy and explains the
    underlying principle.
    """
    # The standard genetic code mapping codons to amino acids.
    genetic_code = {
        'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu',
        'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
        'UAU': 'Tyr', 'UAC': 'Tyr', 'UAA': 'Stop', 'UAG': 'Stop',
        'UGU': 'Cys', 'UGC': 'Cys', 'UGA': 'Stop', 'UGG': 'Trp',
        'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
        'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
        'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
        'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
        'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile',
        'AUG': 'Met',
        'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
        'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
        'AGU': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
        'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
        'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
        'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
        'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly',
    }

    # Sequence from Answer Choice A
    rna_sequence = "GAUACGUACGAU"
    print(f"Analyzing Sequence: 5'-{rna_sequence}-3'\n")

    # Step 1: Break the sequence into codons
    codons = [rna_sequence[i:i+3] for i in range(0, len(rna_sequence), 3)]
    print(f"Codons: {codons}")

    # Step 2: Translate codons to amino acids
    amino_acids = [genetic_code[codon] for codon in codons]
    print(f"Translated Amino Acid Sequence: {amino_acids}\n")

    # Step 3 & 4: Analyze degeneracy and explain the condition
    print("Analysis:")
    print("The sequence contains four codons: GAU, ACG, UAC, and GAU.")
    print("When translated, this produces the amino acid sequence: Aspartic Acid (Asp), Threonine (Thr), Tyrosine (Tyr), Aspartic Acid (Asp).")
    print("The amino acid Aspartic Acid (Asp) is encoded by the codon GAU at two different positions in the polypeptide chain.")
    print("\nExplanation of the Condition:")
    print("The condition provided is the 'third position wobble effect'. This is the correct and fundamental reason for the degeneracy of the genetic code.")
    print("The wobble effect describes how the third nucleotide of a codon can often be changed without affecting the amino acid being coded for.")
    print("For example, Aspartic Acid is coded by both GAU and GAC.")
    print("This 'wobble' allows some amino acids like Leucine, Serine, and Arginine to be coded by up to 6 different codons, which is the 'maximum amino acid degeneracy'.")
    print("Other answer choices are incorrect because their sequences are of invalid length (B, C, D) or their stated conditions are scientifically unrelated to codon degeneracy (B, D, E).")
    print("Therefore, sequence A paired with the 'third position wobble effect' is the best answer.")


if __name__ == '__main__':
    analyze_genetic_code_degeneracy()
<<<A>>>