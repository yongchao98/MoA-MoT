import textwrap

def analyze_degeneracy():
    """
    Analyzes an RNA sequence to demonstrate the principle of genetic code degeneracy,
    specifically the wobble effect, as described in Choice A.
    """
    # Standard genetic code (RNA codons to Amino Acids)
    genetic_code = {
        'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu',
        'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
        'UAU': 'Tyr', 'UAC': 'Tyr', 'UAA': 'Stop', 'UAG': 'Stop',
        'UGU': 'Cys', 'UGC': 'Cys', 'UGA': 'Stop', 'UGG': 'Trp',
        'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
        'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
        'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
        'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
        'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile', 'AUG': 'Met',
        'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
        'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
        'AGU': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
        'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
        'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
        'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
        'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'
    }

    # Create a reverse mapping from amino acids to their codons
    amino_acid_codons = {}
    for codon, aa in genetic_code.items():
        if aa != 'Stop': # We are interested in amino acids
            if aa not in amino_acid_codons:
                amino_acid_codons[aa] = []
            amino_acid_codons[aa].append(codon)

    # Sequence from Choice A
    sequence = "GAUACGUACGAU"

    print(f"Analysis for Choice A\n" + "="*20)
    print(f"Sequence: 5'-{sequence}-3'\n")

    # Step 1: Translate the sequence
    codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
    amino_acids = [genetic_code[c] for c in codons]
    translated_sequence = "-".join(amino_acids)

    print(f"The sequence is broken into codons: {codons}")
    print(f"The translated amino acid sequence is: {translated_sequence}\n")

    # Step 2: Explain the degeneracy and wobble effect for each unique amino acid
    print("Explanation of Degeneracy (Wobble Effect):\n")
    unique_amino_acids = sorted(list(set(amino_acids)))
    for aa in unique_amino_acids:
        all_codons_for_aa = amino_acid_codons[aa]
        first_two_bases = all_codons_for_aa[0][:2]
        explanation = (f"The amino acid {aa} is encoded by the codons: {', '.join(all_codons_for_aa)}.\n"
                       f"The first two bases ('{first_two_bases}') are constant, while the third base can vary. "
                       f"This is the 'third position wobble effect'.")
        print(textwrap.fill(explanation, width=70))
        print("-" * 25)

    print("\nConclusion: Choice A correctly pairs a sequence demonstrating degeneracy with its fundamental cause, the third position wobble effect.")


if __name__ == '__main__':
    analyze_degeneracy()