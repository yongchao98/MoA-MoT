import collections

def analyze_degeneracy():
    """
    Analyzes an RNA sequence to demonstrate the principle of genetic code degeneracy.
    """
    # Standard RNA codon table
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

    # Invert the genetic code to map amino acids to their codons
    aa_to_codons = collections.defaultdict(list)
    for codon, aa in genetic_code.items():
        aa_to_codons[aa].append(codon)

    # Sequence from Option A
    rna_sequence = "GAUACGUACGAU"
    
    # Process the sequence
    codons = [rna_sequence[i:i+3] for i in range(0, len(rna_sequence), 3)]
    amino_acids = [genetic_code[c] for c in codons]
    unique_amino_acids = sorted(list(set(amino_acids)))

    # Print the analysis
    print(f"Analysis of Sequence: 5'-{rna_sequence}-3'")
    print("-" * 30)
    print(f"Codons found: {codons}")
    print(f"Translated amino acid sequence: {'-'.join(amino_acids)}")
    print("\nExamining the degeneracy of the amino acids in the sequence:")

    for aa in unique_amino_acids:
        synonymous_codons = aa_to_codons[aa]
        print(f"- {aa} ({genetic_code[synonymous_codons[0]]}) is coded by {len(synonymous_codons)} codons: {', '.join(synonymous_codons)}")

    print("\nConclusion:")
    print("The primary reason for multiple codons specifying the same amino acid is the 'wobble' in the third base of the codon.")
    print("This means the third nucleotide can often change without altering the amino acid that is coded.")
    print("Option A correctly pairs a sequence with this fundamental principle, making it the best answer.")


if __name__ == '__main__':
    analyze_degeneracy()