def analyze_rna_sequence(rna_sequence):
    """
    Translates an RNA sequence into an amino acid sequence and explains the degeneracy.
    """
    codon_map = {
        'AUA':'Ile', 'AUC':'Ile', 'AUU':'Ile', 'AUG':'Met',
        'ACA':'Thr', 'ACC':'Thr', 'ACG':'Thr', 'ACU':'Thr',
        'AAC':'Asn', 'AAU':'Asn', 'AAA':'Lys', 'AAG':'Lys',
        'AGC':'Ser', 'AGU':'Ser', 'AGA':'Arg', 'AGG':'Arg',
        'CUA':'Leu', 'CUC':'Leu', 'CUG':'Leu', 'CUU':'Leu',
        'CCA':'Pro', 'CCC':'Pro', 'CCG':'Pro', 'CCU':'Pro',
        'CAC':'His', 'CAU':'His', 'CAA':'Gln', 'CAG':'Gln',
        'CGA':'Arg', 'CGC':'Arg', 'CGG':'Arg', 'CGU':'Arg',
        'GUA':'Val', 'GUC':'Val', 'GUG':'Val', 'GUU':'Val',
        'GCA':'Ala', 'GCC':'Ala', 'GCG':'Ala', 'GCU':'Ala',
        'GAC':'Asp', 'GAU':'Asp', 'GAA':'Glu', 'GAG':'Glu',
        'GGA':'Gly', 'GGC':'Gly', 'GGG':'Gly', 'GGU':'Gly',
        'UCA':'Ser', 'UCC':'Ser', 'UCG':'Ser', 'UCU':'Ser',
        'UUC':'Phe', 'UUU':'Phe', 'UUA':'Leu', 'UUG':'Leu',
        'UAC':'Tyr', 'UAU':'Tyr', 'UAA':'Stop', 'UAG':'Stop',
        'UGC':'Cys', 'UGU':'Cys', 'UGA':'Stop', 'UGG':'Trp',
    }

    print(f"Analyzing RNA Sequence: 5'-{rna_sequence}-3'")
    print("-" * 30)

    # Split the sequence into codons
    codons = [rna_sequence[i:i+3] for i in range(0, len(rna_sequence), 3)]
    
    amino_acids = []
    print("Translation Process:")
    for codon in codons:
        amino_acid = codon_map.get(codon, 'Unknown')
        amino_acids.append(amino_acid)
        print(f"Codon {codon} -> Amino Acid {amino_acid}")
        
    final_peptide = "-".join(amino_acids)
    print(f"\nFinal Amino Acid Sequence: {final_peptide}")
    print("-" * 30)

    # Explanation of degeneracy
    print("\nAnalysis of Degeneracy:")
    print("The sequence is GAU-ACG-UAC-GAU, which translates to Asp-Thr-Tyr-Asp.")
    print("The amino acid Aspartic Acid (Asp) appears twice.")
    print("This demonstrates redundancy in the genetic code.")
    print("\nThe associated condition, 'third position wobble effect', is the correct scientific explanation.")
    print("This principle states that the nucleotide in the third position of a codon can often be changed without altering the amino acid.")
    print("For example, for Aspartic Acid:")
    print("Codon GAU -> Aspartic Acid")
    print("Codon GAC -> Aspartic Acid")
    print("The change from U to C in the third position does not change the amino acid, demonstrating the wobble effect.")
    print("\nConclusion: Option A provides a sequence demonstrating redundancy and the correct scientific condition explaining it.")

# The sequence from option A
sequence_A = "GAUACGUACGAU"
analyze_rna_sequence(sequence_A)
<<<A>>>