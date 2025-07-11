import textwrap

def analyze_genetic_degeneracy():
    """
    Analyzes RNA sequences and conditions to determine the best explanation
    for amino acid degeneracy.
    """
    # --- Step 1: Define Genetic Code and Degeneracy Data ---
    genetic_code = {
        'AUA': 'I', 'AUC': 'I', 'AUU': 'I', 'AUG': 'M', 'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
        'AAC': 'N', 'AAU': 'N', 'AAA': 'K', 'AAG': 'K', 'AGC': 'S', 'AGU': 'S', 'AGA': 'R', 'AGG': 'R',
        'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
        'CAC': 'H', 'CAU': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',
        'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
        'GAC': 'D', 'GAU': 'D', 'GAA': 'E', 'GAG': 'E', 'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
        'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S', 'UUC': 'F', 'UUU': 'F', 'UUA': 'L', 'UUG': 'L',
        'UAC': 'Y', 'UAU': 'Y', 'UAA': '_', 'UAG': '_', 'UGC': 'C', 'UGU': 'C', 'UGG': 'W', 'UGA': '_'
    }

    amino_acid_degeneracy = {
        'I': 3, 'M': 1, 'T': 4, 'N': 2, 'K': 2, 'S': 6, 'R': 6, 'L': 6, 'P': 4,
        'H': 2, 'Q': 2, 'V': 4, 'A': 4, 'D': 2, 'E': 2, 'G': 4, 'F': 2, 'Y': 2,
        'C': 2, 'W': 1, '_': 3
    }

    amino_acid_map = {
        'A': 'Alanine', 'R': 'Arginine', 'N': 'Asparagine', 'D': 'Aspartic Acid',
        'C': 'Cysteine', 'E': 'Glutamic Acid', 'Q': 'Glutamine', 'G': 'Glycine',
        'H': 'Histidine', 'I': 'Isoleucine', 'L': 'Leucine', 'K': 'Lysine',
        'M': 'Methionine', 'F': 'Phenylalanine', 'P': 'Proline', 'S': 'Serine',
        'T': 'Threonine', 'W': 'Tryptophan', 'Y': 'Tyrosine', 'V': 'Valine', '_': 'Stop'
    }

    # --- Step 2: Analyze Option A as the primary candidate ---
    sequence_a = "GAUACGUACGAU"
    print(f"Analysis of Option A: Sequence 5'-{sequence_a}-3' and Condition 'third position wobble effect'")
    print("-" * 70)

    # --- Step 3: Break sequence into codons and translate ---
    codons_a = textwrap.wrap(sequence_a, 3)
    amino_acids_a = [genetic_code[c] for c in codons_a]
    peptide_a = "-".join([amino_acid_map[aa] for aa in amino_acids_a])

    print("1. Translation and Degeneracy of Amino Acids:")
    max_degeneracy_in_seq = 0
    max_deg_aa_info = ""

    for i in range(len(codons_a)):
        codon = codons_a[i]
        aa_short = amino_acids_a[i]
        aa_long = amino_acid_map[aa_short]
        degeneracy = amino_acid_degeneracy[aa_short]
        print(f"   Codon: {codon} -> Amino Acid: {aa_long} ({aa_short}). This amino acid is {degeneracy}-fold degenerate.")
        if degeneracy > max_degeneracy_in_seq:
            max_degeneracy_in_seq = degeneracy

    print(f"\nThe resulting peptide is: {peptide_a}")
    print(f"The maximum degeneracy level for an amino acid in this sequence is {max_degeneracy_in_seq}-fold (for Threonine).")
    print("-" * 70)
    
    # --- Step 4: Evaluate the sequence and condition together ---
    print("2. Evaluation:")
    print("The question asks for the best combination of a sequence AND a condition.")
    print("\n*   Sequence Evaluation:")
    print("    - The sequence contains the codon 'ACG' which codes for Threonine.")
    print("    - Threonine is a 4-fold degenerate amino acid, meaning 4 different codons (ACU, ACC, ACA, ACG) specify it.")
    print("    - This demonstrates that the sequence contains elements of high degeneracy.")
    print("\n*   Condition Evaluation:")
    print("    - The condition given is the 'third position wobble effect'.")
    print("    - This is the fundamental, correct scientific principle explaining why multiple codons can code for the same amino acid. The first two bases of the codon form strong pairs with the tRNA anticodon, while the third base has more flexible ('wobble') pairing rules.")
    
    print("\n*   Comparison with other options (e.g., Option E):")
    print("    - Option E's sequence (AUC GCA GCU AGC) is a stronger example because it translates to Ile-Ala-Ala-Ser, explicitly using two different codons (GCA, GCU) for Alanine.")
    print("    - However, Option E's condition ('use of alternative splicing') is completely unrelated to codon degeneracy. Alternative splicing is a process of pre-mRNA modification and does not affect translation.")
    print("-" * 70)

    # --- Step 5: Final Conclusion ---
    print("3. Conclusion:")
    print("Option A is the best answer because it correctly pairs a sequence containing an amino acid with high degeneracy (Threonine) with the precise scientific principle (third position wobble effect) that causes this degeneracy.")
    print("While other sequences might show degeneracy more explicitly, their associated conditions are incorrect, making them invalid choices.")

analyze_genetic_degeneracy()