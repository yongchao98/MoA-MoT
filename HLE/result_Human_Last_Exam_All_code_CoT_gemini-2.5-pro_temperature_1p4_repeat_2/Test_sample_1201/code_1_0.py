import textwrap
from collections import Counter

def solve_genetic_degeneracy():
    """
    Analyzes RNA sequences to determine which exhibits maximum amino acid degeneracy
    and is paired with a correct supporting condition.
    """

    # 1. Define the standard genetic code (RNA codons to Amino Acid names)
    genetic_code = {
        'UUU': 'Phenylalanine', 'UUC': 'Phenylalanine', 'UUA': 'Leucine', 'UUG': 'Leucine',
        'CUU': 'Leucine', 'CUC': 'Leucine', 'CUA': 'Leucine', 'CUG': 'Leucine',
        'AUU': 'Isoleucine', 'AUC': 'Isoleucine', 'AUA': 'Isoleucine', 'AUG': 'Methionine',
        'GUU': 'Valine', 'GUC': 'Valine', 'GUA': 'Valine', 'GUG': 'Valine',
        'UCU': 'Serine', 'UCC': 'Serine', 'UCA': 'Serine', 'UCG': 'Serine',
        'CCU': 'Proline', 'CCC': 'Proline', 'CCA': 'Proline', 'CCG': 'Proline',
        'ACU': 'Threonine', 'ACC': 'Threonine', 'ACA': 'Threonine', 'ACG': 'Threonine',
        'GCU': 'Alanine', 'GCC': 'Alanine', 'GCA': 'Alanine', 'GCG': 'Alanine',
        'UAU': 'Tyrosine', 'UAC': 'Tyrosine', 'UAA': 'Stop', 'UAG': 'Stop',
        'CAU': 'Histidine', 'CAC': 'Histidine', 'CAA': 'Glutamine', 'CAG': 'Glutamine',
        'AAU': 'Asparagine', 'AAC': 'Asparagine', 'AAA': 'Lysine', 'AAG': 'Lysine',
        'GAU': 'Aspartic Acid', 'GAC': 'Aspartic Acid', 'GAA': 'Glutamic Acid', 'GAG': 'Glutamic Acid',
        'UGU': 'Cysteine', 'UGC': 'Cysteine', 'UGA': 'Stop', 'UGG': 'Tryptophan',
        'CGU': 'Arginine', 'CGC': 'Arginine', 'CGA': 'Arginine', 'CGG': 'Arginine',
        'AGU': 'Serine', 'AGC': 'Serine', 'AGA': 'Arginine', 'AGG': 'Arginine',
        'GGU': 'Glycine', 'GGC': 'Glycine', 'GGA': 'Glycine', 'GGG': 'Glycine'
    }

    # 2. Calculate the degeneracy for each amino acid
    amino_acids = [aa for aa in genetic_code.values() if aa != 'Stop']
    degeneracy_counts = Counter(amino_acids)

    options = {
        'A': {'seq': 'GAUACGUACGAU', 'cond': 'third position wobble effect.'},
        'B': {'seq': 'GUUUCAGAUUC', 'cond': 'presence of inosines.'},
        'C': {'seq': 'ACGGUCAACGU', 'cond': 'second position pyrimidine.'},
        'D': {'seq': 'CUUAUUGAUGU', 'cond': 'AUA as methionine in mitochondria.'},
        'E': {'seq': 'AUCGCAGCUAGC', 'cond': 'use of alternative splicing.'}
    }

    print("Analyzing each option based on sequence degeneracy and the validity of the condition:\n")

    for key, value in options.items():
        seq = value['seq']
        condition = value['cond']
        
        # 3. Analyze each sequence
        codons = [seq[i:i+3] for i in range(0, len(seq) - len(seq)%3, 3)]
        
        # 4. Evaluate degeneracy
        analysis_lines = []
        aa_list = []
        for codon in codons:
            aa = genetic_code.get(codon, "Unknown")
            if aa != "Unknown" and aa != "Stop":
                degeneracy = degeneracy_counts[aa]
                analysis_lines.append(f"{codon} -> {aa} (Degeneracy: {degeneracy})")
                aa_list.append(aa)

        print(f"--- Option {key} ---")
        print(f"Sequence: 5'-{seq}-3'")
        print(f"Condition: {condition}\n")
        print("Analysis:")
        for line in analysis_lines:
            print(f"  {line}")

        if not aa_list:
            print("  Sequence does not code for any standard amino acids.")
        print("-" * 20 + "\n")

    # 5. Synthesize and Conclude
    conclusion = """
    Conclusion:
    The best answer is B. Here's why:

    1.  **Sequence Degeneracy**: The sequence in B (5'-GUUUCAGAUUC-3') codes for Valine-Serine-Aspartic Acid. Serine (coded by UCA) is one of the most degenerate amino acids, with 6 different codons specifying it. Valine (coded by GUU) is also highly degenerate with 4 codons. This sequence therefore results in maximum amino acid degeneracy. Options D and E also feature 6-fold degenerate amino acids, but their conditions are incorrect.

    2.  **Condition Validity**: The condition is "presence of inosines." Inosine is a modified purine base found in the anticodon of certain tRNA molecules. Its presence in the 'wobble' position of the tRNA anticodon allows it to pair with multiple bases (A, U, or C) in the third position of the mRNA codon. This is a specific molecular mechanism (a "unique subset condition") that directly enables and supports the high degeneracy seen for amino acids like Serine, Valine, Alanine, and others.

    In contrast:
    - A: The condition is correct but general. The sequence itself is less degenerate than B's.
    - C: The condition is not the primary cause of degeneracy.
    - D: The condition describes a non-standard genetic code (mitochondrial) and is not a general principle for maximum degeneracy.
    - E: The condition (alternative splicing) is a completely unrelated biological process.

    Therefore, option B correctly pairs a sequence exhibiting maximum degeneracy with a specific, valid molecular mechanism that supports it.
    """
    print(textwrap.dedent(conclusion))

solve_genetic_degeneracy()
<<<B>>>