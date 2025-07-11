import collections

def analyze_degeneracy():
    """
    Analyzes RNA sequences for amino acid degeneracy and evaluates the associated conditions.
    """
    # Standard genetic code mapping codons to amino acids
    genetic_code = {
        'AUA':'I', 'AUC':'I', 'AUU':'I', 'AUG':'M', 'ACA':'T', 'ACC':'T',
        'ACG':'T', 'ACU':'T', 'AAC':'N', 'AAU':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGU':'S', 'AGA':'R', 'AGG':'R', 'CUA':'L', 'CUC':'L',
        'CUG':'L', 'CUU':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P',
        'CAC':'H', 'CAU':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R',
        'CGG':'R', 'CGU':'R', 'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A', 'GAC':'D', 'GAU':'D',
        'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G',
        'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S', 'UUC':'F', 'UUU':'F',
        'UUA':'L', 'UUG':'L', 'UAC':'Y', 'UAU':'Y', 'UGC':'C', 'UGU':'C',
        'UGG':'W', 'UAA':'Stop', 'UAG':'Stop', 'UGA':'Stop'
    }

    # Calculate degeneracy for each amino acid
    aa_counts = collections.Counter(genetic_code.values())
    amino_acid_degeneracy = dict(aa_counts)

    # Sequences and conditions from answer choices
    choices = {
        "A": {"seq": "GAUACGUACGAU", "cond": "third position wobble effect."},
        "B": {"seq": "GUUUCAGAUUC", "cond": "presence of inosines."},
        "C": {"seq": "ACGGUCAACGU", "cond": "second position pyrimidine."},
        "D": {"seq": "CUUAUUGAUGU", "cond": "AUA as methionine in mitochondria."},
        "E": {"seq": "AUCGCAGCUAGC", "cond": "use of alternative splicing."}
    }

    print("Analyzing Amino Acid Degeneracy for Each Sequence:\n")

    for key, data in choices.items():
        seq = data["seq"]
        # Break sequence into complete 3-letter codons
        codons = [seq[i:i+3] for i in range(0, len(seq) - len(seq) % 3, 3)]
        amino_acids = [genetic_code.get(c, 'X') for c in codons]
        degeneracies = [amino_acid_degeneracy.get(aa, 0) for aa in amino_acids]

        print(f"Option {key}:")
        print(f"  Sequence: {seq}")
        print(f"  Codons -> Amino Acids: {' -> '.join([f'{c}({aa})' for c, aa in zip(codons, amino_acids)])}")
        print(f"  Degeneracy of each Amino Acid: {degeneracies}")
        print(f"  Condition: {data['cond']}")
        print("-" * 30)

    print("\nConclusion:")
    print("1. Sequences in options B, D, and E contain amino acids with the highest degeneracy (6-fold for Serine or Leucine).")
    print("2. The conditions for C, D, and E are either scientifically incorrect or irrelevant to codon degeneracy.")
    print("   - (C) Degeneracy is primarily in the third position, not the second.")
    print("   - (D) The mitochondrial code for AUA *reduces* its degeneracy compared to the standard code.")
    print("   - (E) Alternative splicing is a pre-translation mechanism affecting which exons are included, not the code itself.")
    print("3. This leaves options A and B. Sequence A shows lower degeneracy than B. Sequence B's condition, the presence of inosine in the tRNA anticodon, is a specific mechanism that enhances the 'wobble' and allows a single tRNA to recognize multiple codons, thus directly supporting degeneracy.")
    print("\nTherefore, Option B provides the best combination of a sequence with high degeneracy and a correct, specific supporting condition.")

    # Final answer output as requested
    final_choice_key = "B"
    final_seq = choices[final_choice_key]["seq"]
    final_codons = [final_seq[i:i+3] for i in range(0, len(final_seq) - len(final_seq) % 3, 3)]
    final_amino_acids = [genetic_code.get(c, 'X') for c in final_codons]
    final_degeneracies = [amino_acid_degeneracy.get(aa, 0) for aa in final_amino_acids]
    
    print("\nFinal Answer Detail (Option B):")
    print(f"The amino acids translated from the codons GUU, UCA, and GAU are Valine, Serine, and Aspartic Acid.")
    print(f"Their respective degeneracy values are {final_degeneracies[0]}, {final_degeneracies[1]}, and {final_degeneracies[2]}.")

analyze_degeneracy()
<<<B>>>