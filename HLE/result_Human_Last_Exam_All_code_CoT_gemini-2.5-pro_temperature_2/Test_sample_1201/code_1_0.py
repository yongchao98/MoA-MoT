import collections

def analyze_degeneracy():
    """
    Analyzes RNA sequences to find the one with maximum degeneracy and a correct supporting condition.
    """
    # Standard Genetic Code mapping from codon to amino acid (single letter)
    GENETIC_CODE = {
        'AUA': 'I', 'AUC': 'I', 'AUU': 'I', 'AUG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
        'AAC': 'N', 'AAU': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGU': 'S', 'AGA': 'R', 'AGG': 'R',
        'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
        'CAC': 'H', 'CAU': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',
        'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
        'GAC': 'D', 'GAU': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
        'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S',
        'UUC': 'F', 'UUU': 'F', 'UUA': 'L', 'UUG': 'L',
        'UAC': 'Y', 'UAU': 'Y', 'UGC': 'C', 'UGU': 'C',
        'UGG': 'W', 'UAA': '*', 'UAG': '*', 'UGA': '*'
    }

    # Map single-letter amino acid code to full name
    AMINO_ACID_NAMES = {
        'A': 'Alanine', 'R': 'Arginine', 'N': 'Asparagine', 'D': 'Aspartic Acid',
        'C': 'Cysteine', 'E': 'Glutamic Acid', 'Q': 'Glutamine', 'G': 'Glycine',
        'H': 'Histidine', 'I': 'Isoleucine', 'L': 'Leucine', 'K': 'Lysine',
        'M': 'Methionine', 'F': 'Phenylalanine', 'P': 'Proline', 'S': 'Serine',
        'T': 'Threonine', 'W': 'Tryptophan', 'Y': 'Tyrosine', 'V': 'Valine',
        '*': 'STOP'
    }

    # Calculate degeneracy (number of codons per amino acid)
    amino_acid_degeneracy = collections.Counter(GENETIC_CODE.values())

    # The sequences and conditions from the answer choices
    choices = {
        "A": {"seq": "GAUACGUACGAU", "cond": "third position wobble effect."},
        "B": {"seq": "GUUUCAGAUUC", "cond": "presence of inosines."},
        "C": {"seq": "ACGGUCAACGU", "cond": "second position pyrimidine."},
        "D": {"seq": "CUUAUUGAUGU", "cond": "AUA as methionine in mitochondria."},
        "E": {"seq": "AUCGCAGCUAGC", "cond": "use of alternative splicing."}
    }

    print("Analyzing each option to determine the best sequence-condition pair:\n")

    for key, data in choices.items():
        sequence = data["seq"]
        condition = data["cond"]
        
        print(f"--- Option {key} ---")
        print(f"Sequence: 5'-{sequence}-3'")
        print(f"Condition: {condition}")

        codons = [sequence[i:i+3] for i in range(0, len(sequence), 3) if len(sequence[i:i+3]) == 3]
        
        print("Breakdown:")
        for codon in codons:
            aa_code = GENETIC_CODE.get(codon, "?")
            aa_name = AMINO_ACID_NAMES.get(aa_code, "Unknown")
            degeneracy = amino_acid_degeneracy.get(aa_code, 0)
            # This printout satisfies the "output each number in the final equation" requirement
            print(f"  Codon: {codon} -> Amino Acid: {aa_name} ({aa_code}), Degeneracy: {degeneracy}")

        # Analysis of the condition's validity
        if key == "A":
            print("Evaluation: The 'third position wobble effect' is the correct and fundamental biological mechanism for genetic code degeneracy. The condition is valid.")
        elif key == "E":
            print("Evaluation: 'Alternative splicing' is an unrelated pre-mRNA processing mechanism. The condition is incorrect.")
        else:
            print("Evaluation: The condition is either incorrect (C), irrelevant to the sequence (B), or describes a non-standard exception rather than a general rule (D).")
        print("-" * 20 + "\n")
        
    print("Final Conclusion:")
    print("Option A provides a sequence composed of codons that are all part of degenerate sets, and correctly pairs it with the 'third position wobble effect', which is the primary scientific principle explaining this degeneracy. While sequence E demonstrates degeneracy most clearly, its associated condition ('alternative splicing') is incorrect. Therefore, Option A is the best-supported answer.")


if __name__ == "__main__":
    analyze_degeneracy()