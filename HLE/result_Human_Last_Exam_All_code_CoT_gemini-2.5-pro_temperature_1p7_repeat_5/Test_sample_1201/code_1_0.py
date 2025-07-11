import collections

def analyze_degeneracy():
    """
    Analyzes RNA sequences for amino acid degeneracy and evaluates the correctness
    of the associated biological conditions.
    """

    # Standard genetic code: RNA codon -> Amino Acid (single-letter)
    genetic_code = {
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'AGU': 'S', 'AGC': 'S',
        'UAU': 'Y', 'UAC': 'Y', 'UAA': 'Stop', 'UAG': 'Stop', 'UGA': 'Stop',
        'UGU': 'C', 'UGC': 'C', 'UGG': 'W',
        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I',
        'AUG': 'M',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }

    # Map for amino acid full names
    amino_acid_full_names = {
        'F': 'Phe', 'L': 'Leu', 'I': 'Ile', 'M': 'Met', 'V': 'Val', 'S': 'Ser',
        'P': 'Pro', 'T': 'Thr', 'A': 'Ala', 'Y': 'Tyr', 'H': 'His', 'Q': 'Gln',
        'N': 'Asn', 'K': 'Lys', 'D': 'Asp', 'E': 'Glu', 'C': 'Cys', 'W': 'Trp',
        'R': 'Arg', 'G': 'Gly', 'Stop': 'Stop'
    }

    # Calculate the degeneracy for each amino acid
    amino_acid_degeneracy = collections.Counter(genetic_code.values())

    # The sequences and conditions from the problem
    options = {
        'A': {"seq": "GAUACGUACGAU", "cond": "third position wobble effect."},
        'B': {"seq": "GUUUCAGAUUC", "cond": "presence of inosines."},
        'C': {"seq": "ACGGUCAACGU", "cond": "second position pyrimidine."},
        'D': {"seq": "CUUAUUGAUGU", "cond": "AUA as methionine in mitochondria."},
        'E': {"seq": "AUCGCAGCUAGC", "cond": "use of alternative splicing."}
    }

    print("Analyzing each option based on the standard genetic code:\n")

    for key, data in options.items():
        seq = data['seq']
        condition = data['cond']
        codons = [seq[i:i+3] for i in range(0, len(seq) - len(seq) % 3, 3)]
        
        amino_acids_short = [genetic_code.get(c, '?') for c in codons]
        amino_acids_full = [amino_acid_full_names.get(aa, 'Unknown') for aa in amino_acids_short]
        degeneracy_values = [amino_acid_degeneracy.get(aa, 0) for aa in amino_acids_short]
        total_degeneracy = sum(degeneracy_values)
        
        print(f"--- Option {key} ---")
        print(f"Sequence: {seq}")
        print(f"Condition: {condition}")
        print(f"Codons: {codons}")
        print(f"Translated Amino Acids: {amino_acids_full}")
        print(f"Individual Degeneracy Values: {degeneracy_values}")
        
        # Formatting the final equation string
        equation_str = ' + '.join(map(str, degeneracy_values))
        print(f"Total Degeneracy Score: {equation_str} = {total_degeneracy}")
        print("-" * 20 + "\n")

    print("--- Conclusion ---")
    print("The question asks for the sequence that results in maximum degeneracy AND the condition that supports it.")
    print("1. Numerical Analysis: Sequence E (AUCGCAGCUAGC) yields the highest total degeneracy score of 17 (Ile(3) + Ala(4) + Ala(4) + Ser(6)).")
    print("2. Condition Analysis:")
    print("   - The condition for E ('alternative splicing') is incorrect. Splicing creates different mRNAs but doesn't relate to the genetic code's degeneracy.")
    print("   - The condition for A ('third position wobble effect') is the correct and fundamental scientific principle that explains the degeneracy of the genetic code.")
    print("   - The sequence for A (GAU-ACG-UAC-GAU) is composed entirely of codons from degenerate sets (2-fold, 4-fold, 2-fold, 2-fold), directly illustrating the wobble effect.")
    print("\nBecause the question requires a correct pairing of both sequence and condition, Option A is the best answer. It presents a sequence that is a clear example of the correctly identified fundamental principle of degeneracy, even though another sequence has a slightly higher numerical score when paired with an unrelated concept.")

analyze_degeneracy()
<<<A>>>