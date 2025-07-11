import textwrap

def analyze_degeneracy():
    """
    Analyzes RNA sequences for amino acid degeneracy and evaluates the correctness
    of the associated biological conditions.
    """
    
    # Step 1: Define the standard genetic code (RNA codons to Amino Acid)
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
        'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'
    }

    # Step 2: Calculate degeneracy for each amino acid
    amino_acid_counts = {}
    for codon, aa in genetic_code.items():
        if aa not in amino_acid_counts:
            amino_acid_counts[aa] = 0
        amino_acid_counts[aa] += 1
    
    degeneracy_map = amino_acid_counts

    # Define the answer choices
    sequences = {
        'A': {
            "seq": "GAUACGUACGAU",
            "cond": "third position wobble effect."
        },
        'B': {
            "seq": "GUUUCAGAUUC",
            "cond": "presence of inosines."
        },
        'C': {
            "seq": "ACGGUCAACGU",
            "cond": "second position pyrimidine."
        },
        'D': {
            "seq": "CUUAUUGAUGU",
            "cond": "AUA as methionine in mitochondria."
        },
        'E': {
            "seq": "AUCGCAGCUAGC",
            "cond": "use of alternative splicing."
        }
    }

    print("Analyzing each sequence for amino acid degeneracy...\n")
    
    # Step 3 & 4: Analyze each sequence
    for option, data in sequences.items():
        seq = data["seq"]
        condition = data["cond"]
        
        codons = textwrap.wrap(seq, 3)
        # Ensure we only consider full codons
        if len(seq) % 3 != 0:
            codons = codons[:-1]

        translated_aas = [genetic_code.get(c, 'Unknown') for c in codons]
        degeneracy_values = [degeneracy_map.get(aa, 0) for aa in translated_aas]

        print(f"--- Option {option} ---")
        print(f"Sequence: 5'-{seq}-3'")
        print(f"Condition: {condition}")
        print(f"Codons: {', '.join(codons)}")
        print(f"Amino Acids: {', '.join(translated_aas)}")
        
        # Output the "equation" for total degeneracy
        equation_str = " + ".join(map(str, degeneracy_values))
        total_degeneracy = sum(degeneracy_values)
        print(f"Degeneracy per Amino Acid: {equation_str} = {total_degeneracy}\n")

    # Step 5: Final conclusion based on the analysis
    print("--- Conclusion ---")
    explanation = """
While sequences in options B, D, and E contain amino acids with the maximum possible 6-fold degeneracy (Serine or Leucine), their associated conditions are scientifically incorrect or irrelevant to the mechanism of degeneracy for the given RNA sequence.
- (B) Inosines are found in tRNA anticodons, not mRNA codons.
- (D) The condition refers to a non-standard mitochondrial genetic code, which is not the default.
- (E) Alternative splicing is a process that occurs before translation and does not affect codon degeneracy.

Option A presents a sequence containing 'ACG', which codes for Threonine. Threonine has a 4-fold degeneracy (ACU, ACC, ACA, ACG). The provided condition, 'third position wobble effect', is the precise and correct scientific principle that explains this 4-fold degeneracy, where the first two bases determine the amino acid, and the third base can vary. Therefore, Option A provides the best and most accurate pairing of a sequence and its corresponding degeneracy principle.
"""
    print(textwrap.dedent(explanation))

if __name__ == '__main__':
    analyze_degeneracy()