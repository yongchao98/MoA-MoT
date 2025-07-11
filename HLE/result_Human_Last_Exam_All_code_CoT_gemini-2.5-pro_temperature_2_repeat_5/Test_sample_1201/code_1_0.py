def analyze_rna_degeneracy():
    """
    Analyzes RNA sequences for amino acid degeneracy.
    For each sequence, it translates it to an amino acid peptide,
    checks for demonstrated degeneracy, and reports the degeneracy
    potential of the constituent amino acids.
    """
    # Standard genetic code (RNA codon -> Amino Acid)
    genetic_code = {
        'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu',
        'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
        'UAU': 'Tyr', 'UAC': 'Tyr', 'UAA': 'STOP', 'UAG': 'STOP',
        'UGU': 'Cys', 'UGC': 'Cys', 'UGA': 'STOP', 'UGG': 'Trp',
        'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
        'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
        'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
        'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
        'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile',
        'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
        'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
        'AGU': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
        'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
        'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
        'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
        'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'
    }

    # Calculate degeneracy count for each amino acid
    degeneracy_count = {}
    for codon, aa in genetic_code.items():
        if aa != 'STOP':
            degeneracy_count[aa] = degeneracy_count.get(aa, 0) + 1

    options = {
        "A": {"seq": "GAUACGUACGAU", "condition": "third position wobble effect"},
        "B": {"seq": "GUUUCAGAUUC", "condition": "presence of inosines"},
        "C": {"seq": "ACGGUCAACGU", "condition": "second position pyrimidine"},
        "D": {"seq": "CUUAUUGAUGU", "condition": "AUA as methionine in mitochondria"},
        "E": {"seq": "AUCGCAGCUAGC", "condition": "use of alternative splicing"}
    }

    print("--- Analysis of RNA Sequences for Amino Acid Degeneracy ---\n")

    for choice, data in options.items():
        seq = data["seq"]
        condition = data["condition"]
        
        # Parse sequence into codons, handling incomplete ones
        codons = [seq[i:i+3] for i in range(0, len(seq), 3) if i + 3 <= len(seq)]
        
        # Translate to amino acids
        amino_acids = [genetic_code.get(c, "N/A") for c in codons]
        
        # Map amino acids to the unique codons that coded for them in this sequence
        codon_map = {}
        for i, aa in enumerate(amino_acids):
            if aa not in codon_map:
                codon_map[aa] = set()
            codon_map[aa].add(codons[i])

        print(f"Option {choice}:")
        print(f"  Sequence:  5'-{seq}-3'")
        print(f"  Codons:    {', '.join(codons)}")
        print(f"  Peptide:   {'-'.join(amino_acids)}")
        
        # Check if degeneracy is explicitly demonstrated
        is_degenerate = any(len(c_set) > 1 for c_set in codon_map.values())
        print(f"  Demonstrates Degeneracy in Sequence: {'Yes' if is_degenerate else 'No'}")
        
        if is_degenerate:
            for aa, c_set in codon_map.items():
                if len(c_set) > 1:
                    print(f"    -> Multiple codons for {aa}: {', '.join(sorted(list(c_set)))}")

        print(f"  Condition: {condition}\n")

    print("--- Conclusion ---")
    print("Option E shows the highest degree of degeneracy in its sequence (GCA and GCU both code for Alanine).")
    print("However, the provided condition (alternative splicing) is an incorrect explanation for this phenomenon.")
    print("Option A's condition ('third position wobble effect') is the correct fundamental principle explaining all genetic code degeneracy.")
    print("Therefore, Option A provides the best combination of a relevant example and the correct scientific reason.")

analyze_rna_degeneracy()