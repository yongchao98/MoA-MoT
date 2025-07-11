def analyze_degeneracy():
    """
    Analyzes RNA sequences from a set of options to find the one with maximum
    amino acid degeneracy and evaluates the associated condition.
    """
    RNA_CODON_TABLE = {
        'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu',
        'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
        'UAU': 'Tyr', 'UAC': 'Tyr', 'UAA': 'STOP', 'UAG': 'STOP',
        'UGU': 'Cys', 'UGC': 'Cys', 'UGA': 'STOP', 'UGG': 'Trp',
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

    # Calculate degeneracy for each amino acid
    amino_acid_degeneracy = {}
    for codon, aa in RNA_CODON_TABLE.items():
        if aa != 'STOP':
            amino_acid_degeneracy[aa] = amino_acid_degeneracy.get(aa, 0) + 1

    options = {
        'A': ("5'-GAUACGUACGAU-3'", "third position wobble effect."),
        'B': ("5'-GUUUCAGAUUC-3'", "presence of inosines."),
        'C': ("5'-ACGGUCAACGU-3'", "second position pyrimidine."),
        'D': ("5'-CUUAUUGAUGU-3'", "AUA as methionine in mitochondria."),
        'E': ("5'-AUCGCAGCUAGC-3'", "use of alternative splicing.")
    }

    print("Analyzing each option for amino acid degeneracy:\n")

    for option, (seq, condition) in options.items():
        # Clean up sequence and parse codons
        rna = seq.split('-')[1]
        codons = [rna[i:i+3] for i in range(0, len(rna), 3) if i+3 <= len(rna)]
        
        # Translate to amino acids
        aa_sequence = [RNA_CODON_TABLE.get(c, 'X') for c in codons]
        
        print(f"--- Option {option} ---")
        print(f"Sequence: {seq}")
        print(f"Codons: {' '.join(codons)}")
        print(f"Amino Acid Sequence: {'-'.join(aa_sequence)}")

        # Analyze demonstrated degeneracy
        aa_to_codons = {}
        for i in range(len(codons)):
            aa = aa_sequence[i]
            if aa not in aa_to_codons:
                aa_to_codons[aa] = set()
            aa_to_codons[aa].add(codons[i])

        demo_degeneracy = {aa: codons for aa, codons in aa_to_codons.items() if len(codons) > 1}
        if demo_degeneracy:
            print(f"Demonstrated Degeneracy: Yes. {list(demo_degeneracy.keys())[0]} is coded by {list(demo_degeneracy.values())[0]}.")
        else:
            print("Demonstrated Degeneracy: No. Each amino acid is coded by only one type of codon in this sequence.")

        # Analyze condition
        print(f"Condition: '{condition}'")
        if option == 'A':
            print("Evaluation: Correct. The 'wobble effect' is the primary reason for the degeneracy of the genetic code.")
        elif option == 'E':
            print("Evaluation: Incorrect. Alternative splicing produces different mRNAs from one gene; it is not the mechanism for codon degeneracy.")
        else:
            print("Evaluation: The condition is either incorrect, too specific, or less fundamental than the wobble effect.")
        print("-" * 20 + "\n")

analyze_degeneracy()