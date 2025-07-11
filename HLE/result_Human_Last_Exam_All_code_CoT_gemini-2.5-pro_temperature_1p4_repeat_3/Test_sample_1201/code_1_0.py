import textwrap

def analyze_degeneracy():
    """
    Analyzes RNA sequences to determine which demonstrates maximum degeneracy
    and is supported by the correct unique subset condition.
    """
    codon_table = {
        'AUA': 'Ile', 'AUC': 'Ile', 'AUU': 'Ile',
        'AUG': 'Met',
        'ACA': 'Thr', 'ACC': 'Thr', 'ACG': 'Thr', 'ACU': 'Thr',
        'AAC': 'Asn', 'AAU': 'Asn',
        'AAA': 'Lys', 'AAG': 'Lys',
        'AGC': 'Ser', 'AGU': 'Ser',
        'AGA': 'Arg', 'AGG': 'Arg',
        'GUA': 'Val', 'GUC': 'Val', 'GUG': 'Val', 'GUU': 'Val',
        'GCA': 'Ala', 'GCC': 'Ala', 'GCG': 'Ala', 'GCU': 'Ala',
        'GAC': 'Asp', 'GAU': 'Asp',
        'GAA': 'Glu', 'GAG': 'Glu',
        'GGA': 'Gly', 'GGC': 'Gly', 'GGG': 'Gly', 'GGU': 'Gly',
        'UCA': 'Ser', 'UCC': 'Ser', 'UCG': 'Ser', 'UCU': 'Ser',
        'UUC': 'Phe', 'UUU': 'Phe',
        'UUA': 'Leu', 'UUG': 'Leu',
        'UAC': 'Tyr', 'UAU': 'Tyr',
        'UGC': 'Cys', 'UGU': 'Cys',
        'UGG': 'Trp',
        'CUA': 'Leu', 'CUC': 'Leu', 'CUG': 'Leu', 'CUU': 'Leu',
        'CCA': 'Pro', 'CCC': 'Pro', 'CCG': 'Pro', 'CCU': 'Pro',
        'CAC': 'His', 'CAU': 'His',
        'CAA': 'Gln', 'CAG': 'Gln',
        'CGA': 'Arg', 'CGC': 'Arg', 'CGG': 'Arg', 'CGU': 'Arg',
        'UAA': 'Stop', 'UAG': 'Stop', 'UGA': 'Stop'
    }

    sequences = {
        'A': ("5'-GAUACGUACGAU-3'", "third position wobble effect."),
        'B': ("5'-GUUUCAGAUUC-3'", "presence of inosines."),
        'C': ("5'-ACGGUCAACGU-3'", "second position pyrimidine."),
        'D': ("5'-CUUAUUGAUGU-3'", "AUA as methionine in mitochondria."),
        'E': ("5'-AUCGCAGCUAGC-3'", "use of alternative splicing.")
    }

    print("Analyzing RNA Sequences for Amino Acid Degeneracy...\n")

    for key, (seq_str, condition) in sequences.items():
        # Clean up sequence string to get just the bases
        rna_seq = seq_str.split('-')[1]
        codons = textwrap.wrap(rna_seq, 3)
        
        # We only consider full codons for translation
        full_codons = [c for c in codons if len(c) == 3]
        amino_acids = [codon_table.get(c, 'Unknown') for c in full_codons]
        
        print(f"--- Option {key} ---")
        print(f"Sequence: {seq_str}")
        print(f"Codons: {', '.join(full_codons)}")
        print(f"Amino Acids: {'-'.join(amino_acids)}")
        print(f"Condition: {condition}")
        
        analysis = ""
        # Check for demonstrated degeneracy (multiple codons for the same AA)
        aa_counts = {}
        for i, aa in enumerate(amino_acids):
            if aa not in aa_counts:
                aa_counts[aa] = set()
            aa_counts[aa].add(full_codons[i])

        degeneracy_found = False
        for aa, codon_set in aa_counts.items():
            if len(codon_set) > 1:
                degeneracy_found = True
                analysis += f"Analysis: This sequence demonstrates degeneracy for {aa}, using codons {codon_set}. "
                break
        if not degeneracy_found:
             analysis += "Analysis: This sequence does not use different codons for the same amino acid. "

        # Evaluate the condition
        if key == 'A':
            analysis += "The condition 'third position wobble effect' is the primary mechanism that explains the degeneracy of the genetic code. All amino acids in this sequence (Asp, Thr, Tyr) have codons that exhibit this wobble effect (e.g., Asp is GAU/GAC)."
        elif key == 'E':
            analysis += "The condition 'alternative splicing' is a mechanism for generating protein diversity from a single gene by varying exon inclusion, but it is not the reason for the degeneracy of the genetic code itself."
        else:
             analysis += "The provided condition is either incorrect or not the best explanation."

        print(analysis)
        print("-" * 20 + "\n")

    print("\n--- Conclusion ---")
    print("Option E (Sequence: 5'-AUCGCAGCUAGC-3') is the only one that explicitly demonstrates degeneracy by using two different codons (GCA and GCU) for the same amino acid (Alanine).")
    print("However, the provided condition, 'use of alternative splicing', is mechanistically incorrect as an explanation for code degeneracy.")
    print("\nOption A (Sequence: 5'-GAUACGUACGAU-3') does not demonstrate degeneracy (it repeats the GAU codon), but its amino acids (Asp, Thr, Tyr) are all subject to degeneracy. Crucially, the associated condition, 'third position wobble effect', is the correct and fundamental principle that causes the degeneracy of the genetic code.")
    print("\nIn scientific questions, pairing an example with the correct underlying principle is key. Therefore, Option A is the best answer because it correctly links a sequence composed of codons subject to degeneracy with the precise scientific principle that explains it.")
