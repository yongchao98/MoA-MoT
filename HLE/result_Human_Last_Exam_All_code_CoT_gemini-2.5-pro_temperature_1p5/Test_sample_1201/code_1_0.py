import textwrap

def solve_degeneracy_problem():
    """
    Analyzes RNA sequences and conditions to find the one that best explains maximum amino acid degeneracy.
    """
    # Standard RNA codon table
    codon_table = {
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
        'UGG': 'W', 'UAA': 'Stop', 'UAG': 'Stop', 'UGA': 'Stop'
    }

    # Amino acid degeneracy count (number of codons per amino acid)
    degeneracy = {}
    for codon, aa in codon_table.items():
        if aa != 'Stop':
            degeneracy[aa] = degeneracy.get(aa, 0) + 1

    # Full names for clarity in printing
    aa_names = {
        'A': 'Alanine', 'R': 'Arginine', 'N': 'Asparagine', 'D': 'Aspartic Acid',
        'C': 'Cysteine', 'E': 'Glutamic Acid', 'Q': 'Glutamine', 'G': 'Glycine',
        'H': 'Histidine', 'I': 'Isoleucine', 'L': 'Leucine', 'K': 'Lysine',
        'M': 'Methionine', 'F': 'Phenylalanine', 'P': 'Proline', 'S': 'Serine',
        'T': 'Threonine', 'W': 'Tryptophan', 'Y': 'Tyrosine', 'V': 'Valine'
    }

    options = {
        'A': {'seq': 'GAUACGUACGAU', 'cond': 'third position wobble effect.'},
        'B': {'seq': 'GUUUCAGAUUC', 'cond': 'presence of inosines.'},
        'C': {'seq': 'ACGGUCAACGU', 'cond': 'second position pyrimidine.'},
        'D': {'seq': 'CUUAUUGAUGU', 'cond': 'AUA as methionine in mitochondria.'},
        'E': {'seq': 'AUCGCAGCUAGC', 'cond': 'use of alternative splicing.'}
    }

    print("--- Analysis of Answer Choices ---")

    analysis_results = {}
    for key, val in options.items():
        seq = val['seq']
        codons = [seq[i:i+3] for i in range(0, len(seq) - len(seq)%3, 3)]
        translation = [codon_table.get(c, '?') for c in codons]
        amino_acids = [aa for aa in translation if aa != '?']
        
        analysis = f"Option {key}:\n"
        analysis += f"  Sequence: 5'-{val['seq']}-3'\n"
        analysis += f"  Condition: {val['cond']}\n"
        analysis += f"  Codons: {', '.join(codons)}\n"
        analysis += f"  Translation: {', '.join([aa_names.get(aa, 'N/A') for aa in amino_acids])}\n"
        
        # Evaluation
        if key == 'A':
            analysis += "  Evaluation: The 'wobble effect' is the general mechanism for degeneracy, not a specific condition for *maximum* degeneracy. The amino acids produced (Asp, Thr, Tyr) are not the most degenerate."
        elif key == 'B':
            analysis += "  Evaluation: The sequence is not divisible by three. Inosine is found in tRNA anticodons, not mRNA codons, so the condition does not apply to the given sequence."
        elif key == 'C':
            analysis += textwrap.fill("  Evaluation: This is a strong candidate. The condition 'second position pyrimidine' (specifically Cytosine 'C') points to a powerful rule: all 'XCX' codons (Pro, Thr, Ala) result in 4-fold degenerate amino acids. The sequence contains ACG (Threonine) and GUC (Valine), both of which belong to 4-fold degenerate families, providing a good example of the principle.", width=80, subsequent_indent='  ')
        elif key == 'D':
            analysis += "  Evaluation: While the sequence contains 'CUU' for Leucine (6-fold degenerate), the condition refers to a mitochondrial genetic code exception that *reduces* degeneracy for Isoleucine. This is not a principle for *maximum* degeneracy in the standard code."
        elif key == 'E':
            analysis += "  Evaluation: While the sequence contains codons for highly degenerate amino acids like Alanine (4-fold) and Serine (6-fold), the condition 'alternative splicing' is a mechanism of gene regulation that generates protein diversity, and is completely unrelated to codon degeneracy."
        
        print(analysis)
        print("-" * 30)
        
    print("\n--- Conclusion ---")
    conclusion = "Option C provides the best combination of a sequence illustrating high degeneracy and a condition that correctly identifies a specific, systematic cause. The condition 'second position pyrimidine' points to the fact that having a Cytosine (a pyrimidine) in the second position of a codon consistently leads to 4-fold degeneracy (e.g., ACG for Threonine). This is a unique and powerful rule for high degeneracy within the genetic code."
    print(textwrap.fill(conclusion, width=80))

    print("\n--- Final Answer Breakdown ---")
    chosen_seq = options['C']['seq']
    codons = [chosen_seq[i:i+3] for i in range(0, len(chosen_seq) - len(chosen_seq)%3, 3)]
    
    print(f"Sequence: 5'-{options['C']['seq']}-3'")
    print(f"Codons: {', '.join(codons)}")
    
    for codon in codons:
        aa_code = codon_table[codon]
        aa_name = aa_names[aa_code]
        aa_degeneracy = degeneracy[aa_code]
        print(f"Codon {codon} -> {aa_name} (This amino acid is {aa_degeneracy}-fold degenerate)")

solve_degeneracy_problem()
<<<C>>>