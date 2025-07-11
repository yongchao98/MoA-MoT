import textwrap

def analyze_genetic_degeneracy():
    """
    Analyzes an RNA sequence to explain amino acid degeneracy based on the
    third-position wobble effect.
    """
    # --- Data Definitions ---
    # The RNA sequence from the correct answer choice
    sequence = "GAUACGUACGAU"

    # A comprehensive RNA codon table mapping codons to single-letter amino acid codes
    codon_table = {
        'AUA': 'I', 'AUC': 'I', 'AUU': 'I', 'AUG': 'M', 'ACA': 'T', 'ACC': 'T',
        'ACG': 'T', 'ACU': 'T', 'AAC': 'N', 'AAU': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGU': 'S', 'AGA': 'R', 'AGG': 'R', 'CUA': 'L', 'CUC': 'L',
        'CUG': 'L', 'CUU': 'L', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
        'CAC': 'H', 'CAU': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGA': 'R', 'CGC': 'R',
        'CGG': 'R', 'CGU': 'R', 'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A', 'GAC': 'D', 'GAU': 'D',
        'GAA': 'E', 'GAG': 'E', 'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
        'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S', 'UUC': 'F', 'UUU': 'F',
        'UUA': 'L', 'UUG': 'L', 'UAC': 'Y', 'UAU': 'Y', 'UGC': 'C', 'UGU': 'C',
        'UGG': 'W', 'UAA': 'Stop', 'UAG': 'Stop', 'UGA': 'Stop'
    }

    # Map for full amino acid names for better readability
    amino_acid_full_names = {
        'D': 'Aspartic Acid (Asp)', 'T': 'Threonine (Thr)', 'Y': 'Tyrosine (Tyr)'
    }

    # Pre-calculate degeneracy for each amino acid
    degeneracy_map = {}
    for codon, aa_code in codon_table.items():
        if aa_code not in degeneracy_map:
            degeneracy_map[aa_code] = 0
        degeneracy_map[aa_code] += 1

    # --- Analysis ---
    print(f"Analyzing Sequence A: 5'-{sequence}-3'")
    print(f"Condition: third position wobble effect.\n")

    print("Step 1: The RNA sequence is parsed into codons.")
    codons = textwrap.wrap(sequence, 3)
    print(f"Codons: {', '.join(codons)}\n")

    print("Step 2: Each codon is translated into an amino acid.")
    polypeptide_codes = [codon_table[c] for c in codons]
    polypeptide_names = [amino_acid_full_names[code] for code in polypeptide_codes]
    print(f"Polypeptide: {' - '.join(polypeptide_names)}\n")

    print("Step 3: The degeneracy for each amino acid is calculated.")
    print("This demonstrates the 'wobble effect', where changes in the third codon position often do not change the specified amino acid.\n")
    
    print("Final Equation (showing degeneracy counts for each amino acid in the sequence):")
    equation_parts = []
    for code in polypeptide_codes:
        degeneracy_count = degeneracy_map[code]
        aa_name = amino_acid_full_names[code]
        equation_parts.append(f"{aa_name} ({degeneracy_count} codons)")

    # Print the "equation" as requested
    final_equation = " + ".join(equation_parts)
    print(final_equation)
    print("\nThe numbers in the equation represent how many distinct codons can code for each amino acid in the resulting polypeptide:")
    print("2 (for Asp) + 4 (for Thr) + 2 (for Tyr) + 2 (for Asp)")

if __name__ == '__main__':
    analyze_genetic_degeneracy()