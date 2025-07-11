def find_best_synthesis_technique():
    """
    Analyzes and recommends the best synthesis technique for a long peptide
    containing an unnatural amino acid (UAA).
    """
    # Peptide characteristics
    peptide_length_aa = 100
    unnatural_amino_acid = "azido phenylalanine"

    print(f"Problem: Synthesize a {peptide_length_aa} amino acid peptide with {unnatural_amino_acid}.\n")

    print("--- Evaluating Options ---")

    # Option 1: Direct Solid-Phase Peptide Synthesis (SPPS)
    print("1. Direct Solid-Phase Peptide Synthesis (SPPS):")
    print(f"   - Assessment: Not suitable. For a {peptide_length_aa} aa peptide, the cumulative coupling efficiency would be extremely low, leading to minimal yield of the correct full-length product.\n")

    # Option 2: Native Chemical Ligation (NCL)
    print("2. Native Chemical Ligation (NCL):")
    print("   - Method: Synthesize two smaller fragments (e.g., 50 aa each) via SPPS and then chemically ligate them.")
    print("   - Assessment: A strong and valid approach. It overcomes the length limitation of direct SPPS. The unnatural amino acid can be easily incorporated into one of the fragments.\n")

    # Option 3: Genetic Code Expansion
    print("3. Genetic Code Expansion:")
    print("   - Method: Use a recombinant host (like E. coli) engineered with an orthogonal aminoacyl-tRNA synthetase and tRNA pair. This system incorporates the unnatural amino acid at a designated stop codon during ribosomal protein synthesis.")
    print("   - Assessment: Highly powerful and efficient for long peptides. It leverages the cell's natural machinery for high-yield, high-fidelity synthesis.\n")
    
    print("--- Conclusion ---")
    best_technique = "Genetic Code Expansion"
    explanation = (
        "While NCL is a strong candidate, Genetic Code Expansion is generally the most helpful and scalable technique "
        f"for producing a peptide of {peptide_length_aa} amino acids with a site-specifically incorporated unnatural amino acid. "
        "It avoids the challenges of chemical synthesis for long chains altogether."
    )

    print(f"Most Helpful Technique: {best_technique}")
    print(f"Reasoning: {explanation}")


if __name__ == '__main__':
    find_best_synthesis_technique()