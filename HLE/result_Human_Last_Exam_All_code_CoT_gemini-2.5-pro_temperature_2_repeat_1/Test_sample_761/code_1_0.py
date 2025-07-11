def compare_polymers():
    """
    Compares the structural features of homopolynucleotides and polysaccharides.
    """
    print("--- Structural Comparison: Homopolynucleotides vs. Polysaccharides ---\n")

    # --- Polysaccharide Structural Features ---
    print("1. Polysaccharide (e.g., Glycogen, Starch, Cellulose):")
    print("   - Monomer Unit: Monosaccharide (e.g., glucose).")
    print("   - Monomer Structure: A simple sugar, typically a 6-carbon ring (C6H12O6).")
    print("   - Backbone Linkage: Glycosidic bond. This is an ether linkage (C-O-C) formed between two monosaccharides.")
    print("   - Key Atomic Composition: Carbon (C), Hydrogen (H), Oxygen (O).\n")


    # --- Homopolynucleotide Structural Features ---
    print("2. Homopolynucleotide (e.g., poly(A) RNA, single-stranded DNA with only 'T'):")
    print("   - Monomer Unit: Nucleotide.")
    print("   - Monomer Structure: Composed of three distinct parts:")
    print("     a) A phosphate group (PO4).")
    print("     b) A pentose (5-carbon) sugar (Deoxyribose in DNA or Ribose in RNA).")
    print("     c) A nitrogenous base (In a homopolynucleotide, this base is the same, e.g., Adenine).")
    print("   - Backbone Linkage: Phosphodiester bond. This bond connects the 5' carbon of one nucleotide's sugar to the 3' carbon of the next via the phosphate group (Sugar-O-P-O-Sugar).")
    print("   - Key Atomic Composition: Carbon (C), Hydrogen (H), Oxygen (O), Nitrogen (N), and crucially, Phosphorus (P).\n")

    # --- Conclusion ---
    print("--- Conclusion ---")
    print("Polynucleotides and polysaccharides are fundamentally different structures.")
    print("- Their monomers (nucleotide vs. monosaccharide) are completely different.")
    print("- Their backbone linkages (phosphodiester vs. glycosidic) are completely different.")
    print("- Polynucleotides contain Phosphorus (P) and Nitrogen (N) in their core structure, while polysaccharides do not.")
    print("\nTherefore, polynucleotides are not, structurally, polysaccharides.")

if __name__ == '__main__':
    compare_polymers()
