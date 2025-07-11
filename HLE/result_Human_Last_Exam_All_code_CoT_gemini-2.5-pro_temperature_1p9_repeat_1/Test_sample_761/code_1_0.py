def explain_polymer_structures():
    """
    This script compares the structures of polysaccharides and polynucleotides
    to answer whether they are structurally similar.
    """

    # --- Step 1: Define the building blocks (monomers) ---
    # A polysaccharide is made of simple sugar units (monosaccharides).
    monosaccharide = "(Sugar)"
    
    # A polynucleotide is made of nucleotide units.
    # Each nucleotide has three parts: a phosphate, a sugar, and a base.
    nucleotide = "(Phosphate-Sugar-Base)"

    # --- Step 2: Define the bonds linking the monomers ---
    polysaccharide_bond = "-[Glycosidic Bond]-"
    polynucleotide_bond = "-[Phosphodiester Bond]-"

    # --- Step 3: "Build" a small polymer of each type to visualize ---
    # A chain of 3 units for clarity.
    
    # Equation for a polysaccharide (e.g., starch, cellulose)
    polysaccharide_chain = (
        monosaccharide + polysaccharide_bond +
        monosaccharide + polysaccharide_bond +
        monosaccharide
    )

    # Equation for a homopolynucleotide (e.g., Poly-A tail of mRNA)
    # Even in a "homo-" polynucleotide, the repeating unit is still the complex nucleotide.
    polynucleotide_chain = (
        nucleotide + polynucleotide_bond +
        nucleotide + polynucleotide_bond +
        nucleotide
    )

    # --- Step 4: Print the comparison and conclusion ---
    print("Let's compare the structure of a Polysaccharide and a Polynucleotide.\n")
    
    print("Structural Equation of a Polysaccharide:")
    print(polysaccharide_chain)
    print(f"Its repeating monomer is a {monosaccharide}.\n")

    print("Structural Equation of a Polynucleotide:")
    print(polynucleotide_chain)
    print(f"Its repeating monomer is a {nucleotide}.\n")

    print("--- Conclusion ---")
    print("As the structures show, these two polymers are fundamentally different:")
    print("1. Monomers: The repeating unit of a polysaccharide is a simple sugar. The repeating unit of a polynucleotide is a much more complex nucleotide (sugar + phosphate + base).")
    print("2. Backbone: The backbone of a polysaccharide is sugar-sugar-sugar. The backbone of a polynucleotide is a sugar-phosphate-sugar-phosphate chain.")
    print("3. Bonds: The monomers are joined by different chemical bonds (Glycosidic vs. Phosphodiester).")
    print("\nTherefore, a polynucleotide (or a homopolynucleotide) is NOT, structurally, a polysaccharide.")

# Execute the function to get the answer.
explain_polymer_structures()