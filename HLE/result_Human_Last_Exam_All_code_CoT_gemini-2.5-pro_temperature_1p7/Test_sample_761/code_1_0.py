def compare_biopolymers():
    """
    Prints a detailed comparison between polynucleotides and polysaccharides.
    """
    print("Are polynucleotides, structurally, polysaccharides? The answer is No.\n")
    print("Here is a detailed structural comparison:\n")

    # --- Polynucleotide Description ---
    print("--- Polynucleotides (e.g., DNA, RNA) ---")
    print("1. Monomer (Building Block):")
    print("   - Nucleotide")
    print("   - A nucleotide consists of three parts:")
    print("     a) A Phosphate Group (contains phosphorus)")
    print("     b) A 5-Carbon Sugar (Deoxyribose or Ribose)")
    print("     c) A Nitrogenous Base (contains nitrogen)")
    print("\n2. Backbone Structure:")
    print("   - A repeating chain of Sugar-Phosphate-Sugar-Phosphate...")
    print("\n3. Linking Bond:")
    print("   - Phosphodiester Bond (links a sugar to a phosphate group)")
    print("\n4. Elemental Composition:")
    print("   - Contains Carbon, Hydrogen, Oxygen, Phosphorus, and Nitrogen (CHOPN)")
    print("-" * 40)

    # --- Polysaccharide Description ---
    print("\n--- Polysaccharides (e.g., Starch, Cellulose) ---")
    print("1. Monomer (Building Block):")
    print("   - Monosaccharide (a simple sugar, e.g., glucose)")
    print("\n2. Backbone Structure:")
    print("   - A repeating chain of Sugar-Sugar-Sugar...")
    print("\n3. Linking Bond:")
    print("   - Glycosidic Bond (links a sugar to another sugar)")
    print("\n4. Elemental Composition:")
    print("   - Typically contains only Carbon, Hydrogen, and Oxygen (CHO)")
    print("-" * 40)

    # --- Conclusion ---
    print("\nConclusion:")
    print("While both are polymers and both contain sugar units, the presence of phosphate groups and nitrogenous bases in the fundamental structure of polynucleotides makes them a completely different class of molecule from polysaccharides.")

if __name__ == "__main__":
    compare_biopolymers()
