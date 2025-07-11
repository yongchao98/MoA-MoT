def identify_reagents():
    """
    Identifies the reagents A and B from the provided chemical reaction scheme.
    """
    reagent_A = "Hydrazine (H2N-NH2)"
    reagent_B = "n-Propylamine (CH3CH2CH2NH2)"

    print("The reagents for the reaction scheme are identified as follows:\n")
    
    print(f"Reagent A: {reagent_A}")
    print(f"Reagent B: {reagent_B}")
    
    print("\n--- Rationale ---")
    print("Reaction A (1 -> 2):")
    print("The starting material (1) is a trioxatriangulenium cation. It reacts to form compound (2) by replacing an oxygen-containing ring with a nitrogen-containing ring that has an amino (-NH2) group. Hydrazine is the chemical that provides the necessary H2N-N fragment for this transformation.")
    
    print("\nReaction B (2 -> 3):")
    print("Compound (2) is converted to a quinacridinium derivative (3). This step involves adding a propyl group to a nitrogen atom. n-Propylamine is the source of the N-propyl fragment. It reacts to form the second nitrogen heterocycle in the final product's core structure.")

identify_reagents()