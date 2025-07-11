def identify_reagents():
    """
    Identifies and explains the reagents A and B in the provided reaction scheme.
    """
    reagent_A = {
        "name": "Hydrazine",
        "formula": "H₂N-NH₂"
    }

    reagent_B = {
        "name": "Propylamine (or n-Propylamine)",
        "formula": "CH₃CH₂CH₂-NH₂"
    }

    print("--- Reagent Identification ---")

    # Explanation for Reagent A
    print("\n[Step 1]: Conversion of Compound 1 to Compound 2")
    print(f"Reagent A is: {reagent_A['name']} ({reagent_A['formula']})")
    print("Explanation:")
    print("The transformation from compound 1 to 2 involves the replacement of an oxygen atom in one of the heterocyclic rings with a nitrogen atom that is bonded to an amino group (-NH₂).")
    print("This reaction is a characteristic conversion of a pyrylium-type salt into an N-aminopyridinium-type salt, which is achieved by using hydrazine as the reagent.")

    # Explanation for Reagent B
    print("\n[Step 2]: Conversion of Compound 2 to Compound 3")
    print(f"Reagent B is: {reagent_B['name']} ({reagent_B['formula']})")
    print("Explanation:")
    print("The transformation from compound 2 to 3 involves two changes:")
    print("  1. A second oxygen atom is replaced by a nitrogen atom bearing a propyl group.")
    print("  2. The N-NH₂ group is reduced to an N-H group (deamination).")
    print("Propylamine provides the N-propyl group for the ring substitution. The deamination reaction can occur concurrently under the reaction conditions (e.g., heating).")

    print("\n--- Summary ---")
    print("Reagent A: Hydrazine")
    print("Reagent B: Propylamine")

identify_reagents()