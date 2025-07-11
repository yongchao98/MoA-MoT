def identify_reagents():
    """
    This function identifies and prints the reagents A and B for the given chemical synthesis.
    The reaction transforms a polycyclic oxygen-containing cation (1) into an N-amino derivative (2),
    and then into a final N-propyl, N-hydro quinacridinium-like structure (3).
    """

    # Step A: 1 -> 2
    # The transformation introduces an N-NH2 group in place of a ring oxygen.
    # This is a characteristic reaction with hydrazine.
    reagent_A = {
        "name": "Hydrazine",
        "formula": "H2N-NH2"
    }

    # Step B: 2 -> 3
    # This step introduces an N-propyl group, for which propylamine is the logical source.
    # The conversion of the existing N-NH2 group to N-H likely occurs under the
    # same reaction conditions.
    reagent_B = {
        "name": "Propylamine",
        "formula": "CH3CH2CH2NH2"
    }

    print("--- Reagent Identification ---")
    print("Reaction A: To convert compound 1 to 2, an oxygen atom is replaced by an N-NH2 group.")
    print(f"Reagent A is: {reagent_A['name']} (Formula: {reagent_A['formula']})")
    print("\nReaction B: To convert compound 2 to 3, another oxygen atom is replaced by an N-propyl group, and the N-NH2 group becomes N-H.")
    print(f"Reagent B is: {reagent_B['name']} (Formula: {reagent_B['formula']})")
    print("----------------------------")

if __name__ == "__main__":
    identify_reagents()