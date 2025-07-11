def identify_reagents():
    """
    This function identifies the reagents A and B from the provided chemical reaction scheme.
    """
    # Reagent A: Converts compound 1 to 2
    # The transformation is the replacement of a ring oxygen (-O-) with an N-amino group (-N-NH2).
    # This is a characteristic reaction of pyrylium salts with hydrazine.
    reagent_A = {
        "name": "Hydrazine",
        "formula": "H2N-NH2"
    }

    # Reagent B: Converts compound 2 to 3
    # The transformation involves replacing another ring oxygen (-O-) with an N-propyl group (-N-CH2CH2CH3).
    # This is a characteristic reaction of pyrylium salts with a primary amine.
    reagent_B = {
        "name": "n-Propylamine",
        "formula": "CH3CH2CH2NH2"
    }

    print("--- Analysis of the Reaction Scheme ---")
    print("\nStep 1: Identifying Reagent A (1 -> 2)")
    print(f"The reaction introduces an N-NH2 group by replacing an oxygen atom. This is achieved using {reagent_A['name']}.")
    print(f"Reagent A is {reagent_A['name']}, with the chemical formula: {reagent_A['formula']}.")

    print("\nStep 2: Identifying Reagent B (2 -> 3)")
    print(f"The reaction introduces an N-propyl group by replacing another oxygen atom. This is achieved using {reagent_B['name']}.")
    print(f"Reagent B is {reagent_B['name']}, with the chemical formula: {reagent_B['formula']}.")

# Execute the function to print the solution
identify_reagents()