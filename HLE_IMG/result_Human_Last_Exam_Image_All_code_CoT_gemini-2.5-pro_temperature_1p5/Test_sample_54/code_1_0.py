def identify_reagents():
    """
    Identifies and prints the reagents A and B for the given chemical transformations.
    """
    reagent_A = {
        "name": "Hydrazine",
        "formula": "H2N-NH2"
    }

    reagent_B = {
        "name": "Propanal",
        "formula": "CH3CH2CHO"
    }

    print("The problem shows a two-step synthesis starting from compound 1.")
    print("Step 1: Compound 1 -> Compound 2, using Reagent A.")
    print("Step 2: Compound 2 -> Compound 3, using Reagent B.")
    print("\nBased on the analysis of the structural changes:")
    print(f"Reagent A is: {reagent_A['name']} ({reagent_A['formula']})")
    print(f"Reagent B is: {reagent_B['name']} ({reagent_B['formula']})")

identify_reagents()