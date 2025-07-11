def find_reagents():
    """
    This script determines the reagents A and B in the provided chemical reaction scheme.

    The reaction proceeds in two steps:
    1. Compound 1 reacts with Reagent A to form Compound 2.
    2. Compound 2 reacts with Reagent B to form Compound 3.
    """

    # Analysis for Reagent A
    # The first step converts a xanthylium-type cation (1, with a ring O+) to an
    # N-aminoacridinium cation (2, with a ring N+-NH2).
    # This transformation is a standard reaction with hydrazine.
    reagent_A = "Hydrazine"
    formula_A = "H2N-NH2"

    # Analysis for Reagent B
    # The second step is a complex transformation from an acridinium derivative (2) to a
    # quinacridinium salt (3). This involves a skeletal rearrangement and the introduction
    # of an N-propyl group.
    # n-Propylamine is the most logical source for the N-propyl group and can participate
    # in the necessary ring-forming reactions.
    reagent_B = "n-Propylamine"
    formula_B = "CH3CH2CH2NH2"

    print("Based on the chemical transformations shown:")
    print(f"Reagent A is: {reagent_A} (Formula: {formula_A})")
    print(f"Reagent B is: {reagent_B} (Formula: {formula_B})")

find_reagents()