def identify_reagents():
    """
    This script identifies and prints the reagents A and B for the given chemical synthesis.
    """
    # Reagent A: The transformation from 1 to 2 involves replacing a bridging -O- with a bridging -N-NH2.
    # The source for the N-NH2 group is hydrazine.
    reagent_A_name = "Hydrazine"
    reagent_A_formula = "H2N-NH2"

    # Reagent B: The transformation from 2 to 3 involves two changes:
    # 1. Replacing another bridging -O- with a bridging -N-propyl group (-N-CH2CH2CH3).
    # 2. Converting the -N-NH2 group to a -N-H group.
    # The source for the N-propyl group is n-propylamine. The second change likely
    # happens under the reaction conditions.
    reagent_B_name = "n-Propylamine"
    reagent_B_formula = "CH3CH2CH2NH2"

    print("Reagent A is: {}".format(reagent_A_name))
    print("The chemical formula for A is: {}".format(reagent_A_formula))
    print("\n")
    print("Reagent B is: {}".format(reagent_B_name))
    print("The chemical formula for B is: {}".format(reagent_B_formula))

identify_reagents()