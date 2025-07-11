def find_reagents():
    """
    This function identifies and prints the reagents A and B for the given reaction scheme.
    """

    # Reagent A: Converts compound 1 to 2 by introducing an N-NH2 group.
    # This is a nucleophilic substitution reaction using hydrazine.
    reagent_A = "Hydrazine (H2N-NH2)"

    # Reagent B: Converts compound 2 to 3, a quinacridinium salt with an N-propyl group.
    # This is a known condensation-cyclization reaction with an aldehyde.
    # The N-propyl group in the product comes from the aldehyde R-CHO where R is a propyl group.
    # The aldehyde is CH3CH2CH2-CHO, which is butanal.
    reagent_B = "Butanal (CH3CH2CH2CHO)"

    print("--- Identified Reagents ---")
    print("Reagent A is:", reagent_A)
    print("Reagent B is:", reagent_B)

# Execute the function to display the results.
find_reagents()