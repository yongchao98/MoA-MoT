def identify_reagents():
    """
    This script identifies the reagents A and B from the provided chemical reaction scheme.
    """

    # Reagent A: Converts compound 1 to 2
    # The reaction shows one of the oxygen atoms in the pyrylium-like core of compound 1
    # being replaced by an N-NH2 group to form compound 2.
    # This is a standard reaction of a pyrylium salt with hydrazine.
    reagent_A_name = "Hydrazine"
    reagent_A_formula = "H2N-NH2"

    # Reagent B: Converts compound 2 to 3
    # The reaction shows another oxygen atom in the core of compound 2 being replaced
    # by an N-propyl group (N-CH2CH2CH3) to form compound 3.
    # This indicates a reaction with a primary amine, specifically propylamine.
    reagent_B_name = "Propylamine"
    reagent_B_formula = "CH3CH2CH2NH2"

    print("The reagents for the reaction scheme are identified as follows:")
    print("-" * 50)
    print("Reagent A:")
    print(f"Name: {reagent_A_name}")
    print(f"Formula: {reagent_A_formula}")
    print("Reasoning: Reagent A introduces the N-NH2 group by reacting with the pyrylium salt moiety in compound 1.")
    print("-" * 50)
    print("Reagent B:")
    print(f"Name: {reagent_B_name}")
    print(f"Formula: {reagent_B_formula}")
    print("Reasoning: Reagent B introduces the N-propyl group by reacting with another pyrylium salt moiety in compound 2.")
    print("-" * 50)

identify_reagents()