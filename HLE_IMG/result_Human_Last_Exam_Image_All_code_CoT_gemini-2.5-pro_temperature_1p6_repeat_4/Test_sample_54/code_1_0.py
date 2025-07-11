def find_reagents():
    """
    Identifies the reagents A and B in the provided chemical reaction scheme.
    """
    
    # Analysis of the reaction from compound 1 to 2
    # The transformation is a replacement of a pyrylium oxygen atom with an N-NH2 group.
    # This is a standard reaction with hydrazine.
    reagent_A_name = "Hydrazine"
    reagent_A_formula = "H2N-NH2"

    # Analysis of the reaction from compound 2 to 3
    # A new ring containing an N-propyl group is formed, and an N-NH2 group is deaminated.
    # The N-propyl group must come from the reagent, pointing to n-propylamine.
    reagent_B_name = "n-Propylamine"
    reagent_B_formula = "CH3CH2CH2NH2"

    print(f"Reagent A is: {reagent_A_name} ({reagent_A_formula})")
    print(f"Reagent B is: {reagent_B_name} ({reagent_B_formula})")

find_reagents()