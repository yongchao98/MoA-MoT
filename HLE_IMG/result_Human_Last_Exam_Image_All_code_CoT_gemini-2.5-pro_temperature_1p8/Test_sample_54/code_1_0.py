def solve_reaction():
    """
    Identifies the reagents A and B in the provided chemical reaction scheme.
    """

    # Reagent A converts the pyrylium salt (1) to an N-aminopyridinium salt (2).
    # This reaction involves the replacement of a ring oxygen atom with an N-NH2 group.
    # The reagent for this transformation is hydrazine.
    reagent_A_name = "Hydrazine"
    reagent_A_formula = "H2N-NH2"

    # Reagent B converts compound (2) to the N-propyl quinacridinium salt (3).
    # This is a complex rearrangement and condensation reaction that incorporates a propyl group.
    # The source of the propyl group is n-propylamine.
    reagent_B_name = "n-Propylamine"
    reagent_B_formula = "CH3CH2CH2NH2"

    print("Identification of Reagents:")
    print("-" * 30)
    print(f"Reagent A is: {reagent_A_name}")
    print(f"Chemical Formula of A: {reagent_A_formula}")
    print("-" * 30)
    print(f"Reagent B is: {reagent_B_name}")
    print(f"Chemical Formula of B: {reagent_B_formula}")
    print("-" * 30)

solve_reaction()