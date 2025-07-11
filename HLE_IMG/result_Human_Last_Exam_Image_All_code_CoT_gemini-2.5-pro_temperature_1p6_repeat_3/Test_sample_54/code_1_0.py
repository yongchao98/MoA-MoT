def solve_reaction():
    """
    This function identifies and prints the reagents A and B for the given chemical reaction scheme.
    """
    # Reagent A is responsible for converting compound 1 to 2 by introducing a N-N group.
    reagent_A_name = "Hydrazine"
    reagent_A_formula = "N2H4"

    # Reagent B is responsible for converting compound 2 to 3 by rearranging the
    # nitrogen-containing ring and adding a propyl group.
    reagent_B_name = "Propylamine"
    reagent_B_formula = "CH3CH2CH2NH2"

    print(f"Reagent A is: {reagent_A_name} ({reagent_A_formula})")
    print(f"Reagent B is: {reagent_B_name} ({reagent_B_formula})")

solve_reaction()