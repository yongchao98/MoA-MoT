def solve_reagents():
    """
    Identifies reagents A and B in the provided chemical reaction scheme.
    """

    # Step 1: Identify Reagent A
    # The reaction from 1 to 2 involves the substitution of a ring oxygen atom
    # with an N-NH2 group. This is a characteristic reaction of pyrylium-like salts
    # with hydrazine.
    reagent_A_name = "Hydrazine"
    reagent_A_formula = "H2N-NH2"

    # Step 2: Identify Reagent B
    # The reaction from 2 to 3 involves two transformations:
    # a) Substitution of another ring oxygen atom with an N-propyl group. This
    #    points to propylamine as the nucleophile.
    # b) Conversion of the N-NH2 group to an N-H group (deamination).
    # Both transformations are achieved in the step labeled "B". The most
    # plausible single reagent that accounts for the addition of the propylamino
    # functionality is propylamine.
    reagent_B_name = "Propylamine"
    reagent_B_formula = "CH3CH2CH2NH2"

    print(f"Reagent A is: {reagent_A_name} ({reagent_A_formula})")
    print(f"Reagent B is: {reagent_B_name} ({reagent_B_formula})")

solve_reagents()