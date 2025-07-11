def solve_e2_reaction():
    """
    This function analyzes the E2 reaction of (1S,2R)-1-bromo-2-methylcyclohexane
    with potassium tert-butoxide and prints the name of the major product.
    """

    substrate = "(1S,2R)-1-bromo-2-methylcyclohexane"
    base = "Potassium tert-butoxide (a strong, bulky base)"
    reaction_type = "E2 elimination"

    # E2 reactions on cyclohexane rings require a trans-diaxial arrangement
    # of the leaving group (Br) and the proton (H) to be removed.

    # Analysis of the substrate's chair conformations:
    # (1S,2R) means the Br and the methyl group are trans.
    # Conformer 1: Br (equatorial), Methyl (equatorial). This is the major, more stable conformer.
    #              However, the equatorial Br is not anti-periplanar to any beta-protons.
    #              Therefore, E2 reaction cannot occur from this conformer.
    # Conformer 2: Br (axial), Methyl (axial). This is the minor, less stable conformer.
    #              The reaction must proceed through this conformation.

    # Identify beta-protons in the reactive (all-axial) conformer:
    # Proton at C2: Is in an equatorial position because the methyl group is axial. It is NOT trans-diaxial to the Br.
    #               Therefore, elimination to form 1-methylcyclohexene (Zaitsev product) is NOT possible.
    # Proton at C6: Has an axial proton. This axial proton IS trans-diaxial to the axial Br at C1.
    #               Therefore, elimination can ONLY occur by removing the axial proton from C6.

    # The double bond forms between C1 and C6.
    product_name = "3-methylcyclohexene"

    print(f"Substrate: {substrate}")
    print(f"Reagent: {base}")
    print(f"Reaction Type: {reaction_type}")
    print("\nReasoning:")
    print("The E2 reaction requires the leaving group (Br) and a beta-proton to be in a trans-diaxial orientation.")
    print("For (1S,2R)-1-bromo-2-methylcyclohexane, this only occurs in the less stable chair conformer where both Br and the methyl group are axial.")
    print("In this conformer, only the axial proton on C6 is trans-diaxial to the bromine.")
    print("Elimination towards C2 to form the more substituted alkene (1-methylcyclohexene) is stereochemically impossible.")
    print("Therefore, the reaction is forced to form the double bond between C1 and C6.")
    print("\nFinal Product Name:")
    print(product_name)

solve_e2_reaction()