def solve_wittig_reaction():
    """
    This function determines the product of a specific Wittig reaction
    and prints the reaction equation.
    """

    # 1. Define the reactants by their names.
    aldehyde = "Pivalaldehyde (2,2-dimethylpropanal)"
    wittig_reagent = "(2-(2-chlorophenyl)ethylidene)triphenyl-l5-phosphane"

    # 2. Determine the products based on the Wittig reaction mechanism.
    # The reaction couples the aldehyde's (CH3)3C-CH= group (from C=O)
    # with the ylide's =CH-CH2-(2-ClC6H4) group (from C=P).
    # The major alkene product is named based on IUPAC rules, considering
    # the unstabilized nature of the ylide favors the Z-isomer.
    alkene_product = "(Z)-1-(2-chlorophenyl)-4,4-dimethylpent-2-ene"
    byproduct = "Triphenylphosphine oxide"

    # 3. Assemble and print the final chemical equation.
    # The numbers from the IUPAC names are included in the final output.
    print(f"The reaction equation is:")
    print(f"{aldehyde} + {wittig_reagent} ---> {alkene_product} + {byproduct}")
    print("\nThe principal organic product is:")
    print(alkene_product)

solve_wittig_reaction()