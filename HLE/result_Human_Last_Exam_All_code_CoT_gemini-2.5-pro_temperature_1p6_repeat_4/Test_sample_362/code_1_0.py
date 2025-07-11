def solve_wittig_reaction():
    """
    This script determines the product of a Wittig reaction between
    pivalaldehyde and (2-(2-chlorophenyl)ethylidene)triphenyl-l5-phosphane.
    """

    # 1. Define the reactants
    aldehyde_name = "Pivalaldehyde (2,2-dimethylpropanal)"
    aldehyde_structure = "(CH3)3C-CHO"
    
    ylide_name = "(2-(2-chlorophenyl)ethylidene)triphenyl-l5-phosphane"
    ylide_structure = "(2-Cl-C6H4)-CH2-CH=P(C6H5)3"

    # 2. Describe the reaction
    print("--- Wittig Reaction Analysis ---")
    print(f"Aldehyde: {aldehyde_name} | Structure: {aldehyde_structure}")
    print(f"Ylide: {ylide_name} | Structure: {ylide_structure}\n")

    # 3. Explain the core transformation
    print("The Wittig reaction forms a new carbon-carbon double bond.")
    print("The carbonyl group (C=O) of the aldehyde is replaced by the ylide's carbon group (=CR'R'').")
    print("Aldehyde part providing (CH3)3C-CH= group.")
    print("Ylide part providing =CH-CH2-(C6H4Cl) group.\n")

    # 4. Determine the product structure and stereochemistry
    # The ylide is unstabilized, leading predominantly to the Z-isomer.
    product_stereochemistry = "(Z)"
    
    # 5. Determine the IUPAC name of the product
    # The numbers in the name are 1, 2, 4, 4.
    product_base_name = "1-(2-chlorophenyl)-4,4-dimethylpent-2-ene"
    final_product_name = f"{product_stereochemistry}-{product_base_name}"
    
    byproduct_name = "triphenylphosphine oxide"
    byproduct_structure = "O=P(C6H5)3"

    # 6. Display the final reaction equation
    print("--- Final Reaction Equation ---")
    reaction_string = (
        f"{aldehyde_structure} + {ylide_structure} ---> "
        f"{final_product_name} + {byproduct_structure}"
    )
    print(reaction_string)

solve_wittig_reaction()
