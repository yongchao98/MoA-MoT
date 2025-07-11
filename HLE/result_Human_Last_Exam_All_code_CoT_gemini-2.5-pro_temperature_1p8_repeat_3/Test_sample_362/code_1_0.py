def solve_wittig_reaction():
    """
    This function determines and prints the product of a specified Wittig reaction.
    """
    # 1. Define the reactants based on their names.
    # Aldehyde: pivalaldehyde is (CH3)3C-CHO
    # Ylide: (2-(2-chlorophenyl)ethylidene)triphenyl-l5-phosphane is Ph3P=CH-CH2-(2-Cl-Ph)
    # The ylide is non-stabilized, which favors the Z-alkene product.
    aldehyde_formula = "(CH3)3C-CHO"
    aldehyde_name = "pivalaldehyde"
    
    ylide_formula = "Ph3P=CH-CH2-(C6H4Cl)"
    ylide_name = "(2-(2-chlorophenyl)ethylidene)triphenylphosphorane"

    # 2. Define the products of the reaction.
    # The reaction combines the (CH3)3C-CH group from the aldehyde and the CH-CH2-(C6H4Cl) group from the ylide.
    # Since the ylide is non-stabilized, the Z (cis) isomer is the major product.
    alkene_product_formula = "(Z)-(CH3)3C-CH=CH-CH2-(C6H4Cl)"
    alkene_product_name = "(Z)-1-(2-chlorophenyl)-4,4-dimethylpent-2-ene"
    
    byproduct_formula = "Ph3P=O"
    byproduct_name = "triphenylphosphine oxide"

    # 3. Print the reaction equation and the name of the product.
    print("The Wittig Reaction:")
    # Print the equation with each reactant and product formula
    reaction_equation = (f"  {aldehyde_formula} ({aldehyde_name}) + "
                         f"{ylide_formula} ({ylide_name}) --> "
                         f"{alkene_product_formula} + {byproduct_formula} ({byproduct_name})")
    
    print(reaction_equation)
    print("\n--------------------------------------------------\n")
    print(f"The major organic product is: {alkene_product_name}")

solve_wittig_reaction()