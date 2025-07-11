def solve_wittig_reaction():
    """
    This function analyzes the specified Wittig reaction, identifies the products,
    and prints the results in a structured format.
    """
    # 1. Define the reactants and their properties
    aldehyde_name = "pivalaldehyde"
    aldehyde_formula = "C5H10O"
    aldehyde_structure_fragment = "(CH3)3C-CH"

    ylide_name = "(2-(2-chlorophenyl)ethylidene)triphenyl-l5-phosphane"
    ylide_formula = "C26H22ClP"
    ylide_structure_fragment = "CH-CH2-(2-chlorophenyl)"
    
    # 2. Determine the products based on the Wittig reaction mechanism
    alkene_name = "(Z)-1-(2-chlorophenyl)-4,4-dimethylpent-2-ene"
    alkene_formula = "C13H17Cl"
    
    byproduct_name = "triphenylphosphine oxide"
    byproduct_formula = "C18H15OP"

    # 3. Print the analysis and results
    print("--- Wittig Reaction Analysis ---")
    print("\nThe reaction is a Wittig reaction between an aldehyde and a phosphonium ylide.")
    
    print("\nReactants:")
    print(f"  1. Aldehyde: {aldehyde_name} ({aldehyde_formula})")
    print(f"  2. Ylide: {ylide_name} ({ylide_formula})")

    print("\nProducts:")
    print(f"  1. Major Alkene Product: {alkene_name} ({alkene_formula})")
    print("     (The Z-isomer is favored because the reaction uses an unstabilized ylide.)")
    print(f"  2. Byproduct: {byproduct_name} ({byproduct_formula})")

    print("\n--- Final Reaction Equation ---")
    print("The balanced chemical equation for the reaction is:")
    # The following print statement outputs the full equation including each number in the molecular formulas.
    print(f"{aldehyde_formula} + {ylide_formula} -> {alkene_formula} + {byproduct_formula}")

# Execute the function to display the solution
solve_wittig_reaction()
