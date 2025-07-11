def solve_wittig_reaction():
    """
    Analyzes a Wittig reaction to determine its product.
    """
    # Step 1: Define the reactants
    aldehyde_name = "pivalaldehyde (2,2-dimethylpropanal)"
    aldehyde_structure = "(CH3)3C-CHO"
    aldehyde_fragment = "(CH3)3C-CH="

    ylide_name = "(2-(2-chlorophenyl)ethylidene)triphenyl-l5-phosphane"
    # Ph3P is triphenylphosphine, -(o-ClC6H4) is an ortho-chlorophenyl group
    ylide_structure = "Ph3P=CH-CH2-(o-ClC6H4)"
    ylide_fragment = "=CH-CH2-(o-ClC6H4)"
    
    # Step 2: Define the byproduct
    byproduct_name = "triphenylphosphine oxide"
    byproduct_structure = "Ph3P=O"

    # Step 3: Combine the fragments to form the product
    product_structure = aldehyde_fragment + ylide_fragment.lstrip('=')

    # Step 4: Analyze the product
    product_iupac_name = "1-(2-chlorophenyl)-4,4-dimethylpent-2-ene"
    product_formula = "C13H17Cl"
    
    # Calculate molar mass
    atomic_weights = {'C': 12.011, 'H': 1.008, 'Cl': 35.453}
    molar_mass = (13 * atomic_weights['C'] + 
                  17 * atomic_weights['H'] + 
                  1 * atomic_weights['Cl'])

    # Step 5: Print the full analysis
    print("Wittig Reaction Analysis")
    print("-" * 45)
    print("Reactant 1 (Aldehyde):")
    print(f"  Name:      {aldehyde_name}")
    print(f"  Structure: {aldehyde_structure}")

    print("\nReactant 2 (Wittig Reagent):")
    print(f"  Name:      {ylide_name}")
    print(f"  Structure: {ylide_structure}")

    print("-" * 45)
    print("The Wittig Reaction Equation:")
    # This line shows each component of the final reaction "equation"
    print(f"{aldehyde_structure} + {ylide_structure} --> {product_structure} + {byproduct_structure}")
    print("-" * 45)

    print("Major Organic Product:")
    print(f"  IUPAC Name: {product_iupac_name}")
    print(f"  Structure:  {product_structure}")
    print(f"  Formula:    {product_formula}")
    print(f"  Molar Mass: {molar_mass:.3f} g/mol")

# Execute the function to print the solution
solve_wittig_reaction()
