def wittig_reaction_product():
    """
    This function outlines the Wittig reaction between pivalaldehyde and a specified ylide,
    and prints the resulting chemical equation.
    """
    # Define the names of the reactants
    aldehyde = "pivalaldehyde"
    ylide = "(2-(2-chlorophenyl)ethylidene)triphenyl-l5-phosphane"

    # Define the names of the products
    main_product = "1-(2-chlorophenyl)-4,4-dimethylpent-2-ene"
    byproduct = "triphenylphosphine oxide"

    # The stoichiometry for this reaction is 1:1 -> 1:1.
    # The numbers in the final product's IUPAC name are essential for defining its structure.
    print("--- Wittig Reaction Equation ---")
    print(f"1 {aldehyde}  +  1 {ylide}  --->  1 {main_product}  +  1 {byproduct}")

    # Provide more detail on the main product
    print("\nDetailed Product Information:")
    print(f"The main organic product is: {main_product}")
    print("Note: Since a new C=C double bond is formed, the product is typically a mixture of (E) and (Z) stereoisomers.")

wittig_reaction_product()