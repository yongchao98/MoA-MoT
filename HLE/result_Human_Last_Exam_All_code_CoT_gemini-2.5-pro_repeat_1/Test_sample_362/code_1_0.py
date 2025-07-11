def solve_wittig_reaction():
    """
    This function identifies and prints the products of the specified Wittig reaction.
    """
    # Define the reactants by name
    aldehyde = "Pivalaldehyde"
    wittig_reagent = "(2-(2-chlorophenyl)ethylidene)triphenyl-l5-phosphane"

    # Determine the products based on the Wittig reaction mechanism
    alkene_product_name = "1-(2-chlorophenyl)-4,4-dimethylpent-2-ene"
    byproduct = "Triphenylphosphine oxide"

    # Print the reaction equation
    print("The reaction is:")
    print(f"{aldehyde} + {wittig_reagent}  --->  {alkene_product_name} + {byproduct}\n")

    # As requested, output the numbers from the IUPAC name of the main organic product
    print("The product of the reaction is the alkene:")
    print(alkene_product_name)
    
    # The numbers in the IUPAC name are the locants for the substituents and the double bond.
    print("\nThe numbers in the final product's name are: 1, 2, 4, 4")

solve_wittig_reaction()