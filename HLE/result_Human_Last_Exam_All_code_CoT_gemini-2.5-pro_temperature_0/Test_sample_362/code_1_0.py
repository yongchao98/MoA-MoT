def solve_wittig_reaction():
    """
    Determines and explains the product of a Wittig reaction.
    """
    # 1. Define the reactants
    aldehyde = "Pivalaldehyde ((CH3)3C-CHO)"
    ylide = "(2-(2-chlorophenyl)ethylidene)triphenyl-l5-phosphane (Ph3P=CH-CH2(C6H4Cl))"

    # 2. Define the products
    product_name = "1-(2-chlorophenyl)-4,4-dimethylpent-2-ene"
    product_structure = "(CH3)3C-CH=CH-CH2-(C6H4Cl)"
    byproduct = "Triphenylphosphine oxide (Ph3P=O)"

    # 3. Print the reaction summary
    print("The Wittig Reaction:")
    print(f"Reactant 1 (Aldehyde): {aldehyde}")
    print(f"Reactant 2 (Ylide): {ylide}")
    print("-" * 30)
    print("The reaction replaces the carbonyl oxygen (=O) of the aldehyde with the carbon group (=CH-R) of the ylide.")
    print("\nOverall Reaction Equation:")
    print(f"(CH3)3C-CHO + Ph3P=CH-CH2(C6H4Cl)  --->  {product_structure} + {byproduct}")
    print("-" * 30)

    # 4. Print the final products
    print(f"Main Product Name: {product_name}")
    print(f"By-product Name: {byproduct}")
    print("\n--- Explanation of the Numbers in the Product Name ---")
    print(f"Product Name: {product_name}")
    print("1: The '(2-chlorophenyl)' group is attached to carbon number 1 of the main chain.")
    print("2: The 'chloro' group is on carbon number 2 of the 'phenyl' ring.")
    print("4,4: Two 'methyl' groups are attached to carbon number 4 of the main chain.")
    print("2: The alkene double bond (C=C) starts at carbon number 2 of the 'pentene' main chain.")

# Execute the function to print the solution
solve_wittig_reaction()