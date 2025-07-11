def solve_wittig_reaction():
    """
    This script determines and displays the product of a Wittig reaction
    between pivalaldehyde and a specific phosphorus ylide.
    """

    # 1. Define the reactants by name
    aldehyde = "Pivalaldehyde ((CH3)3C-CHO)"
    wittig_reagent = "(2-(2-chlorophenyl)ethylidene)triphenyl-l5-phosphane"
    wittig_reagent_structure = "((2-Cl-Ph)-CH2-CH=PPh3)"

    # 2. Define the products
    # The reaction replaces the C=O and C=PPh3 bonds with a C=C bond.
    # The main organic product's IUPAC name is determined from its structure.
    main_product_name = "1-(2-chlorophenyl)-4,4-dimethylpent-2-ene"
    main_product_structure = "((CH3)3C-CH=CH-CH2-(2-Cl-Ph))"
    byproduct = "Triphenylphosphine oxide (O=PPh3)"

    # 3. Print the full reaction equation
    print("The Wittig Reaction:")
    print(f"Reactant 1 (Aldehyde): {aldehyde}")
    print(f"Reactant 2 (Ylide): {wittig_reagent} with structure {wittig_reagent_structure}")
    print("\n" + "="*50 + "\n")
    print("The reaction combines the carbon backbone of the aldehyde with that of the ylide, forming an alkene.")
    print("\nFinal Equation:")
    print(f"{aldehyde}  +  {wittig_reagent}")
    print("      |")
    print("      V")
    print(f"{main_product_name}  +  {byproduct}")
    print("\nStructural representation of the main product:")
    print(main_product_structure)


solve_wittig_reaction()
