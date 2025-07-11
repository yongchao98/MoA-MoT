def wittig_reaction_product():
    """
    This script identifies and prints the product of a specific Wittig reaction.
    """
    # Reactants
    aldehyde = "pivalaldehyde [(CH3)3C-CHO]"
    wittig_reagent = "(2-(2-chlorophenyl)ethylidene)triphenyl-l5-phosphane [Ph3P=CH-CH2-(2-Cl-Ph)]"

    # Products
    main_product_name = "1-(2-chlorophenyl)-4,4-dimethylpent-2-ene"
    main_product_structure = "(CH3)3C-CH=CH-CH2-(2-Cl-Ph)"
    by_product = "triphenylphosphine oxide [Ph3P=O]"

    # Print the reaction summary
    print("Wittig Reaction:")
    print(f"Aldehyde: {aldehyde}")
    print(f"Wittig Reagent: {wittig_reagent}")
    print("-" * 20)
    print("The reaction forms a new alkene by replacing the carbonyl oxygen with the ylide's carbon group.")
    print("-" * 20)
    print("Products:")
    print(f"Main Product Name: {main_product_name}")
    print(f"Main Product Structure: {main_product_structure}")
    print(f"By-product: {by_product}")

wittig_reaction_product()