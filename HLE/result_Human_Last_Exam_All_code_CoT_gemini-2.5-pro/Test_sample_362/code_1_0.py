def solve_wittig_reaction():
    """
    This function analyzes a Wittig reaction and prints the resulting chemical equation.
    """

    # 1. Define the reactants by name and structure
    reactant_aldehyde_name = "pivalaldehyde"
    reactant_aldehyde_structure = "(CH3)3C-CHO"
    reactant_ylide_name = "(2-(2-chlorophenyl)ethylidene)triphenyl-l5-phosphane"
    reactant_ylide_structure = "(2-Cl-C6H4)-CH2-CH=P(Ph)3"

    # 2. Explain the reaction mechanism
    # In a Wittig reaction, the oxygen of the aldehyde is swapped with the
    # carbon group of the ylide to form an alkene and triphenylphosphine oxide.
    # Aldehyde part: (CH3)3C-CH=
    # Ylide part: =CH-CH2-(C6H4Cl)

    # 3. Form the products
    product_alkene_structure = "(CH3)3C-CH=CH-CH2-(C6H4Cl)"
    product_oxide_name = "Triphenylphosphine oxide"
    product_oxide_structure = "O=P(Ph)3"

    # 4. Determine the IUPAC name of the alkene product
    # The longest carbon chain containing the double bond has 5 carbons (pentene).
    # Numbering the chain to give the double bond and substituents the lowest locants:
    # 1-(2-chlorophenyl)-4,4-dimethylpent-2-ene
    product_alkene_name = "1-(2-chlorophenyl)-4,4-dimethylpent-2-ene"

    # 5. Print the full reaction equation
    print("The Wittig reaction is as follows:")
    print(f"{reactant_aldehyde_name} + {reactant_ylide_name} -> {product_alkene_name} + {product_oxide_name}")
    print("\nIn terms of structure:")
    print(f"{reactant_aldehyde_structure} + {reactant_ylide_structure} -> {product_alkene_structure} + {product_oxide_structure}")

    # 6. Output the numbers from the final product's name as requested
    print("\nThe final product in the equation is the alkene:")
    print(product_alkene_name)
    print("\nThe numbers in the IUPAC name of this product are:")
    print("1 (from the '1-(2-chlorophenyl)' part)")
    print("2 (from the '2-chlorophenyl' and 'pent-2-ene' parts)")
    print("4 (from the '4,4-dimethyl' part, appears twice)")
    print("4 (from the '4,4-dimethyl' part, appears twice)")


solve_wittig_reaction()