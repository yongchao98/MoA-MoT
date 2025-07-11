def solve_wittig_reaction():
    """
    Analyzes the Wittig reaction between pivalaldehyde and a given ylide,
    and prints the resulting chemical equation and product details.
    """

    # 1. Define Reactants
    pivalaldehyde = {
        "name": "pivalaldehyde",
        "structure": "(CH3)3C-CHO"
    }

    wittig_reagent = {
        "name": "(2-(2-chlorophenyl)ethylidene)triphenyl-l5-phosphane",
        "structure": "(2-Cl-C6H4)-CH2-CH=P(Ph)3"
    }

    # 2. Determine Product Structures
    # The C=O from the aldehyde and the C=P from the ylide are exchanged.
    alkene_product_structure = "(CH3)3C-CH=CH-CH2-(C6H4Cl-2)"
    side_product_structure = "O=P(Ph)3"

    # 3. Determine the IUPAC name and stereochemistry of the major product
    # The numbers in the name are: 1, 2, 4, 4, 2
    product_name = "(E)-1-(2-chlorophenyl)-4,4-dimethylpent-2-ene"

    # 4. Print the final analysis
    print("### Wittig Reaction Analysis ###")
    print(f"\nReactant 1 (Aldehyde): {pivalaldehyde['name']}")
    print(f"Structure: {pivalaldehyde['structure']}")
    
    print(f"\nReactant 2 (Ylide): {wittig_reagent['name']}")
    print(f"Structure: {wittig_reagent['structure']}")

    print("\n" + "="*40)
    print("Final Reaction Equation:")
    print("="*40)
    
    # Printing the full equation. The stoichiometric coefficients are all 1.
    final_equation = (f"{pivalaldehyde['structure']}  +  {wittig_reagent['structure']}  --->  "
                      f"{alkene_product_structure}  +  {side_product_structure}")
    print(final_equation)

    print("\n" + "="*40)
    print("Major Product Information:")
    print("="*40)
    print(f"Structure: {alkene_product_structure}")
    print(f"IUPAC Name: {product_name}")
    print("\nNote on Stereochemistry:")
    print("The reaction produces the (E)-isomer as the major product. This is because both reactants contain bulky groups (t-butyl and 2-chlorobenzyl).")
    print("Placing these groups on opposite sides of the new double bond minimizes steric repulsion, making the (E) product more stable and favored.")
    
    print("\nThe numbers in the IUPAC name of the final product are: 1, 2, 4, 4, 2.")

# Execute the function
solve_wittig_reaction()