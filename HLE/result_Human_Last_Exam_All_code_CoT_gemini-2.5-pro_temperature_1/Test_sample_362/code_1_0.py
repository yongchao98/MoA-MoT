def solve_wittig_reaction():
    """
    This function determines the product of a Wittig reaction and prints the
    details of the chemical equation.
    """
    # Reactant 1: The aldehyde
    aldehyde_name = "Pivalaldehyde"
    aldehyde_smiles = "CC(C)(C)C=O" # (CH3)3C-CHO

    # Reactant 2: The Wittig reagent (phosphorus ylide)
    wittig_reagent_name = "(2-(2-chlorophenyl)ethylidene)triphenyl-l5-phosphane"
    # The reactive part of the ylide is the carbon group double-bonded to phosphorus.
    # Structure: =CH-CH2-(C6H4Cl), where C6H4Cl is a 2-chlorophenyl group.
    ylide_reactive_fragment_smiles = "CCc1ccccc1Cl"

    # --- The Wittig Reaction Mechanism ---
    # The reaction replaces the carbonyl oxygen (=O) from the aldehyde with the
    # carbon group (=CR'R'') from the ylide to form a new alkene (C=C).

    # Step 1: Get the aldehyde's carbon skeleton by conceptually removing the oxygen.
    aldehyde_carbon_fragment = aldehyde_smiles.replace("=O", "=")

    # Step 2: The ylide's carbon fragment is what will be added.
    ylide_carbon_fragment = ylide_reactive_fragment_smiles

    # Step 3: Combine the two fragments to form the final alkene product.
    alkene_product_smiles = aldehyde_carbon_fragment + ylide_carbon_fragment
    alkene_product_name = "1-(2-chlorophenyl)-4,4-dimethylpent-2-ene"

    # Step 4: The other reaction product is triphenylphosphine oxide.
    oxide_product_name = "Triphenylphosphine oxide"
    oxide_product_smiles = "O=P(c1ccccc1)(c2ccccc2)c3ccccc3"

    # --- Output the results ---
    print("--- Wittig Reaction Analysis ---")
    print(f"\nReactants:")
    print(f"1. Aldehyde: {aldehyde_name} ({aldehyde_smiles})")
    print(f"2. Wittig Reagent: {wittig_reagent_name}")

    print(f"\nProducts:")
    print(f"1. Alkene: {alkene_product_name} ({alkene_product_smiles})")
    print(f"2. Oxide: {oxide_product_name} ({oxide_product_smiles})")

    print("\n--- Final Chemical Equation ---")
    print(f"{aldehyde_name} + {wittig_reagent_name}  --->  {alkene_product_name} + {oxide_product_name}")

    print("\nAs requested, here are the numbers from the main product's IUPAC name:")
    print("Product: 1-(2-chlorophenyl)-4,4-dimethylpent-2-ene")
    print("Numbers in the name: 1, 2, 4, 4")
    print("These numbers specify the positions of the '2-chlorophenyl' group (on carbon 1), the double bond (starting at carbon 2), and the two 'methyl' groups (both on carbon 4).")

solve_wittig_reaction()