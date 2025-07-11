def identify_reaction_products():
    """
    This script identifies the two major products (A and B) from the reaction
    of styrene with tert-butyl peroxybenzoate, catalyzed by Fe(OTf)3.

    The reaction is an oxidative difunctionalization where a tert-butoxy group
    and a benzoyloxy group add across the double bond of styrene, forming
    two constitutional isomers.
    """

    # --- Reactant Information ---
    styrene_name = "Styrene"
    peroxide_name = "tert-butyl peroxybenzoate"

    # --- Product Information ---
    # The two products, A and B, are the following regioisomers.
    # The assignment to 'A' or 'B' is arbitrary.

    # Isomer 1: Formed by addition of tBuO• first, followed by benzoyloxy group.
    product_1_name = "2-(tert-butoxy)-1-phenylethyl benzoate"
    product_1_smiles = "CC(C)(C)OCC(c1ccccc1)OC(=O)c2ccccc2"
    # Structure: Ph-CH(OOCPh)-CH2-OtBu

    # Isomer 2: Formed by addition of PhCOO• first, followed by tert-butoxy group.
    product_2_name = "2-(tert-butoxy)-2-phenylethyl benzoate"
    product_2_smiles = "c1ccccc1C(=O)OCC(c2ccccc2)OC(C)(C)C"
    # Structure: Ph-CH(OtBu)-CH2-OOCPh

    # --- Print the results ---
    print("Analysis of the Chemical Reaction")
    print("="*40)
    print(f"The reaction is: {styrene_name} + {peroxide_name} -> A + B")
    print("\nThe balanced chemical equation with stoichiometric numbers is:")
    print(f"1 {styrene_name} + 1 {peroxide_name} -> 1 A + 1 B")
    print("\n--- Product Identification ---")
    print("The two major products, A and B, are the following constitutional isomers:")

    print("\nProduct 1 (e.g., A):")
    print(f"  IUPAC Name: {product_1_name}")
    print(f"  SMILES Notation: {product_1_smiles}")

    print("\nProduct 2 (e.g., B):")
    print(f"  IUPAC Name: {product_2_name}")
    print(f"  SMILES Notation: {product_2_smiles}")
    print("="*40)


if __name__ == "__main__":
    identify_reaction_products()