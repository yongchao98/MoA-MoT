def identify_products():
    """
    Identifies and describes the three chemical products A, B, and C based on the reaction description.
    Structures are represented by their names and SMILES strings.
    """

    print("The predicted structures of the three products A, B, and C are:\n")

    # --- Product A ---
    product_A_name = "Product A: methyl 1-(4-(acetylamino)butyl)pyrrole-2-carboxylate"
    product_A_smiles = "CC(=O)NCCCCn1cccc1C(=O)OC"
    print(f"{product_A_name}")
    print(f"SMILES: {product_A_smiles}")
    print("Description: A pyrrole ring with a methyl carboxylate group at position 2 and a 4-(acetylamino)butyl group on the nitrogen atom (position 1).\n")

    # --- Product B ---
    product_B_name = "Product B: hexahydropyrrolizin-3-one"
    product_B_smiles = "O=C1CC2N(C1)CCC2"
    print(f"{product_B_name}")
    print(f"SMILES: {product_B_smiles}")
    print("Description: A fused bicyclic system made of two 5-membered rings sharing a nitrogen (a pyrrolizidine skeleton), with a carbonyl group at position 3, forming a lactam.\n")

    # --- Product C ---
    product_C_name = "Product C: N-acetylpyrrolidine"
    product_C_smiles = "CC(=O)N1CCCC1"
    print(f"{product_C_name}")
    print(f"SMILES: {product_C_smiles}")
    print("Description: A five-membered pyrrolidine ring with an acetyl group attached to the nitrogen atom.\n")

if __name__ == "__main__":
    identify_products()