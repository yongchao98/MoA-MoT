def identify_products():
    """
    Identifies the two major products (A and B) from the reaction of
    styrene with tert-butyl peroxybenzoate catalyzed by Fe(OTf)3.
    """
    
    # Product A: 2-(tert-butoxy)-1-phenylethyl benzoate
    # Formed from the initial addition of the tert-butoxyl radical to styrene.
    product_A_name = "2-(tert-butoxy)-1-phenylethyl benzoate"
    product_A_smiles = "CC(C)(C)OCC(c1ccccc1)OC(=O)c2ccccc2"

    # Product B: 1-(tert-butoxy)-2-phenylethyl benzoate
    # Formed from the initial addition of the benzoyloxy radical to styrene.
    product_B_name = "1-(tert-butoxy)-2-phenylethyl benzoate"
    product_B_smiles = "CC(C)(C)OC(c1ccccc1)COC(=O)c2ccccc2"

    print("The two major products, A and B, are constitutional isomers formed by the addition of tert-butoxy and benzoyloxy groups across the double bond of styrene.\n")
    
    print("Product A:")
    print(f"  Name: {product_A_name}")
    print(f"  SMILES: {product_A_smiles}")
    print("\nDescription: A phenyl group and a benzoate group are attached to one carbon, and a tert-butoxy group is on the adjacent carbon.\n")

    print("Product B:")
    print(f"  Name: {product_B_name}")
    print(f"  SMILES: {product_B_smiles}")
    print("\nDescription: A phenyl group and a tert-butoxy group are attached to one carbon, and a benzoate group is on the adjacent carbon.")

if __name__ == "__main__":
    identify_products()