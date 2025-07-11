def identify_products():
    """
    Identifies and describes the two major products from the reaction of
    styrene with tert-butyl peroxybenzoate catalyzed by Fe(OTf)3.
    """

    print("The reaction produces two major regioisomeric products, A and B.")
    print("The order is arbitrary as it is not specified which GCMS peak corresponds to which compound.\n")

    # Product 1
    product_1_name = "1-(benzoyloxy)-2-(tert-butoxy)-1-phenylethane"
    product_1_smiles = "c1ccccc1C(OC(=O)c2ccccc2)COCC(C)(C)C"
    print("Product A:")
    print(f"  Name: {product_1_name}")
    print(f"  SMILES: {product_1_smiles}")
    print("  Description: In this isomer, the benzoyloxy group is attached to the benzylic carbon and the tert-butoxy group is on the terminal carbon.\n")

    # Product 2
    product_2_name = "1-(tert-butoxy)-2-(benzoyloxy)-1-phenylethane"
    product_2_smiles = "c1ccccc1C(OC(C)(C)C)COC(=O)c2ccccc2"
    print("Product B:")
    print(f"  Name: {product_2_name}")
    print(f"  SMILES: {product_2_smiles}")
    print("  Description: In this isomer, the tert-butoxy group is attached to the benzylic carbon and the benzoyloxy group is on the terminal carbon.")

if __name__ == "__main__":
    identify_products()