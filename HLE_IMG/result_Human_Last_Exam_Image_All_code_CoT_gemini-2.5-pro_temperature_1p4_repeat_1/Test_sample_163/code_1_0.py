def identify_products():
    """
    Identifies and prints the names of the two major products from the reaction
    of styrene with tert-butyl peroxybenzoate catalyzed by Fe(OTf)3.
    """
    product_A_name = "Product A: 2-(tert-butoxy)-1-phenylethyl benzoate"
    product_B_name = "Product B: 2-(benzoyloxy)-1-(tert-butoxy)-1-phenylethane"

    print("The reaction is an iron-catalyzed oxidative difunctionalization of styrene.")
    print("The two major products, A and B, are constitutional isomers.")
    print("-" * 50)
    print(product_A_name)
    print(product_B_name)
    print("-" * 50)
    print("Structure of Product A (SMILES): c1ccccc1C(OC(=O)c2ccccc2)COC(C)(C)C")
    print("Structure of Product B (SMILES): c1ccccc1C(OC(C)(C)C)COC(=O)c2ccccc2")


if __name__ == "__main__":
    identify_products()