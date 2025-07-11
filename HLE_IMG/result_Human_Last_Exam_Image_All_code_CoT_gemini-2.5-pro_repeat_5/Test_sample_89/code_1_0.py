def get_product_structures():
    """
    This function identifies and prints the structures of the three products A, B, and C
    based on the provided reaction descriptions.
    The structures are presented by their chemical names and SMILES strings.
    """

    # Product A is formed via a cycloaddition, fragmentation, and subsequent reactions.
    # The final product is an acetylated amine.
    product_A_name = "N-(but-3-en-1-yl)acetamide"
    product_A_smiles = "CC(=O)NCCC=C"

    # Product B is identified as the bicyclic ketone from the Michael addition pathway.
    product_B_name = "Tetrahydro-3H-pyrrolizin-3-one"
    # Note: "Tetrahydro" can be ambiguous. This SMILES represents hexahydropyrrolizin-3-one.
    product_B_smiles = "O=C1CN2CCCC2C1"

    # Product C is identified as the acetylated pyrrolidine from the Dakin-West type pathway.
    product_C_name = "N-acetylpyrrolidine"
    product_C_smiles = "CC(=O)N1CCCC1"

    print("The structures of the three products are:\n")

    print(f"Product A: {product_A_name}")
    print(f"SMILES representation for A: {product_A_smiles}\n")

    print(f"Product B: {product_B_name}")
    print(f"SMILES representation for B: {product_B_smiles}\n")

    print(f"Product C: {product_C_name}")
    print(f"SMILES representation for C: {product_C_smiles}")

if __name__ == "__main__":
    get_product_structures()