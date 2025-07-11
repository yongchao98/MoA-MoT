def get_product_structures():
    """
    This function provides the structures of products A, B, and C
    in the form of their names and SMILES strings.
    """

    # Product A: N-(3,4-dihydro-2H-pyrrol-2-yl)acetamide
    # Formed via cycloaddition, fragmentation to an isocyanate, then hydrolysis,
    # decarboxylation, and acetylation.
    product_A_name = "Product A: N-(3,4-dihydro-2H-pyrrol-2-yl)acetamide"
    product_A_smiles = "CC(=O)NC1=NCCC1"

    # Products B and C: N-acetyl-N-(3,4-dihydro-2H-pyrrol-2-yl)acetamide
    # Both pathways B (Michael Addition) and C (Dakin-West type) are described
    # to yield the same "tethered imide". The most plausible structure for this
    # is the di-acetylated amine of the dihydropyrrole tether.
    product_B_name = "Product B: N-acetyl-N-(3,4-dihydro-2H-pyrrol-2-yl)acetamide"
    product_B_smiles = "CC(=O)N(C1=NCCC1)C(=O)C"

    product_C_name = "Product C: N-acetyl-N-(3,4-dihydro-2H-pyrrol-2-yl)acetamide"
    product_C_smiles = "CC(=O)N(C1=NCCC1)C(=O)C"

    # Print the results
    print("Based on the reaction description, the structures of the three products are:")
    print("-" * 70)
    print(product_A_name)
    print("SMILES representation: ", product_A_smiles)
    print("-" * 70)
    print(product_B_name)
    print("SMILES representation: ", product_B_smiles)
    print("-" * 70)
    print(product_C_name)
    print("SMILES representation: ", product_C_smiles)
    print("-" * 70)
    print("Note: Products B and C are identical.")


if __name__ == "__main__":
    get_product_structures()
