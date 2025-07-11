def get_product_structures():
    """
    Prints the deduced structures for products A, B, and C.
    The structures are represented by their names and SMILES strings.
    """

    # Structure for Product A
    print("Product A:")
    print("Description: Formed via Huisgen cycloaddition followed by rearrangement/fragmentation.")
    print("Deduced Name: Methyl 2-(3,4-dihydro-2H-pyrrol-5-yl)pyridine-3-carboxylate")
    print("SMILES: COC(=O)c1cnccc1C1=NCCC1")
    print("-" * 30)

    # Structure for Product B
    print("Product B:")
    print("Description: Formed via Michael addition followed by rearrangement and fragmentation.")
    print("Deduced Name: Hexahydro-3H-pyrrolizin-3-one")
    print("SMILES: O=C1CN2CCCC2C1")
    print("-" * 30)

    # Structure for Product C
    print("Product C:")
    print("Description: Formed via a Dakin-West type reaction with acetic anhydride.")
    print("Deduced Name: 1-Acetyl-2-acetylpyrrolidine")
    print("SMILES: CC(=O)N1CCCC1C(=O)C")
    print("-" * 30)

if __name__ == "__main__":
    get_product_structures()