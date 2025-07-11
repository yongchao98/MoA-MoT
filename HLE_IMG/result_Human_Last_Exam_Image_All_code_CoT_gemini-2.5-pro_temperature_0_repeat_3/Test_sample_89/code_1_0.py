def get_product_structures():
    """
    This function returns the SMILES strings for products A, B, and C.
    SMILES (Simplified Molecular-Input Line-Entry System) is a standard way
    to represent chemical structures using short ASCII strings.
    """

    # Product A: 2,3-Dihydro-1H-pyrrolo[1,2-a]pyrrole-7-carboxamide
    # Formed via Huisgen cycloaddition, followed by fragmentation and aminolysis.
    # It is a primary amide with a dihydropyrrolizine core.
    product_A_smiles = "NC(=O)c1cn2c(c1)CCC2"

    # Product B: A maleimide ring fused to a pyrrolidine ring.
    # Formed via Michael addition of the intermediate to methyl propiolate,
    # followed by cyclization and fragmentation.
    product_B_smiles = "O=C1C=CC(=O)N2C1CCC2"

    # Product C: A succinimide ring fused to a pyrrolidine ring.
    # Formed via Dakin-West acylation of the intermediate, followed by
    # cyclization and fragmentation. It is the saturated analog of B's core structure.
    product_C_smiles = "O=C1CC(=O)N2C1CCC2"

    print("The deduced structures for products A, B, and C are:")
    print("Represented as SMILES strings:")
    print(f"Product A: {product_A_smiles}")
    print(f"Product B: {product_B_smiles}")
    print(f"Product C: {product_C_smiles}")

if __name__ == "__main__":
    get_product_structures()