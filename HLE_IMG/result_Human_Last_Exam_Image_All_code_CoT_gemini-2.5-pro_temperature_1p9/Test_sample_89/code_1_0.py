def get_product_structures():
    """
    This function provides the structures of products A, B, and C based on the reaction description.
    The structures are given as their common names and SMILES strings.
    """
    
    # Product A: N-acetylpyrrolidine
    # Derived from the proline fragment via the Huisgen cycloaddition pathway followed by
    # fragmentation, hydration, decarboxylation, and acetylation.
    product_A_name = "N-acetylpyrrolidine"
    product_A_smiles = "CC(=O)N1CCCC1"

    # Product B: Hexahydropyrrolizin-3-one
    # Derived from the dihydropyrrole and methyl propiolate fragments via the Michael addition pathway.
    # This corresponds to the named "bicyclic tetrahydro-3H-pyrrolizin-3-one".
    product_B_name = "Hexahydropyrrolizin-3-one"
    product_B_smiles = "O=C1CN2C(CCC2)C1"

    # Product C: N-acetyl-2-pyrrolidone
    # An imide-like N-acyl lactam derived from the proline fragment via the Dakin-West-type pathway.
    product_C_name = "N-acetyl-2-pyrrolidone"
    product_C_smiles = "CC(=O)N1CCCC1=O"

    print("Structure of Product A:")
    print(f"Name: {product_A_name}")
    print(f"SMILES: {product_A_smiles}\n")
    
    print("Structure of Product B:")
    print(f"Name: {product_B_name}")
    print(f"SMILES: {product_B_smiles}\n")

    print("Structure of Product C:")
    print(f"Name: {product_C_name}")
    print(f"SMILES: {product_C_smiles}")

# Execute the function to display the results
get_product_structures()