def get_product_structures():
    """
    This function defines and prints the structures of products A, B, and C
    in SMILES format based on the provided reaction descriptions.
    """
    # Product A: 2-(Acetylamino)-1-pyrroline
    # Derived from the pyrroline moiety via a described cycloaddition/fragmentation/acetylation pathway.
    product_A_name = "2-(Acetylamino)-1-pyrroline"
    product_A_smiles = "CC(=O)NC1=NCCC1"

    # Product B: N-(1-pyrrolin-2-yl)maleimide
    # Described as a "tethered imide" from Michael addition. The pyrroline is the tether.
    product_B_name = "N-(1-pyrrolin-2-yl)maleimide"
    product_B_smiles = "O=C1C=CC(=O)N1C2=NCCC2"

    # Product C: N-Acetylpyrrolidine
    # Derived from the proline moiety reacting with acetic anhydride.
    product_C_name = "N-Acetylpyrrolidine"
    product_C_smiles = "CC(=O)N1CCCC1"

    print("The structures of the three products are:")
    print(f"Product A: {product_A_name}")
    print(f"SMILES: {product_A_smiles}\n")
    print(f"Product B: {product_B_name}")
    print(f"SMILES: {product_B_smiles}\n")
    print(f"Product C: {product_C_name}")
    print(f"SMILES: {product_C_smiles}")

get_product_structures()