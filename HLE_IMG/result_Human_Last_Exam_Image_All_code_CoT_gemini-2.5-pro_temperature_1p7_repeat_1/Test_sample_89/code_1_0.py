def get_product_structures():
    """
    This function provides the chemical structures of products A, B, and C
    in SMILES format based on the described reaction pathways.
    """

    # Product A: methyl 5-((S)-1-acetylpyrrolidin-2-yl)-1H-pyrrole-3-carboxylate
    # Formed from the Huisgen cycloaddition pathway.
    product_A_name = "methyl 5-((S)-1-acetylpyrrolidin-2-yl)-1H-pyrrole-3-carboxylate"
    product_A_smiles = "COC(=O)c1cc([C@@H]2CCCN2C(=O)C)cn1"

    # Product B: 1-acetyl-1,5-dihydro-2H-pyrrol-2-one
    # This is the "tethered imide" mentioned, a major product from the fragmentation pathways.
    product_B_name = "1-acetyl-1,5-dihydro-2H-pyrrol-2-one"
    product_B_smiles = "CC(=O)N1C(=O)C=CC1"

    # Product C: 1-((S)-pyrrolidin-2-yl)ethan-1-one
    # This is the "acetyl pyrrolidine" fragment formed during the reaction with acetic anhydride.
    product_C_name = "1-((S)-pyrrolidin-2-yl)ethan-1-one"
    product_C_smiles = "CC(=O)[C@@H]1CCCN1"
    
    print("Structures of the three products (A, B, C) in SMILES format:\n")
    print(f"Product A: {product_A_name}")
    print(f"SMILES: {product_A_smiles}\n")
    
    print(f"Product B: {product_B_name}")
    print(f"SMILES: {product_B_smiles}\n")
    
    print(f"Product C: {product_C_name}")
    print(f"SMILES: {product_C_smiles}")

if __name__ == "__main__":
    get_product_structures()