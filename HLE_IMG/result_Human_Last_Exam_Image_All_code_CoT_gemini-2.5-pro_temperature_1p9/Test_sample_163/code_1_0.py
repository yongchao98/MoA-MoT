def identify_reaction_products():
    """
    Identifies and describes the two major products from the Fe-catalyzed
    reaction of styrene with tert-butyl peroxybenzoate.
    """

    # --- Reactant Information ---
    styrene_smiles = "c1ccccc1C=C"
    peroxide_smiles = "CC(C)(C)OOC(=O)c1ccccc1"

    # --- Product Information ---
    # The two products, A and B, are constitutional isomers.
    
    # Product A: formed by adding -OtBu to C2 and -OCOPh to C1 of styrene
    product_A = {
        "name": "2-(tert-butoxy)-1-phenylethyl benzoate",
        "structure": "Ph-CH(OCOPh)-CH2-OtBu",
        "smiles": "c1ccc(cc1)C(C[O]C(C)(C)C)OC(=O)c2ccccc2"
    }

    # Product B: formed by adding -OCOPh to C2 and -OtBu to C1 of styrene
    product_B = {
        "name": "2-tert-butoxy-2-phenylethyl benzoate",
        "structure": "Ph-CH(OtBu)-CH2-OCOPh",
        "smiles": "c1ccc(cc1)C(COC(=O)c2ccccc2)OC(C)(C)C"
    }

    # --- Print Results ---
    print("The reaction yields two major products, A and B, which are constitutional isomers.\n")
    print("The reaction can be represented by the following addition equations:")
    
    # Printing the "equation" with numbers as requested
    print("\nFormation of Product A:")
    print(f"1 Styrene + 1 Peroxide -> 1 Product A")
    print(f"1 ({styrene_smiles}) + 1 ({peroxide_smiles}) -> 1 ({product_A['smiles']})\n")

    print("Formation of Product B:")
    print(f"1 Styrene + 1 Peroxide -> 1 Product B")
    print(f"1 ({styrene_smiles}) + 1 ({peroxide_smiles}) -> 1 ({product_B['smiles']})\n")
    
    print("-" * 30)
    print("Product A Identification:")
    print(f"  Name:      {product_A['name']}")
    print(f"  Structure: {product_A['structure']}")
    
    print("\nProduct B Identification:")
    print(f"  Name:      {product_B['name']}")
    print(f"  Structure: {product_B['structure']}")


# Run the function to display the solution
identify_reaction_products()