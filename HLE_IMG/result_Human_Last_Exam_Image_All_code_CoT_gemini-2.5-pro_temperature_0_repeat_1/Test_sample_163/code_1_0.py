def print_reaction_summary():
    """
    Prints the summary of the chemical reaction, identifying reactants and products.
    """
    # Reactants
    styrene_name = "Styrene"
    styrene_smiles = "c1ccccc1C=C"
    peroxide_name = "tert-Butyl peroxybenzoate"
    peroxide_smiles = "c1ccccc1C(=O)OOC(C)(C)C"

    # Products (A and B are the two regioisomers)
    product_A_name = "2-(tert-butoxy)-1-phenylethyl benzoate"
    product_A_smiles = "c1ccccc1C(OC(=O)c2ccccc2)COC(C)(C)C"
    
    product_B_name = "2-(benzoyloxy)-1-(tert-butoxy)-1-phenylethane"
    product_B_smiles = "c1ccccc1C(OC(C)(C)C)COC(=O)c2ccccc2"

    print("Reaction Summary:")
    print("=================\n")
    print("Reactants:")
    print(f"- {styrene_name} (SMILES: {styrene_smiles})")
    print(f"- {peroxide_name} (SMILES: {peroxide_smiles})\n")
    
    print("Reaction Equation:")
    print(f"{styrene_name} + {peroxide_name} --[Fe(OTf)3, 80 C]--> Product A + Product B\n")

    print("Products:")
    print("The two major products, A and B, are regioisomers:\n")
    print(f"Product A: {product_A_name}")
    print(f"SMILES: {product_A_smiles}\n")
    
    print(f"Product B: {product_B_name}")
    print(f"SMILES: {product_B_smiles}")

if __name__ == "__main__":
    print_reaction_summary()