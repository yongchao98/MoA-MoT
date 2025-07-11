try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    print("RDKit is not installed. Please install it using: pip install rdkit")
    exit()

def identify_compound_3():
    """
    This function simulates the three-step reaction sequence starting from terpinolene
    to identify the final product, Compound 3.
    """
    # Step 0: Define the starting material, terpinolene
    # SMILES for terpinolene: 1-methyl-4-(propan-2-ylidene)cyclohex-1-ene
    terpinolene_smiles = 'CC1=CCC(C(=C)C)CC1'
    terpinolene = Chem.MolFromSmiles(terpinolene_smiles)
    if terpinolene is None:
        print("Error: Could not parse terpinolene SMILES string.")
        return

    print("--- Simulating the Reaction Sequence ---")
    print(f"Starting Material (Terpinolene): {Chem.MolToSmiles(terpinolene)}\n")

    # Step 1: Epoxidation of the more substituted (exocyclic) double bond
    # The reaction SMARTS targets a carbon in a ring bonded to a non-ring carbon of a double bond.
    # [C;R:1]=[C;!R:2] targets the exocyclic C=C bond.
    rxn1_smarts = '[C;R:1]=[C;!R:2]([C:3])([C:4])>>[C:1]1O[C:2]1([C:3])([C:4])'
    rxn1 = AllChem.ReactionFromSmarts(rxn1_smarts)
    products1 = rxn1.RunReactants((terpinolene,))
    compound1 = products1[0][0]
    Chem.SanitizeMol(compound1)
    
    print("Step 1: Terpinolene + m-CPBA -> Compound 1 (Epoxide)")
    print(f"Equation: {Chem.MolToSmiles(terpinolene)} -> {Chem.MolToSmiles(compound1)}\n")

    # Step 2: Conversion of epoxide to thiirane
    rxn2_smarts = '[C:1]1O[C:2]1>>[C:1]1S[C:2]1'
    rxn2 = AllChem.ReactionFromSmarts(rxn2_smarts)
    products2 = rxn2.RunReactants((compound1,))
    compound2 = products2[0][0]
    Chem.SanitizeMol(compound2)
    
    print("Step 2: Compound 1 + (CH3)2NCHS -> Compound 2 (Thiirane)")
    print(f"Equation: {Chem.MolToSmiles(compound1)} -> {Chem.MolToSmiles(compound2)}\n")

    # Step 3: Reductive desulfurization of the thiirane back to an alkene
    rxn3_smarts = '[C:1]1S[C:2]1>>[C:1]=[C:2]'
    rxn3 = AllChem.ReactionFromSmarts(rxn3_smarts)
    products3 = rxn3.RunReactants((compound2,))
    compound3 = products3[0][0]
    Chem.SanitizeMol(compound3)

    print("Step 3: Compound 2 + LiAlH4 -> Compound 3 (Final Product)")
    print(f"Equation: {Chem.MolToSmiles(compound2)} -> {Chem.MolToSmiles(compound3)}\n")

    # --- Final Identification ---
    # Comparing molecules using InChI is more robust than comparing SMILES strings
    terpinolene_inchi = Chem.MolToInchi(terpinolene)
    compound3_inchi = Chem.MolToInchi(compound3)

    print("--- Conclusion ---")
    if terpinolene_inchi == compound3_inchi:
        final_product_name = "Terpinolene"
        print(f"The final product, Compound 3, is identical to the starting material.")
        print(f"What is compound 3? It is {final_product_name}.")
    else:
        final_product_name = f"Unknown (SMILES: {Chem.MolToSmiles(compound3)})"
        print(f"The final product is not terpinolene. It is: {final_product_name}")
    
    return final_product_name

if __name__ == "__main__":
    identify_compound_3()