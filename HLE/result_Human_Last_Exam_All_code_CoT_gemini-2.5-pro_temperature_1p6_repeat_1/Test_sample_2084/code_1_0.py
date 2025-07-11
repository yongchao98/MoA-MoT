import sys
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
except ImportError:
    print("RDKit library not found. Please install it using 'pip install rdkit-pypi'")
    sys.exit(1)

def solve_reaction():
    """
    Solves the chemical problem by deducing reactant and product structures,
    calculating their molar masses, and identifying the product with the higher molar mass.
    """
    # Step 1 & 3: Define the deduced reactant molecule based on the problem's hints.
    # Reactant: 2-methyl-7-phenyl-1,4,6-trioxaspiro[4.4]nonane
    # This structure is consistent with a spiro-orthoester containing methyl and phenyl groups.
    reactant_smiles = "c1ccc(C2CCOC3(O2)OC(C)CO3)cc1"
    reactant_mol = Chem.MolFromSmiles(reactant_smiles)
    if reactant_mol is None:
        print(f"Error: Could not parse reactant SMILES: {reactant_smiles}")
        return

    # Reactant H2O
    water_smiles = "O"
    water_mol = Chem.MolFromSmiles(water_smiles)

    # Step 4: Determine the hydrolysis product.
    # The first hydrolysis product is a hydroxy-ester, which has the highest molar mass.
    # It is formed from propane-1,2-diol and 4-hydroxy-4-phenylbutanoic acid.
    # Structure: 2-hydroxypropyl 4-hydroxy-4-phenylbutanoate
    product_smiles = "CC(O)COC(=O)CCC(O)c1ccccc1"
    product_mol = Chem.MolFromSmiles(product_smiles)
    if product_mol is None:
        print(f"Error: Could not parse product SMILES: {product_smiles}")
        return
        
    # Step 5: Calculate molar masses to verify the product is heavier.
    reactant_mass = Descriptors.MolWt(reactant_mol)
    water_mass = Descriptors.MolWt(water_mol)
    product_mass = Descriptors.MolWt(product_mol)

    # Output the results, including the "equation"
    print("Reaction Analysis:")
    print("Reactant + H2O -> Higher Molar Mass Product\n")
    
    print(f"Reactant: {reactant_smiles}")
    print(f"Molar Mass: {reactant_mass:.2f} g/mol\n")
    
    print(f"Water: {water_smiles}")
    print(f"Molar Mass: {water_mass:.2f} g/mol\n")
    
    print(f"Product: {product_smiles}")
    print(f"Molar Mass: {product_mass:.2f} g/mol\n")

    # Final check
    if product_mass > reactant_mass:
        print("As expected, the hydrolysis product has a higher molar mass than the reactant.")
        print("\nThe product with the higher molar mass is:")
        print(product_smiles)
    else:
        print("An error occurred in the analysis, the product is not heavier.")

if __name__ == "__main__":
    solve_reaction()