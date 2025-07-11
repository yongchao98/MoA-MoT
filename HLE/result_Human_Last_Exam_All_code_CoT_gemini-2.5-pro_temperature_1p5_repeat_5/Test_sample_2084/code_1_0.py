from rdkit import Chem
from rdkit.Chem import Descriptors

def solve_reaction_product_mass():
    """
    Calculates the molar masses of the potential products of a hydrolysis reaction
    and identifies the product with the higher mass.
    """
    # Step 1: Define the SMILES strings for the plausible products.
    # Product 1: Acetophenone
    acetophenone_smiles = "CC(=O)c1ccccc1"
    # Product 2: 1,2,6-Hexanetriol
    triol_smiles = "OCC(O)CCCCO"

    # Step 2: Create molecule objects from the SMILES strings.
    acetophenone_mol = Chem.MolFromSmiles(acetophenone_smiles)
    triol_mol = Chem.MolFromSmiles(triol_smiles)

    # Step 3: Calculate the molar mass for each product.
    acetophenone_mw = Descriptors.MolWt(acetophenone_mol)
    triol_mw = Descriptors.MolWt(triol_mol)
    
    # Step 4: Display the equation and the molar masses.
    print(f"The hydrolysis reaction yields two products:")
    print(f"1. Acetophenone (SMILES: {acetophenone_smiles})")
    print(f"2. 1,2,6-Hexanetriol (SMILES: {triol_smiles})")
    print("\n--- Molar Mass Calculation ---")
    print(f"Molar Mass of Acetophenone: {acetophenone_mw:.2f} g/mol")
    print(f"Molar Mass of 1,2,6-Hexanetriol: {triol_mw:.2f} g/mol")

    # Step 5: Compare the molar masses and determine the heavier product.
    print("\n--- Conclusion ---")
    if triol_mw > acetophenone_mw:
        heavier_product_name = "1,2,6-Hexanetriol"
        heavier_product_smiles = triol_smiles
    else:
        heavier_product_name = "Acetophenone"
        heavier_product_smiles = acetophenone_smiles
    
    print(f"The product with the higher molar mass is {heavier_product_name}.")
    print(f"Its SMILES string is: {heavier_product_smiles}")
    
    # Final answer in the specified format
    print(f"\n<<<{heavier_product_smiles}>>>")

solve_reaction_product_mass()