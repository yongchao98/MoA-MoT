# First, please install rdkit if you haven't yet:
# pip install rdkit

from rdkit import Chem
from rdkit.Chem import Descriptors

def solve_reaction_product():
    """
    Solves for the product with the higher molar mass from a presumed hydrolysis reaction.
    """
    # Step 1: Define the SMILES strings for the likely products based on analysis.
    # The invalid reactant SMILES suggests a ketal made from acetophenone and an ether-diol.
    # The simplest ether-diol is diethylene glycol.
    # Product 1: Acetophenone
    smiles_ketone = "CC(=O)c1ccccc1"
    name_ketone = "Acetophenone"
    
    # Product 2: Diethylene Glycol
    smiles_diol = "OCCOCCO"
    name_diol = "Diethylene Glycol"

    # Step 2: Create RDKit molecule objects from the SMILES strings.
    mol_ketone = Chem.MolFromSmiles(smiles_ketone)
    mol_diol = Chem.MolFromSmiles(smiles_diol)

    if not mol_ketone or not mol_diol:
        print("Error: Could not create molecule from one of the SMILES strings.")
        return

    # Step 3: Calculate the molar mass (Molecular Weight) for each product.
    # This represents the "equation" part of the problem, comparing the two numbers.
    mass_ketone = Descriptors.MolWt(mol_ketone)
    mass_diol = Descriptors.MolWt(mol_diol)
    
    print(f"The reaction produces two main products. We will compare their molar masses.")
    print("-" * 40)
    print(f"Product 1: {name_ketone}")
    print(f"SMILES: {smiles_ketone}")
    print(f"Molar Mass: {mass_ketone:.4f}")
    print("-" * 40)
    print(f"Product 2: {name_diol}")
    print(f"SMILES: {smiles_diol}")
    print(f"Molar Mass: {mass_diol:.4f}")
    print("-" * 40)

    # Step 4: Compare the masses and determine the heavier product.
    if mass_ketone > mass_diol:
        heavier_product_smiles = smiles_ketone
        heavier_product_name = name_ketone
    else:
        heavier_product_smiles = smiles_diol
        heavier_product_name = name_diol
        
    print(f"The product with the higher molar mass is: {heavier_product_name}")
    print(f"\nFinal Answer (SMILES string):")
    print(heavier_product_smiles)

solve_reaction_product()