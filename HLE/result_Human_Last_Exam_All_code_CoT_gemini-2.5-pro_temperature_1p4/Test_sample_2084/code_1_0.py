# First, ensure you have rdkit installed:
# pip install rdkit

from rdkit import Chem
from rdkit.Chem import Descriptors

def solve_reaction_product_mass():
    """
    Identifies the products of a hypothesized ketal hydrolysis and finds the one with the higher molar mass.
    """
    # The provided SMILES string for the reactant is invalid.
    # Based on the fragments (methyl, phenyl, ether, ketal) and reaction conditions (acid hydrolysis),
    # we deduce the reaction is the hydrolysis of the ketal formed from acetophenone and diethylene glycol.
    
    # Products of the hydrolysis reaction:
    # 1. Acetophenone
    # 2. Diethylene glycol

    product1_name = "Acetophenone"
    product1_smiles = "CC(=O)c1ccccc1"
    
    product2_name = "Diethylene glycol"
    product2_smiles = "OCCOCCO"

    # Create RDKit molecule objects from SMILES strings
    mol1 = Chem.MolFromSmiles(product1_smiles)
    mol2 = Chem.MolFromSmiles(product2_smiles)

    # Calculate molar masses
    # The "numbers in the final equation" will be the molar masses of the products.
    mass1 = Descriptors.MolWt(mol1)
    mass2 = Descriptors.MolWt(mol2)

    print("The acid-catalyzed hydrolysis is expected to yield two products:\n")
    
    print(f"Product 1: {product1_name}")
    print(f"SMILES: {product1_smiles}")
    print(f"Molar Mass: {mass1:.2f} g/mol")
    
    print("\n" + "-"*30 + "\n")

    print(f"Product 2: {product2_name}")
    print(f"SMILES: {product2_smiles}")
    print(f"Molar Mass: {mass2:.2f} g/mol")
    
    print("\n" + "="*30 + "\n")

    # Compare molar masses and determine the final answer
    if mass1 > mass2:
        higher_mass_product_name = product1_name
        higher_mass_product_smiles = product1_smiles
    else:
        higher_mass_product_name = product2_name
        higher_mass_product_smiles = product2_smiles
        
    print(f"Comparing the molar masses ({mass1:.2f} g/mol vs {mass2:.2f} g/mol), the product with the higher molar mass is {higher_mass_product_name}.")
    print(f"\nThe SMILES string for the product with the higher molar mass is: {higher_mass_product_smiles}")

solve_reaction_product_mass()