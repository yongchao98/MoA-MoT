# First, ensure you have rdkit installed:
# pip install rdkit

from rdkit import Chem
from rdkit.Chem import Descriptors

def find_heavier_product():
    """
    Identifies the heavier product from the acid-catalyzed hydrolysis of a ketal
    derived from acetophenone and a diol.
    """

    # SMILES strings for the two hydrolysis products
    # Product 1: Acetophenone
    smiles_acetophenone = "CC(=O)c1ccccc1"
    
    # Product 2: 1,4-Butanediol (assumed diol component)
    smiles_diol = "OCCCCO"

    # Create RDKit molecule objects
    mol_acetophenone = Chem.MolFromSmiles(smiles_acetophenone)
    mol_diol = Chem.MolFromSmiles(smiles_diol)

    if mol_acetophenone is None or mol_diol is None:
        print("Error: Could not parse one of the SMILES strings.")
        return

    # Calculate molar masses
    mw_acetophenone = Descriptors.MolWt(mol_acetophenone)
    mw_diol = Descriptors.MolWt(mol_diol)
    
    print(f"Product 1: Acetophenone")
    print(f"SMILES: {smiles_acetophenone}")
    print(f"Molar Mass: {mw_acetophenone:.2f} g/mol")
    print("-" * 20)
    print(f"Product 2: 1,4-Butanediol")
    print(f"SMILES: {smiles_diol}")
    print(f"Molar Mass: {mw_diol:.2f} g/mol")
    print("-" * 20)

    # Compare molar masses and identify the heavier product
    if mw_acetophenone > mw_diol:
        heavier_product_smiles = smiles_acetophenone
        heavier_product_name = "Acetophenone"
    else:
        heavier_product_smiles = smiles_diol
        heavier_product_name = "1,4-Butanediol"

    print(f"The product with the higher molar mass is {heavier_product_name}.")
    print("Its SMILES string is:")
    print(heavier_product_smiles)

if __name__ == "__main__":
    find_heavier_product()
