# First, ensure you have the RDKit library installed.
# If not, you can install it via pip:
# pip install rdkit-pypi

import re
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula
except ImportError:
    print("Error: The 'rdkit' library is required. Please install it using 'pip install rdkit-pypi'")
    exit()

def identify_compound_A():
    """
    Identifies and characterizes Compound A from the reaction sequence.
    The reaction is a Strecker synthesis on an imine.
    Compound A is 2-((cyano)(phenylamino)methyl)pyridin-3-ol.
    """
    
    # The structure of Compound A is represented by the following SMILES string.
    # Oc1c(C(C#N)Nc2ccccc2)nccc1
    # O=hydroxyl, c1...c1=aromatic ring, C=carbon, #N=nitrile, N=amine, n=nitrogen in ring
    compound_a_smiles = "Oc1c(C(C#N)Nc2ccccc2)nccc1"
    
    # Create an RDKit molecule object from the SMILES string
    molecule = Chem.MolFromSmiles(compound_a_smiles)
    
    if molecule is None:
        print("Error: Invalid SMILES string. Could not create molecule.")
        return

    print("Compound A has been identified as: 2-((cyano)(phenylamino)methyl)pyridin-3-ol")
    print("-" * 50)
    
    # Calculate and print the molecular formula
    molecular_formula = CalcMolFormula(molecule)
    print(f"Molecular Formula: {molecular_formula}")
    
    print("\nThe numbers in the final molecular formula are:")
    # Use regular expressions to parse the formula and print the count of each atom
    atom_counts = re.findall(r'([A-Z][a-z]*)(\d*)', molecular_formula)
    
    for element, count in atom_counts:
        # If count is an empty string, it means the count is 1
        atom_count = int(count) if count else 1
        print(f"Number of {element} atoms: {atom_count}")

    # Calculate and print the molecular weight
    molecular_weight = Descriptors.MolWt(molecule)
    print(f"\nMolecular Weight: {molecular_weight:.4f} g/mol")

if __name__ == '__main__':
    identify_compound_A()
