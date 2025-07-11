# First, ensure you have the RDKit library installed:
# pip install rdkit

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

def get_compound_A_properties():
    """
    Calculates and prints the properties of the final product, Compound A.
    The structure is determined from the reaction mechanism: imine formation
    followed by a Strecker-type cyanide addition.
    """
    # The SMILES string for the final product,
    # 2-(3-hydroxypyridin-2-yl)-2-(phenylamino)acetonitrile, is constructed below.
    compound_A_smiles = 'N#CC(Nc1ccccc1)c1c(O)ccccn1'
    
    # Create a molecule object from the SMILES string
    mol = Chem.MolFromSmiles(compound_A_smiles)
    
    if mol is None:
        print(f"Error: Could not generate molecule from SMILES string: {compound_A_smiles}")
        return
        
    # Calculate key properties of the molecule
    molecular_formula = CalcMolFormula(mol)
    molar_mass = Descriptors.MolWt(mol)
    iupac_name = "2-(3-hydroxypyridin-2-yl)-2-(phenylamino)acetonitrile"

    print("--- Analysis of Final Product: Compound A ---")
    print(f"Reaction Pathway: Imine formation followed by Strecker-type cyanide addition.")
    print(f"\nStructure (SMILES): {compound_A_smiles}")
    print(f"IUPAC Name: {iupac_name}")
    print(f"\nMolecular Formula: {molecular_formula}")
    print(f"Molar Mass: {molar_mass:.4f} g/mol")
    
    # As requested to "output each number in the final equation", this section
    # details the elemental composition (atom counts) of the final product, Compound A.
    print("\nElemental Composition of Compound A:")
    
    # A robust way to count atoms, including implicit hydrogens
    mol_with_hs = Chem.AddHs(mol)
    atom_counts = {}
    for atom in mol_with_hs.GetAtoms():
        symbol = atom.GetSymbol()
        atom_counts[symbol] = atom_counts.get(symbol, 0) + 1
    
    # Print the counts for each element
    for element in sorted(atom_counts.keys()):
        print(f"  - Number of {element} atoms: {atom_counts[element]}")

# Execute the function to display the results
get_compound_A_properties()
