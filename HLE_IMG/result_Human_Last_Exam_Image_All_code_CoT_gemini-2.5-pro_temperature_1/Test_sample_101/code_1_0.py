# First, ensure you have RDKit installed:
# pip install rdkit

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from collections import defaultdict

def identify_compound_A():
    """
    Identifies and characterizes the final product, Compound A, of the given reaction.
    The reaction is a two-step synthesis:
    1. Imine formation from 3-hydroxy-pyridine-2-carbaldehyde and aniline.
    2. Strecker-type addition of cyanide to the imine.
    The final product is an alpha-aminonitrile.
    """

    # SMILES string for the final product, Compound A: 2-(cyano(phenylamino)methyl)pyridin-3-ol
    compound_A_smiles = "N#CC(Nc1ccccc1)c1ncccc1O"

    # Create an RDKit molecule object from the SMILES string
    mol = Chem.MolFromSmiles(compound_A_smiles)

    if mol is None:
        print("Error: Could not parse the SMILES string for Compound A.")
        return

    # --- Characterize Compound A ---
    
    # 1. Chemical Name
    chemical_name = "2-(cyano(phenylamino)methyl)pyridin-3-ol"
    
    # 2. Molecular Formula
    molecular_formula = rdMolDescriptors.CalcMolFormula(mol)

    # 3. Molecular Weight
    molecular_weight = Descriptors.MolWt(mol)

    # 4. Atom Counts (to fulfill the "output each number" requirement)
    atom_counts = defaultdict(int)
    for atom in mol.GetAtoms():
        atom_counts[atom.GetSymbol()] += 1

    # --- Print the results ---
    print("Identification of Compound A:")
    print("-" * 30)
    print(f"Chemical Name: {chemical_name}")
    print(f"SMILES: {compound_A_smiles}")
    print(f"Molecular Formula: {molecular_formula}")
    print(f"Molecular Weight: {molecular_weight:.2f}")
    print("\nAtom Counts (the numbers in the final 'equation'):")
    for element, count in sorted(atom_counts.items()):
        print(f"Number of {element} atoms: {count}")
    print("-" * 30)


if __name__ == "__main__":
    identify_compound_A()