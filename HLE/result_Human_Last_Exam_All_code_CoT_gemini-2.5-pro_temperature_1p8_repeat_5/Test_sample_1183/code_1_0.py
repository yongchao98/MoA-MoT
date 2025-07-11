# First, please make sure you have rdkit and pyscf installed:
# pip install rdkit pyscf

from pyscf import gto
from rdkit import Chem
from rdkit.Chem import AllChem

def get_molecule_symmetry(smiles_string):
    """
    Calculates and prints the point group symmetry of a molecule
    from its SMILES string.
    """
    # Create a molecule object from the SMILES string
    mol_rdkit = Chem.MolFromSmiles(smiles_string)
    if mol_rdkit is None:
        print("Error: Invalid SMILES string.")
        return

    # Add hydrogens to the molecule
    mol_rdkit_h = Chem.AddHs(mol_rdkit)

    # Embed the molecule to generate a 3D conformation
    # A random seed is used for reproducibility
    AllChem.EmbedMolecule(mol_rdkit_h, randomSeed=1)
    
    # Optimize the 3D structure using a force field for a more realistic geometry
    try:
        AllChem.UFFOptimizeMolecule(mol_rdkit_h)
    except Exception as e:
        print(f"Warning: RDKit optimization failed, proceeding with initial structure. Reason: {e}")

    # Extract atom symbols and coordinates
    symbols = [atom.GetSymbol() for atom in mol_rdkit_h.GetAtoms()]
    coordinates_angstrom = mol_rdkit_h.GetConformer().GetPositions()

    # Use PySCF to determine the symmetry
    # Create a PySCF Mole object
    mol_pyscf = gto.Mole()
    mol_pyscf.atom = [[symbols[i], tuple(coordinates_angstrom[i])] for i in range(len(symbols))]
    mol_pyscf.unit = 'Angstrom'
    
    # Enable symmetry detection
    mol_pyscf.symmetry = True
    
    # Build the molecule, which triggers the symmetry analysis
    mol_pyscf.build(verbose=0) # verbose=0 suppresses detailed PySCF output

    # Print the final result
    print(f"SMILES string: {smiles_string}")
    print(f"The point group of the molecule is: {mol_pyscf.topgroup}")

# The SMILES string for the molecule in question
smiles = "C#Cc1cc2ccc3c(C#C)cc4ccc5c(C#C)cc6ccc1c7c2c3c4c5c67"

# Run the function to find and print the symmetry group
get_molecule_symmetry(smiles)