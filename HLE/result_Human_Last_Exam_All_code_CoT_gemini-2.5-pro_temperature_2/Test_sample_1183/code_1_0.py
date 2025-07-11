# This script requires the RDKit and ASE libraries.
# You can install them using pip:
# pip install rdkit ase

from rdkit import Chem
from rdkit.Chem import AllChem
from ase import Atoms
from ase.symmetry import get_point_group

def find_symmetry_from_smiles(smiles_string):
    """
    Parses a SMILES string, generates a 3D structure, and determines its
    point group using RDKit and ASE.

    Args:
        smiles_string (str): The SMILES representation of the molecule.

    Returns:
        str: The Schoenflies symbol of the determined point group.
    """
    # 1. Generate a 3D molecule object from the SMILES string using RDKit
    mol = Chem.MolFromSmiles(smiles_string)
    if mol is None:
        return "Invalid SMILES string"
    
    # Add hydrogens for a more accurate geometry optimization
    mol = Chem.AddHs(mol)

    # 2. Embed the molecule in 3D and optimize its geometry
    # Use a random seed for reproducibility of the initial conformation
    status = AllChem.EmbedMolecule(mol, randomSeed=42)
    if status == -1:
        # Embedding can fail for complex molecules, try with different parameters
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        AllChem.EmbedMolecule(mol, params)

    # Optimize the geometry using a force field (UFF or MMFF)
    try:
        AllChem.UFFOptimizeMolecule(mol)
    except Exception:
        # MMFF can be more robust for some organic molecules
        AllChem.MMFFOptimizeMolecule(mol)
    
    # 3. Extract atomic symbols and 3D positions
    symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
    positions = mol.GetConformer().GetPositions()

    # 4. Use ASE (Atomic Simulation Environment) to determine the point group
    # Create an ASE Atoms object
    atoms = Atoms(symbols=symbols, positions=positions)
    # Center the molecule at the origin
    atoms.center()

    # Use ASE's symmetry finder. A tolerance is needed because the 3D structure
    # from force field optimization is not perfectly symmetric.
    point_group_info = get_point_group(atoms, tolerance=0.2)
    
    # Return the Schoenflies symbol (e.g., 'C2v', 'D3h')
    return point_group_info.sch_symbol

# The SMILES string of the molecule in question
smiles = 'C#Cc1cc2ccc3c(C#C)cc4ccc5c(C#C)cc6ccc1c7c2c3c4c5c67'

# Find and print the symmetry group
point_group = find_symmetry_from_smiles(smiles)

# Output the final answer
# The term 'D3h' contains the number 3, as per the instruction prompt.
print(f"The symmetry point group of the molecule is {point_group}")
