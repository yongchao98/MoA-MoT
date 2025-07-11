#
# This script determines the molecular point group from a SMILES string.
# It requires the 'rdkit' and 'ase' libraries. You can install them using pip:
#
# pip install rdkit ase
#

import sys

def find_symmetry_from_smiles(smiles_string):
    """
    Calculates the point group of a molecule from its SMILES string.

    Args:
        smiles_string (str): The SMILES representation of the molecule.

    Returns:
        str: The Schoenflies symbol of the point group, or an error message.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from ase import Atoms
        from ase.symmetry import get_point_group
    except ImportError:
        return "Error: This script requires 'rdkit' and 'ase'. Please run: pip install rdkit ase"

    # Step 1: Create a 3D molecule object from the SMILES string using RDKit.
    mol = Chem.MolFromSmiles(smiles_string)
    if mol is None:
        return f"Error: Could not parse the SMILES string: {smiles_string}"
    
    # Add hydrogen atoms, which are often implicit in SMILES.
    mol = Chem.AddHs(mol)

    # Step 2: Generate a 3D conformation and optimize its geometry.
    # The ETKDGv3 algorithm is good for complex, rigid systems. A random seed ensures reproducibility.
    params = AllChem.ETKDGv3()
    params.randomSeed = 42 # for reproducibility
    if AllChem.EmbedMolecule(mol, params) == -1:
         return "Error: RDKit failed to generate a 3D conformation."

    # Optimize the geometry with the Universal Force Field (UFF).
    try:
        AllChem.UFFOptimizeMolecule(mol)
    except Exception:
        # Optimization might fail on some complex structures, but we can often proceed.
        pass

    # Step 3: Convert the RDKit molecule to an ASE (Atomic Simulation Environment) Atoms object.
    try:
        conformer = mol.GetConformer(0)
        positions = conformer.GetPositions()
        symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
        ase_molecule = Atoms(symbols=symbols, positions=positions)
    except Exception as e:
        return f"Error during conversion to ASE object: {e}"

    # Step 4: Use ASE to determine the point group.
    # A tolerance is used to account for small deviations from perfect symmetry in the optimized structure.
    try:
        point_group = get_point_group(ase_molecule, tolerance=0.1)
        return point_group.sch_symbol
    except Exception as e:
        return f"Error during symmetry analysis: {e}"

if __name__ == '__main__':
    # The SMILES string for the molecule in question.
    smiles = 'C#Cc1cc2ccc3c(C#C)cc4ccc5c(C#C)cc6ccc1c7c2c3c4c5c67'
    
    # Get the symmetry group.
    symmetry_group = find_symmetry_from_smiles(smiles)
    
    # Print the final result.
    print(f"The symmetry group of the molecule is: {symmetry_group}")
