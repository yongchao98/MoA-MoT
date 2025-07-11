# This script requires the rdkit and pymatgen libraries.
# You can install them by running the following command in your terminal:
# pip install rdkit pymatgen

import sys

def find_molecule_symmetry(smiles_string):
    """
    Determines the symmetry point group of a molecule from its SMILES string.

    Args:
        smiles_string (str): The SMILES representation of the molecule.

    Returns:
        str: The Schönflies symbol of the point group, or an error message.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        return "Error: RDKit library not found. Please run 'pip install rdkit'."

    try:
        import pymatgen.core as mg
        from pymatgen.symmetry.analyzer import PointGroupAnalyzer
    except ImportError:
        return "Error: Pymatgen library not found. Please run 'pip install pymatgen'."

    # 1. Create RDKit molecule object from SMILES
    mol = Chem.MolFromSmiles(smiles_string)
    if mol is None:
        return f"Error: Could not parse the SMILES string: {smiles_string}"

    # 2. Add hydrogens and generate a 3D structure
    mol = Chem.AddHs(mol)
    
    # Use a reproducible random seed for embedding
    status = AllChem.EmbedMolecule(mol, randomSeed=42)
    if status == -1:
        return "Error: RDKit failed to generate 3D coordinates for the molecule."
        
    # Optimize the structure using a force field
    try:
        AllChem.UFFOptimizeMolecule(mol)
    except Exception as e:
        # Optimization might fail on complex molecules, but we can often proceed
        print(f"Warning: RDKit UFF optimization raised an error, but attempting symmetry analysis anyway. Details: {e}", file=sys.stderr)


    # 3. Extract atomic data for Pymatgen
    try:
        conformer = mol.GetConformer()
        atom_symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
        positions = conformer.GetPositions()
    except ValueError:
        return "Error: Could not retrieve conformer information from the RDKit molecule."

    # 4. Create Pymatgen molecule object
    pymatgen_mol = mg.Molecule(atom_symbols, positions)

    # 5. Analyze symmetry and get the point group
    # A tolerance value is used to account for small deviations from perfect symmetry
    # in the computationally generated structure.
    pga = PointGroupAnalyzer(pymatgen_mol, tolerance=0.3)
    point_group_symbol = pga.get_point_group_symbol()
    
    return point_group_symbol

if __name__ == '__main__':
    # The molecule provided by the user
    smiles = 'C#Cc1cc2ccc3c(C#C)cc4ccc5c(C#C)cc6ccc1c7c2c3c4c5c67'
    
    # Calculate and print the result
    symmetry_group = find_molecule_symmetry(smiles)
    
    print(f"The SMILES string of the molecule is:")
    print(smiles)
    print("\nThe determined symmetry point group is:")
    print(symmetry_group)
