#
# This script determines the symmetry point group of a molecule from its SMILES string.
#
# Required libraries: rdkit and pyscf
# You can install them via pip:
# pip install rdkit pyscf
#

import re
import sys

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    print("Error: RDKit library not found.")
    print("Please install it using: 'pip install rdkit'")
    sys.exit(1)

try:
    from pyscf import gto
except ImportError:
    print("Error: PySCF library not found.")
    print("Please install it using: 'pip install pyscf'")
    sys.exit(1)

def find_symmetry_from_smiles(smiles_string):
    """
    Calculates the point group symmetry of a molecule from its SMILES string.

    Args:
        smiles_string (str): The SMILES representation of the molecule.

    Returns:
        str: The Schoenflies symbol for the molecule's point group.
             Returns an error message string if calculation fails.
    """
    # Step 1: Create a molecule object from the SMILES string
    mol = Chem.MolFromSmiles(smiles_string)
    if mol is None:
        return "Error: Invalid SMILES string provided."

    # Step 2: Add hydrogens to complete the valence
    mol = Chem.AddHs(mol)

    # Step 3: Generate a 3D conformation. A random seed is used for reproducibility.
    status = AllChem.EmbedMolecule(mol, randomSeed=1234)
    if status == -1:
        return "Error: Could not generate 3D conformation."
    
    # Step 4: Optimize the geometry using the MMFF94 force field for a realistic structure
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception:
        # Proceeding with unoptimized geometry is possible, but less reliable
        print("Warning: MMFF optimization failed. Symmetry result may be less accurate.")

    # Step 5: Extract atomic symbols and coordinates for PySCF
    conformer = mol.GetConformer()
    symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
    coords = conformer.GetPositions()
    
    # PySCF expects atoms as a list of tuples: [(symbol, (x, y, z)), ...]
    pyscf_atoms = list(zip(symbols, coords))

    # Step 6: Use PySCF to determine the symmetry group
    mol_pyscf = gto.Mole()
    mol_pyscf.atom = pyscf_atoms
    mol_pyscf.basis = 'sto-3g'  # A basis set is required, but not critical for symmetry finding
    mol_pyscf.symmetry = True
    
    # Set a tolerance for symmetry detection. Structures from force fields
    # are not perfectly symmetric, so a larger tolerance than the default is needed.
    # 0.1 Angstrom is a reasonable value.
    mol_pyscf.symmetry_tol = 0.1
    
    # Build the molecule object, which triggers the symmetry analysis
    mol_pyscf.build(0, 0) # Verbosity and output level set to 0

    return mol_pyscf.groupname


# The input SMILES string for the molecule
smiles = 'C#Cc1cc2ccc3c(C#C)cc4ccc5c(C#C)cc6ccc1c7c2c3c4c5c67'

# Per the user request, find all numbers in the SMILES string
numbers_in_smiles = re.findall(r'\d', smiles)
print(f"Numbers found in the SMILES string: {', '.join(numbers_in_smiles)}")

# Determine the symmetry group
symmetry_group = find_symmetry_from_smiles(smiles)
print(f"The determined symmetry group of the molecule is: {symmetry_group}")

# Final answer in the required format
final_answer = symmetry_group
print(f"<<<{final_answer}>>>")