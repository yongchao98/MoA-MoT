# Import necessary libraries
# Ensure you have rdkit and pyscf installed:
# pip install rdkit pyscf
import sys
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

# Suppress RDKit warnings for cleaner output
from rdkit import rdBase
rdBase.DisableLog('rdApp.warning')

try:
    from pyscf import gto
except ImportError:
    print("PySCF library not found. Please install it using: pip install pyscf")
    sys.exit(1)

def find_molecule_symmetry(smiles_string):
    """
    Calculates the point group symmetry of a molecule from its SMILES string.

    Args:
        smiles_string (str): The SMILES representation of the molecule.

    Returns:
        str: The Schonflies symbol for the point group, or an error message.
    """
    try:
        # 1. Create a molecule object from the SMILES string
        mol_rdkit = Chem.MolFromSmiles(smiles_string)
        if mol_rdkit is None:
            return "Error: Invalid SMILES string."

        # 2. Add hydrogens and generate a 3D conformer
        mol_rdkit_h = Chem.AddHs(mol_rdkit)
        
        # Use a robust embedding method (ETKDGv2)
        params = AllChem.ETKDGv2()
        params.randomSeed = 42 # for reproducibility
        AllChem.EmbedMolecule(mol_rdkit_h, params)
        
        # 3. Optimize the geometry using a force field to get a symmetric structure
        try:
            AllChem.MMFFOptimizeMolecule(mol_rdkit_h, mmffVariant='MMFF94s')
        except Exception as e:
            # Some complex molecules might fail optimization, but we can proceed
            pass
            
        conf = mol_rdkit_h.GetConformer()
        symbols = [atom.GetSymbol() for atom in mol_rdkit_h.GetAtoms()]
        coords_angstrom = conf.GetPositions()

        # 4. Use PySCF to determine the point group
        mol_pyscf = gto.Mole()
        # Create the atom string in PySCF format: [['C', (x, y, z)], ['H', (x, y, z)], ...]
        mol_pyscf.atom = [[symbols[i], tuple(coords_angstrom[i])] for i in range(len(symbols))]
        
        # A basis set is required, but it does not affect symmetry detection
        mol_pyscf.basis = 'sto-3g'
        
        # Explicitly enable symmetry detection
        mol_pyscf.symmetry = True
        
        # Build the molecule, which triggers symmetry analysis
        mol_pyscf.build(verbose=0) # verbose=0 suppresses pyscf output
        
        return mol_pyscf.groupname

    except ImportError:
        return "Error: RDKit library not found. Please install it."
    except Exception as e:
        return f"An error occurred: {e}"

# The SMILES string for the molecule in question
smiles = 'C#Cc1cc2ccc3c(C#C)cc4ccc5c(C#C)cc6ccc1c7c2c3c4c5c67'

# Calculate and print the symmetry group
symmetry_group = find_molecule_symmetry(smiles)
print(f"The point group of the molecule is: {symmetry_group}")