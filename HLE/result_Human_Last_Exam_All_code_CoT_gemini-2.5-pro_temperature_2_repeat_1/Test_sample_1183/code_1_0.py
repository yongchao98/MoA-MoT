# To run this script, you need to install RDKit and PySCF.
# You can install them using pip:
# pip install rdkit pyscf

from rdkit import Chem
from rdkit.Chem import AllChem
from pyscf import gto

def get_symmetry_from_smiles(smiles_string):
    """
    Parses a SMILES string, generates a 3D structure,
    and determines its point group symmetry using PySCF.
    
    Args:
        smiles_string (str): The SMILES representation of the molecule.

    Returns:
        str: The identified point group symbol or an error message.
    """
    try:
        # Step 1: Parse SMILES to create a molecule object
        mol = Chem.MolFromSmiles(smiles_string)
        if mol is None:
            return "Error: Invalid SMILES string provided."

        # Step 2: Add hydrogens to complete the structure
        mol = Chem.AddHs(mol)

        # Step 3: Generate an initial 3D conformation
        # Using a fixed random seed for reproducibility
        params = AllChem.ETKDGv3()
        params.randomSeed = 42 
        AllChem.EmbedMolecule(mol, params)

        # Step 4: Optimize the geometry using the UFF force field
        # This gives a reasonable structure for symmetry analysis
        AllChem.UFFOptimizeMolecule(mol)

        # Step 5: Extract atomic symbols and coordinates for PySCF
        atoms_pyscf = []
        conformer = mol.GetConformer()
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            pos = conformer.GetAtomPosition(atom.GetIdx())
            atoms_pyscf.append((symbol, (pos.x, pos.y, pos.z)))

        # Step 6: Use PySCF to build the molecule and find its point group
        pyscf_mol = gto.Mole()
        pyscf_mol.atom = atoms_pyscf
        pyscf_mol.basis = 'sto-3g'  # Basis set choice doesn't affect symmetry finding
        pyscf_mol.build(0, 0) # Build molecule with symmetry detection

        # The determined point group is stored in mol.topgroup (in lowercase)
        point_group_lower = pyscf_mol.topgroup
        
        # Format the point group string to standard notation (e.g., 'c3h' -> 'C3h')
        if len(point_group_lower) > 1 and point_group_lower[0] in 'cdst' and not point_group_lower[1].isdigit():
             formatted_group = point_group_lower.capitalize() # E.g., Cs, Ci
        elif len(point_group_lower) > 2 and point_group_lower[1:-1] == 'inf':
            formatted_group = point_group_lower.capitalize().replace('inf','∞') # E.g., Dinfh -> D∞h
        elif len(point_group_lower) > 1:
            formatted_group = point_group_lower[0].upper() + point_group_lower[1:]
        else: # C1
            formatted_group = point_group_lower.capitalize()

        return formatted_group

    except Exception as e:
        return f"An error occurred during the process: {e}"

# The SMILES string for the molecule in question
smiles = "C#Cc1cc2ccc3c(C#C)cc4ccc5c(C#C)cc6ccc1c7c2c3c4c5c67"

# Perform the analysis and get the symmetry group
symmetry_group = get_symmetry_from_smiles(smiles)

# Print the final result
print(f"The symmetry group for the molecule with SMILES:")
print(smiles)
print(f"is: {symmetry_group}")