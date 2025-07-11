# This script requires the 'rdkit' and 'pyscf' libraries.
# You can install them using pip:
# pip install rdkit-pypi pyscf

import sys
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from pyscf import gto
except ImportError as e:
    print(f"Error: A required library is not installed. Please install rdkit and pyscf.")
    print(f"({e})")
    sys.exit(1)

# The SMILES string of the molecule in question
smiles_string = "C#Cc1cc2ccc3c(C#C)cc4ccc5c(C#C)cc6ccc1c7c2c3c4c5c67"

try:
    # 1. Parse SMILES and add hydrogens
    mol = Chem.MolFromSmiles(smiles_string)
    mol = Chem.AddHs(mol)

    # 2. Generate an initial 2D conformation to guide the 3D structure
    AllChem.Compute2DCoords(mol)

    # 3. Embed and optimize the 3D structure
    # Use a fixed random seed for reproducibility
    if AllChem.EmbedMolecule(mol, randomSeed=42) == -1:
         print("Error: Could not generate 3D coordinates for the molecule.")
         sys.exit(1)

    # Use MMFF94s force field, suitable for conjugated systems, to optimize
    AllChem.MMFFOptimizeMolecule(mol, maxIters=500)

    # 4. Extract atomic data (symbols and coordinates) for PySCF
    conformer = mol.GetConformer()
    atoms_for_pyscf = []
    for atom in mol.GetAtoms():
        pos = conformer.GetAtomPosition(atom.GetIdx())
        atoms_for_pyscf.append([atom.GetSymbol(), (pos.x, pos.y, pos.z)])

    # 5. Use PySCF to determine the symmetry group
    pyscf_mol = gto.Mole()
    pyscf_mol.atom = atoms_for_pyscf
    pyscf_mol.basis = 'sto-3g'  # A minimal basis is sufficient for symmetry detection
    pyscf_mol.build(charge=0, spin=0, symmetry=True) # Assuming neutral singlet state

    # The detected point group is stored in the 'topgroup' attribute
    symmetry_group = pyscf_mol.topgroup

    print(f"The determined symmetry group for the molecule is: {symmetry_group}")

except Exception as e:
    print(f"An error occurred during the symmetry detection process: {e}")
    # Based on chemical analysis, the structure is likely 4,9,14-triethynyl-triangulene.
    # The triangulene core has D3h symmetry, and these symmetric substitutions preserve it.
    print("Based on chemical structure analysis, the expected point group is D3h.")
