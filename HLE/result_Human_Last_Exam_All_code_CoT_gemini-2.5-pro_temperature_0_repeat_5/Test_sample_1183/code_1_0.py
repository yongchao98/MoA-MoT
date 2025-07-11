import sys

# This script requires the RDKit and PySCF libraries.
# You can install them using pip:
# pip install rdkit pyscf

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    print("Error: RDKit library not found.")
    print("Please install it using: pip install rdkit")
    sys.exit(1)

try:
    from pyscf import gto
    from pyscf.symm import addons
except ImportError:
    print("Error: PySCF library not found.")
    print("Please install it using: pip install pyscf")
    sys.exit(1)

def find_molecule_point_group(smiles_string: str):
    """
    Calculates the point group of a molecule from its SMILES string.

    Args:
        smiles_string: The SMILES representation of the molecule.

    Returns:
        The Schoenflies symbol of the point group as a string, or None on failure.
    """
    print(f"Analyzing SMILES: {smiles_string}")

    # Step 1: Create molecule from SMILES and add hydrogens for realistic geometry
    mol = Chem.MolFromSmiles(smiles_string)
    if not mol:
        print("Error: Failed to parse SMILES string.")
        return None
    mol = Chem.AddHs(mol)

    # Step 2: Generate a 3D structure and optimize it
    print("Generating 3D structure...")
    # Use a reproducible random seed for embedding
    status = AllChem.EmbedMolecule(mol, randomSeed=42)
    if status == -1:
        print("Warning: RDKit failed to generate an initial 3D conformation. Trying a different algorithm...")
        # Fallback to a more robust embedding algorithm if the default fails
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        status = AllChem.EmbedMolecule(mol, params)
        if status == -1:
            print("Error: Failed to generate 3D structure even with advanced options.")
            return None

    print("Optimizing molecular geometry with UFF force field...")
    try:
        AllChem.UFFOptimizeMolecule(mol)
    except Exception as e:
        print(f"Warning: Geometry optimization failed: {e}. The result may be inaccurate.")

    # Step 3: Prepare atom data for PySCF
    atoms_for_pyscf = []
    conformer = mol.GetConformer()
    for atom in mol.GetAtoms():
        pos = conformer.GetAtomPosition(atom.GetIdx())
        symbol = atom.GetSymbol()
        coords = (pos.x, pos.y, pos.z)
        atoms_for_pyscf.append([symbol, coords])

    # Step 4: Detect symmetry using PySCF
    print("Detecting symmetry group...")
    try:
        # A tolerance is needed because the optimized geometry is not perfectly symmetric.
        # 5e-2 Angstrom is a reasonable tolerance for force field geometries.
        point_group = addons.detect_symm(atoms_for_pyscf, tol=5e-2)
    except Exception as e:
        print(f"An error occurred during symmetry detection: {e}")
        return None

    return point_group

if __name__ == "__main__":
    # The SMILES string of the molecule in question
    smiles = 'C#Cc1cc2ccc3c(C#C)cc4ccc5c(C#C)cc6ccc1c7c2c3c4c5c67'

    # Execute the function and get the result
    symmetry_group = find_molecule_point_group(smiles)

    # Print the final answer
    if symmetry_group:
        print("\n--- RESULT ---")
        print(f"The symmetry group of the molecule is: {symmetry_group}")
        print("--------------")
        # The final answer format requested by the user
        print(f"<<<{symmetry_group}>>>")
    else:
        print("\nCould not determine the symmetry group.")