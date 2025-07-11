import sys
import warnings

def find_molecule_point_group(smiles_string: str):
    """
    This function takes a SMILES string, generates a 3D structure,
    and determines its symmetry point group.
    
    The process involves:
    1. Parsing the SMILES string using RDKit.
    2. Generating a 3D conformer and optimizing it with the UFF force field.
    3. Extracting atom symbols and coordinates.
    4. Using PySCF's symmetry finder to determine the point group.
    
    A tolerance is used in the symmetry finding step to account for
    small numerical inaccuracies from the force field optimization.
    """
    # Suppress non-critical warnings for a cleaner output
    warnings.filterwarnings("ignore")
    from rdkit import rdBase
    rdBase.DisableLog('rdApp.warning')

    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        print("Error: RDKit library not found. Please install it using: pip install rdkit-pypi")
        sys.exit(1)
        
    try:
        from pyscf import symm
    except ImportError:
        print("Error: PySCF library not found. Please install it using: pip install pyscf")
        sys.exit(1)

    # Step 1: Parse the SMILES string and add hydrogens for a complete model
    mol = Chem.MolFromSmiles(smiles_string)
    if mol is None:
        print(f"Error: Could not parse the SMILES string: {smiles_string}")
        return
    
    mol = Chem.AddHs(mol)

    # Step 2: Generate 3D coordinates and optimize the geometry
    # We use a random seed for reproducible results
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    if AllChem.EmbedMolecule(mol, params) == -1:
        print("Error: Could not generate 3D conformer for the molecule.")
        return

    # Optimize the geometry using the UFF force field
    try:
        AllChem.UFFOptimizeMolecule(mol, maxIters=500)
    except Exception as e:
        print(f"Warning: Force field optimization failed, but proceeding with unoptimized geometry. Details: {e}")

    # Step 3: Extract atom symbols and coordinates for PySCF
    atoms_pyscf = []
    try:
        conformer = mol.GetConformer()
        for atom in mol.GetAtoms():
            pos = conformer.GetAtomPosition(atom.GetIdx())
            atoms_pyscf.append([atom.GetSymbol(), (pos.x, pos.y, pos.z)])
    except ValueError:
         print("Error: No valid 3D conformer found after embedding attempt.")
         return

    # Step 4: Use PySCF to find the point group
    # A tolerance is needed because UFF-optimized structures are not perfectly symmetric.
    # We use topgroup to find the highest possible symmetry within this tolerance.
    point_group_raw = symm.topgroup(atoms_pyscf, tol=1e-2)
    
    # Format the point group name to standard notation (e.g., D3h)
    if len(point_group_raw) > 1:
        point_group = point_group_raw[0].upper() + point_group_raw[1:].lower()
    else:
        point_group = point_group_raw.capitalize()
        
    print(f"The molecule with SMILES string '{smiles_string}' has been analyzed.")
    print(f"The determined symmetry point group is: {point_group}")
    
    return point_group

if __name__ == '__main__':
    # The SMILES string provided by the user
    smiles = "C#Cc1cc2ccc3c(C#C)cc4ccc5c(C#C)cc6ccc1c7c2c3c4c5c67"
    find_molecule_point_group(smiles)