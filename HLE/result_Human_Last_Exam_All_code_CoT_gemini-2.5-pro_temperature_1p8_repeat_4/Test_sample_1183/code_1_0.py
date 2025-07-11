import sys
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    print("Error: RDKit is not installed. Please install it by running 'pip install rdkit'")
    sys.exit(1)

try:
    from pymatgen.core import Molecule
    from pymatgen.symmetry.analyzer import PointGroupAnalyzer
except ImportError:
    print("Error: Pymatgen is not installed. Please install it by running 'pip install pymatgen'")
    sys.exit(1)

def find_molecule_symmetry(smiles_string: str):
    """
    Calculates and prints the point group symmetry of a molecule from its SMILES string.

    This function performs the following steps:
    1. Generates a 3D molecular structure from the SMILES string using RDKit.
       This includes adding hydrogens, generating a 3D conformer, and optimizing the geometry.
    2. Converts the RDKit molecule object to a Pymatgen Molecule object.
    3. Analyzes the Pymatgen molecule to determine its symmetry point group.
    4. Prints the final result.

    Args:
        smiles_string: The SMILES representation of the molecule.
    """
    print(f"Analyzing SMILES string: {smiles_string}")

    # Step 1: Generate a 3D molecular structure using RDKit
    mol = Chem.MolFromSmiles(smiles_string)
    if mol is None:
        print("Error: Invalid SMILES string provided.")
        return

    mol_with_hs = Chem.AddHs(mol)

    # Embed the molecule to generate a 3D conformer with a reproducible seed
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    try:
        AllChem.EmbedMolecule(mol_with_hs, params)
    except Exception as e:
        print(f"Could not generate 3D coordinates for the molecule. Error: {e}")
        return

    # Optimize the structure to find a low-energy, symmetric conformation
    try:
        AllChem.MMFFOptimizeMolecule(mol_with_hs)
    except Exception:
        # Fallback to UFF if MMFF fails, which can happen for complex structures
        AllChem.UFFOptimizeMolecule(mol_with_hs)

    # Step 2: Convert to Pymatgen object
    symbols = [atom.GetSymbol() for atom in mol_with_hs.GetAtoms()]
    positions = mol_with_hs.GetConformer(0).GetPositions()
    pmg_mol = Molecule(symbols, positions)

    # Step 3: Analyze symmetry with Pymatgen's PointGroupAnalyzer
    # A tolerance is used to check for symmetry operations. A value of 0.3 is robust
    # for geometries optimized with standard force fields.
    try:
        pga = PointGroupAnalyzer(pmg_mol, tolerance=0.3)
        point_group = pga.get_pointgroup()
        print(f"The determined point group is: {point_group}")
    except Exception as e:
        print(f"Could not determine point group. Pymatgen error: {e}")


# The SMILES string for the molecule in question
molecule_smiles = 'C#Cc1cc2ccc3c(C#C)cc4ccc5c(C#C)cc6ccc1c7c2c3c4c5c67'
find_molecule_symmetry(molecule_smiles)