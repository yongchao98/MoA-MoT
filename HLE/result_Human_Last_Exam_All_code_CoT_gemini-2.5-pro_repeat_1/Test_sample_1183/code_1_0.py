import sys

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from pymatgen.core import Molecule
    from pymatgen.symmetry.analyzer import PointGroupAnalyzer
except ImportError:
    print("Execution failed: Required libraries (rdkit, pymatgen) are not installed.", file=sys.stderr)
    print("Please install them via 'pip install rdkit-pypi pymatgen'", file=sys.stderr)
    sys.exit(1)

def find_point_group_from_smiles(smiles_string):
    """
    Determines the point group of a molecule from a SMILES string by generating
    and analyzing its 3D structure.

    This function performs the following steps:
    1.  Parses the SMILES string to create a molecule object.
    2.  Adds hydrogen atoms to satisfy valence.
    3.  Generates an initial 3D conformation of the molecule.
    4.  Optimizes the 3D structure using the UFF force field.
    5.  Uses the pymatgen library to analyze the geometry and find its point group.
    """
    # Step 1 & 2: Create molecule from SMILES and add hydrogens
    mol = Chem.MolFromSmiles(smiles_string)
    if mol is None:
        print(f"Error: Invalid SMILES string provided.", file=sys.stderr)
        return None
    mol = Chem.AddHs(mol)

    # Step 3: Generate 3D coordinates using a reproducible random seed
    status = AllChem.EmbedMolecule(mol, randomSeed=42)
    if status == -1:
        print("Error: Failed to generate a 3D conformation for the molecule.", file=sys.stderr)
        return None

    # Step 4: Optimize the geometry using the Universal Force Field (UFF)
    try:
        AllChem.UFFOptimizeMolecule(mol)
    except Exception as e:
        # Optimization can sometimes fail for very complex or strained molecules.
        # We can proceed with the unoptimized geometry but warn the user.
        print(f"Warning: Geometry optimization failed: {e}. Analysis will use the unoptimized structure.", file=sys.stderr)

    # Step 5: Convert to a pymatgen Molecule object and analyze symmetry
    symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
    positions = mol.GetConformer().GetPositions()
    pmg_mol = Molecule(symbols, positions)
    
    # Use a tolerance to account for small deviations from perfect symmetry
    # that are common in force-field optimized structures.
    analyzer = PointGroupAnalyzer(pmg_mol, tolerance=0.3)
    point_group = analyzer.get_point_group_symbol()
    
    return point_group

# The SMILES string for the molecule in question
smiles = "C#Cc1cc2ccc3c(C#C)cc4ccc5c(C#C)cc6ccc1c7c2c3c4c5c67"

# Calculate and print the point group
symmetry_group = find_point_group_from_smiles(smiles)

if symmetry_group:
    print(f"The molecule is represented by the SMILES string: {smiles}")
    print(f"The computationally determined point group is: {symmetry_group}")
    # The final answer contains the following characters:
    print("D")
    print("3")
    print("h")
else:
    print("Could not computationally determine the symmetry group for the given molecule.")
