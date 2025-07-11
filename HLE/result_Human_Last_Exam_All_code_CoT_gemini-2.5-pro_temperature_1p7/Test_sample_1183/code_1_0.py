import sys
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    print("Error: RDKit is not installed. Please install it using 'pip install rdkit-pypi'")
    sys.exit(1)

try:
    from pymatgen.core import Molecule
    from pymatgen.symmetry.analyzer import PointGroupAnalyzer
except ImportError:
    print("Error: Pymatgen is not installed. Please install it using 'pip install pymatgen'")
    sys.exit(1)

def get_molecule_point_group(smiles_string: str) -> str:
    """
    Calculates the point group of a molecule from its SMILES string.

    Args:
        smiles_string: The SMILES representation of the molecule.

    Returns:
        The point group symbol as a string.
    """
    # Create a molecule object from SMILES
    mol = Chem.MolFromSmiles(smiles_string)
    if mol is None:
        return "Invalid SMILES string"

    # Add hydrogens to satisfy valence
    mol = Chem.AddHs(mol)

    # Generate a 3D conformation using a random seed for reproducibility
    try:
        AllChem.EmbedMolecule(mol, randomSeed=42)
        # Optimize the geometry using a force field to find a symmetric conformation
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception as e:
        # Fallback to UFF if MMFF fails or is not available
        try:
            AllChem.UFFOptimizeMolecule(mol)
        except Exception as e_uff:
             return f"3D embedding or optimization failed: {e_uff}"


    # Extract atomic species and 3D coordinates
    conformer = mol.GetConformer()
    species = [atom.GetSymbol() for atom in mol.GetAtoms()]
    coords = conformer.GetPositions()

    # Create a Pymatgen Molecule object
    pymatgen_mol = Molecule(species, coords)

    # Analyze the point group symmetry. A tolerance is used to check for
    # positional equivalency after symmetry operations. 0.3 is a robust value.
    pga = PointGroupAnalyzer(pymatgen_mol, tolerance=0.3)
    point_group = pga.get_pointgroup()
    return point_group

# The SMILES string for the molecule in question
smiles = "C#Cc1cc2ccc3c(C#C)cc4ccc5c(C#C)cc6ccc1c7c2c3c4c5c67"

# Calculate and print the point group
symmetry_group = get_molecule_point_group(smiles)

print(f"The symmetry group for the molecule is: {symmetry_group}")