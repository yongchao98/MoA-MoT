import sys
from rdkit import Chem
from rdkit.Chem import AllChem

# Suppress Pymatgen warnings for cleaner output
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="pymatgen")

try:
    from pymatgen.core.structure import Molecule
    from pymatgen.symmetry.analyzer import PointGroupAnalyzer
except ImportError:
    print("Pymatgen library not found. Please install it with 'pip install pymatgen'")
    sys.exit(1)

def get_molecule_point_group(smiles_string: str) -> str:
    """
    Calculates the point group of a molecule from its SMILES string.

    Args:
        smiles_string: The SMILES representation of the molecule.

    Returns:
        The Schoenflies symbol of the molecule's point group.
    """
    # 1. Create a molecule object from the SMILES string
    mol = Chem.MolFromSmiles(smiles_string)
    if mol is None:
        return "Error: Invalid SMILES string."

    # 2. Add hydrogens to get a complete structure
    mol_h = Chem.AddHs(mol)

    # 3. Generate a 3D conformation and optimize it
    # Use a fixed random seed for reproducibility
    AllChem.EmbedMolecule(mol_h, randomSeed=42)
    try:
        AllChem.UFFOptimizeMolecule(mol_h)
    except Exception:
        # Fallback for molecules that UFF might fail on
        pass

    # 4. Extract atomic symbols and coordinates for Pymatgen
    symbols = [atom.GetSymbol() for atom in mol_h.GetAtoms()]
    conformer = mol_h.GetConformer()
    if conformer is None:
        return "Error: Could not generate 3D conformation."
    coords = conformer.GetPositions()

    # 5. Use Pymatgen to find the point group
    # Create a Pymatgen Molecule object
    pymatgen_mol = Molecule(symbols, coords)

    # Use the PointGroupAnalyzer
    # A tolerance is used to account for small deviations from perfect symmetry
    # in the computationally optimized structure.
    pga = PointGroupAnalyzer(pymatgen_mol, tolerance=0.3)
    point_group = pga.get_pointgroup()
    
    return point_group

# The SMILES string for the molecule in question
smiles = "C#Cc1cc2ccc3c(C#C)cc4ccc5c(C#C)cc6ccc1c7c2c3c4c5c67"

# Calculate and print the point group
point_group_symbol = get_molecule_point_group(smiles)
print(f"The symmetry group of the molecule is: {point_group_symbol}")

<<<D3h>>>