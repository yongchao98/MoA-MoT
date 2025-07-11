import sys
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from pymatgen.core import Molecule
    from pymatgen.symmetry.analyzer import PointGroupAnalyzer
except ImportError:
    print("Please install the required libraries: pip install rdkit pymatgen")
    sys.exit(1)

def get_molecule_point_group(smiles_string: str) -> str:
    """
    Calculates the point group of a molecule from its SMILES string.

    Args:
        smiles_string: The SMILES representation of the molecule.

    Returns:
        The Schoenflies symbol of the point group as a string, or an error message.
    """
    # Create a molecule object from the SMILES string
    mol = Chem.MolFromSmiles(smiles_string)
    if mol is None:
        return "Error: Could not parse the SMILES string."

    # Add hydrogens to the molecule
    mol = Chem.AddHs(mol)

    # Generate a 3D conformation. A random seed ensures reproducibility.
    status = AllChem.EmbedMolecule(mol, randomSeed=42)
    if status == -1:
        return "Error: Could not generate a 3D conformation."

    # Optimize the geometry using the MMFF94 force field for a more realistic structure
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception as e:
        # Proceed with unoptimized structure if optimization fails
        print(f"Warning: Could not optimize molecule, using unoptimized structure. Error: {e}", file=sys.stderr)

    # Extract atomic symbols and 3D coordinates
    try:
        symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
        conformer = mol.GetConformer(0)
        positions = conformer.GetPositions()
    except Exception as e:
        return f"Error: Could not extract atomic data from RDKit molecule. {e}"


    # Create a Pymatgen Molecule object. Pymatgen expects coordinates in Angstroms,
    # which is the default unit for RDKit.
    pmg_mol = Molecule(symbols, positions)

    # Use Pymatgen's PointGroupAnalyzer to find the point group.
    # A slightly larger tolerance is used to account for small deviations
    # from perfect symmetry in the force-field optimized structure.
    try:
        analyzer = PointGroupAnalyzer(pmg_mol, tolerance=0.3)
        point_group = analyzer.get_point_group()
        return point_group
    except Exception as e:
        return f"Error: Pymatgen could not determine the point group. {e}"

# The SMILES string for the molecule provided by the user
smiles = "C#Cc1cc2ccc3c(C#C)cc4ccc5c(C#C)cc6ccc1c7c2c3c4c5c67"

# Get the point group and print the result
point_group_symbol = get_molecule_point_group(smiles)

print(f"The SMILES string is: {smiles}")
print(f"The determined point group of the molecule is: {point_group_symbol}")