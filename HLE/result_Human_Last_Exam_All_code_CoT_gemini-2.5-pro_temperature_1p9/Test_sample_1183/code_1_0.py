try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    import pymatgen.core as mg
    from pymatgen.symmetry.analyzer import PointGroupAnalyzer
except ImportError:
    print("This script requires the 'rdkit-pypi' and 'pymatgen' packages.")
    print("Please install them, for example using pip: 'pip install rdkit-pypi pymatgen'")
    exit()

def get_molecule_point_group(smiles_string: str) -> str:
    """
    Calculates the point group of a molecule from its SMILES string.

    Args:
        smiles_string: The SMILES representation of the molecule.

    Returns:
        The Schönflies symbol of the point group.
    """
    # 1. Parse SMILES and create a molecule object
    mol = Chem.MolFromSmiles(smiles_string)
    if mol is None:
        return "Error: Invalid SMILES string."

    # 2. Add hydrogens and generate a 3D conformation
    mol = Chem.AddHs(mol)
    # Use a fixed random seed for reproducibility of the initial 3D structure
    embed_result = AllChem.EmbedMolecule(mol, randomSeed=42)
    if embed_result == -1:
        return "Error: Failed to generate 3D conformation."

    # 3. Optimize the geometry using MMFF94 force field
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception as e:
        # Some complex structures might have issues with MMFF, try UFF as a fallback
        try:
            AllChem.UFFOptimizeMolecule(mol)
        except Exception as e_uff:
             return f"Error: Geometry optimization failed. Details: {e_uff}"

    # 4. Extract atomic data for Pymatgen
    try:
        conformer = mol.GetConformer()
        species = [atom.GetSymbol() for atom in mol.GetAtoms()]
        coords = conformer.GetPositions()

        # 5. Create a Pymatgen molecule object
        pmg_mol = mg.Molecule(species, coords)

        # 6. Analyze the symmetry and get the point group
        # The tolerance is adjusted to robustly find symmetry in computed structures
        pga = PointGroupAnalyzer(pmg_mol, tol=0.3)
        point_group = pga.get_point_group_symbol()
        return point_group
    except Exception as e:
        return f"An error occurred during symmetry analysis: {e}"


# The SMILES string for the molecule in question
smiles = "C#Cc1cc2ccc3c(C#C)cc4ccc5c(C#C)cc6ccc1c7c2c3c4c5c67"

# Calculate and print the symmetry group
point_group_symbol = get_molecule_point_group(smiles)

print(f"The SMILES string is: {smiles}")
print(f"The symmetry group of the molecule is: {point_group_symbol}")

<<<D3h>>>