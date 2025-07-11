# This script requires the 'rdkit' and 'pymatgen' libraries.
# You can install them using pip:
# pip install rdkit pymatgen
#
# Note: Pymatgen may take a few minutes to install.

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from pymatgen.core.structure import Molecule
    from pymatgen.symmetry.analyzer import PointGroupAnalyzer
except ImportError:
    print("This script requires rdkit and pymatgen.")
    print("Please install them using 'pip install rdkit pymatgen'")
    exit()

def get_molecule_point_group(smiles_string: str):
    """
    Calculates the point group of a molecule from its SMILES string.

    Args:
        smiles_string: The SMILES representation of the molecule.
    """
    # Step 1: Generate a 3D structure from the SMILES string using RDKit
    mol = Chem.MolFromSmiles(smiles_string)
    if mol is None:
        print("Error: Invalid SMILES string.")
        return

    # Add hydrogens to the molecule
    mol = Chem.AddHs(mol)

    # Embed the molecule in 3D space. This generates an initial 3D conformer.
    # The random seed ensures the result is reproducible.
    status = AllChem.EmbedMolecule(mol, randomSeed=42)
    if status == -1:
        print("Error: RDKit could not generate a 3D conformer for this complex molecule.")
        return

    # Optimize the geometry using the UFF (Universal Force Field)
    # This helps find a low-energy structure that reveals the molecule's true symmetry.
    try:
        AllChem.UFFOptimizeMolecule(mol)
    except Exception as e:
        # Optimization can sometimes fail, but we can proceed with the unoptimized conformer.
        print(f"Warning: UFF optimization raised an error: {e}. Using unoptimized structure for analysis.")

    # Step 2: Extract atomic information for Pymatgen
    # Get the optimized 3D coordinates and atomic symbols
    conformer = mol.GetConformer()
    positions = conformer.GetPositions()
    symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]

    # Step 3: Use Pymatgen to determine the point group
    # Create a Pymatgen Molecule object from the symbols and coordinates
    pmg_mol = Molecule(species=symbols, coords=positions)

    # Analyze the molecule's symmetry. The tolerance parameter accounts for
    # small deviations from perfect symmetry in the computationally generated structure.
    try:
        analyzer = PointGroupAnalyzer(pmg_mol, tolerance=0.1)
        point_group = analyzer.get_point_group()
        print(f"The symmetry group of the molecule is: {point_group}")
    except Exception as e:
        print(f"An error occurred during Pymatgen symmetry analysis: {e}")


# The SMILES string for the molecule provided in the task
smiles = "C#Cc1cc2ccc3c(C#C)cc4ccc5c(C#C)cc6ccc1c7c2c3c4c5c67"

# Calculate and print the point group for the molecule
get_molecule_point_group(smiles)
