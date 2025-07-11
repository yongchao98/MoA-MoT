import sys
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from pymatgen.core import Molecule
    from pymatgen.symmetry.analyzer import PointGroupAnalyzer
except ImportError:
    print("This script requires the 'rdkit-pypi' and 'pymatgen' libraries.")
    print("Please install them, for example using pip:")
    print("pip install rdkit-pypi pymatgen")
    sys.exit(1)

def find_molecule_symmetry(smiles_string: str):
    """
    Calculates the point group symmetry of a molecule from its SMILES string.

    Args:
        smiles_string: The SMILES representation of the molecule.
    """
    try:
        # Step 1: Create RDKit molecule object and add hydrogens
        mol = Chem.MolFromSmiles(smiles_string)
        if mol is None:
            print(f"Error: Could not parse the SMILES string: {smiles_string}")
            return

        mol = Chem.AddHs(mol)

        # Step 2: Generate 3D coordinates
        # Embed the molecule in 3D space, using a random seed for reproducibility
        status = AllChem.EmbedMolecule(mol, randomSeed=42)
        if status == -1:
            print("Error: Failed to generate 3D coordinates for the molecule.")
            return

        # Optimize the geometry using a force field for a more accurate structure
        AllChem.MMFFOptimizeMolecule(mol)

        # Step 3: Extract information for Pymatgen
        conformer = mol.GetConformer(0)
        symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
        coords = conformer.GetPositions()

        # Step 4: Create Pymatgen molecule and analyze symmetry
        # A tolerance is used to find symmetry in the numerical coordinates.
        # The default of 0.1 Angstrom is usually sufficient.
        pmg_mol = Molecule(symbols, coords)
        analyzer = PointGroupAnalyzer(pmg_mol, tolerance=0.1)
        point_group = analyzer.get_point_group_symbol()

        # Step 5: Print the result
        print(f"Molecule SMILES: {smiles_string}")
        print(f"The symmetry point group is: {point_group}")

    except Exception as e:
        print(f"An error occurred: {e}")

# The SMILES string for the molecule in question
smiles = "C#Cc1cc2ccc3c(C#C)cc4ccc5c(C#C)cc6ccc1c7c2c3c4c5c67"

# Run the analysis
find_molecule_symmetry(smiles)