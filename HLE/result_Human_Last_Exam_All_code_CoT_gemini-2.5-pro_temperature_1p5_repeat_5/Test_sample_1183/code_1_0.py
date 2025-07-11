try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from pymatgen.core import Molecule
    from pymatgen.symmetry.analyzer import PointGroupAnalyzer
except ImportError:
    print("This script requires RDKit and Pymatgen.")
    print("Please install them using: pip install rdkit pymatgen")
    exit()

def find_molecule_symmetry(smiles_string: str):
    """
    Calculates the point group symmetry of a molecule from its SMILES string.

    Args:
        smiles_string: The SMILES representation of the molecule.

    Returns:
        The Schoenflies symbol of the point group, or None on failure.
    """
    # --- Step 1 & 2: Parse SMILES and generate 3D structure ---
    mol = Chem.MolFromSmiles(smiles_string)

    if mol is None:
        print(f"Error: Could not parse the SMILES string: '{smiles_string}'")
        return None

    # Add hydrogens to get the complete chemical structure
    mol = Chem.AddHs(mol)

    # Generate an initial 3D conformation using a reproducible method
    params = AllChem.ETKDGv3()
    params.randomSeed = 0xf00d # for reproducibility
    if AllChem.EmbedMolecule(mol, params) == -1:
        print("Error: Could not generate 3D conformation.")
        return None

    # --- Step 3: Optimize the molecular geometry ---
    try:
        # Use MMFF94, a good general-purpose force field
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    except Exception:
        # Fallback to UFF if MMFF fails
        AllChem.UFFOptimizeMolecule(mol, maxIters=500)

    # --- Step 4: Extract structural data for Pymatgen ---
    try:
        symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
        # Get positions from the optimized conformer
        positions = mol.GetConformer().GetPositions()
        
        # Create a Pymatgen Molecule object
        pmg_mol = Molecule(symbols, positions)

        # --- Step 5: Analyze symmetry and determine point group ---
        # The tolerance parameter handles small deviations from perfect symmetry.
        # 0.3 is a reasonable value for force-field optimized structures.
        analyzer = PointGroupAnalyzer(pmg_mol, tolerance=0.3)
        point_group = analyzer.get_pointgroup_symbol()
        
        return point_group

    except Exception as e:
        print(f"An error occurred during symmetry analysis: {e}")
        return None

# The SMILES string for the molecule in question
smiles = 'C#Cc1cc2ccc3c(C#C)cc4ccc5c(C#C)cc6ccc1c7c2c3c4c5c67'

# Run the analysis
symmetry_group = find_molecule_symmetry(smiles)

# Print the final result
if symmetry_group:
    print(f"The SMILES string is: {smiles}")
    print(f"The determined symmetry point group for the molecule is: {symmetry_group}")
else:
    print("Could not determine the symmetry group.")

<<<D3h>>>